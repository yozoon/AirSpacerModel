#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#ifdef WITH_VIENNAPS
#include "TrenchGeometry.hpp"
#include <SimpleDeposition.hpp>
#include <lsWriteVisualizationMesh.hpp>
#include <psProcess.hpp>
#include <psSmartPointer.hpp>
#endif

#include "CSVReader.hpp"
#include "SplineGridInterpolation.hpp"
#include "Utils.hpp"

template <typename NumericType> class Parameters {
public:
  const NumericType stickingProbability;
  const NumericType trenchBottomWidth;
  const NumericType trenchDepth;
  const NumericType topWidthLimit;
  const NumericType fillingRatioThreshold;

  // Factory pattern to construct an instance of the class from a map
  static Parameters
  fromMap(const std::unordered_map<std::string, std::string> &map) {
    // Local variable instances
    NumericType stickingProbability = 0.01;
    NumericType trenchBottomWidth = 5;
    NumericType trenchDepth = 50;
    NumericType topWidthLimit = 10;
    NumericType fillingRatioThreshold = 0.99;

    // Now assign the items
    Utils::AssignItems(map,
                       Utils::Item{"stickingProbability", stickingProbability,
                                   &Utils::toUnitRange<NumericType>},
                       Utils::Item{"trenchBottomWidth", trenchBottomWidth,
                                   &Utils::toStrictlyPositive<NumericType>},
                       Utils::Item{"trenchDepth", trenchDepth,
                                   &Utils::toStrictlyPositive<NumericType>},
                       Utils::Item{"topWidthLimit", topWidthLimit,
                                   &Utils::toStrictlyPositive<NumericType>},
                       Utils::Item{"fillingRatioThreshold",
                                   fillingRatioThreshold,
                                   &Utils::toUnitRange<NumericType>});

    return Parameters(stickingProbability, trenchBottomWidth, trenchDepth,
                      topWidthLimit, fillingRatioThreshold);
  }

  void print() const {
    std::cout << "Parameters: " << std::endl;
    std::cout << "- stickingProbability: " << stickingProbability << '\n';
    std::cout << "- trenchBottomWidth: " << trenchBottomWidth << '\n';
    std::cout << "- trenchDepth: " << trenchDepth << '\n';
    std::cout << "- topWidthLimit: " << topWidthLimit << '\n';
    std::cout << "- fillingRatioThreshold: " << fillingRatioThreshold << '\n';
  }

private:
  // Private constructor, so that the user is forced to use the factory
  Parameters(NumericType passedStickingProbability,
             NumericType passedTrenchBottomWidth, NumericType passedTrenchDepth,
             NumericType passedTopWidthLimit,
             NumericType passedFillingRatioThreshold)
      : stickingProbability(passedStickingProbability),
        trenchBottomWidth(passedTrenchBottomWidth),
        trenchDepth(passedTrenchDepth), topWidthLimit(passedTopWidthLimit),
        fillingRatioThreshold(passedFillingRatioThreshold) {}
};

int main(const int argc, const char *const *const argv) {
  // Input: sticking probability, initial_bottom_width, desired depth, maximum
  // top width Output: is there a solution at all. If there is a solution:
  // minimum taper angle that ensures trench filling.
  using NumericType = double;

  if (argc != 3) {
    std::cout << "Usage: " << argv[0] << " <dataset csv file> <config file>\n";
    return EXIT_SUCCESS;
  }

  std::string dataFilename = "trenchfill.csv";
  if (argc > 1)
    dataFilename = argv[1];

  // Load the data captured from previous simulations
  CSVReader<NumericType> reader;
  reader.setFilename(dataFilename);
  auto data = reader.getData();

  std::unordered_map<std::string, std::string> config;
  if (argc > 2) {
    config = Utils::readConfigFile(argv[2]);
    if (config.empty()) {
      lsMessage::getInstance().addError("Empty config provided").print();
      return EXIT_FAILURE;
    }
  }
  auto params = Parameters<NumericType>::fromMap(config);
  params.print();

  // Calculate the aspect ratio
  NumericType aspectRatio = params.trenchDepth / params.trenchBottomWidth;

  // The input dimension is provided by the csv file named parameter
  // `InputDim` In our case it is 4 since we use aspectRatio, leftTaperAngle,
  // stickintProbability and timeStep as input parameters.
  int inputDim{};
  auto namedParams = reader.getNamedParameters();
  if (auto id = namedParams.find("InputDimension"); id != namedParams.end()) {
    inputDim = static_cast<int>(std::round(id->second));
  } else {
    lsMessage::getInstance()
        .addError("`InputDimension` not found in provided data CSV file.")
        .print();
  }

  // The dimension of the data that is to be interpolated. In our case this
  // are the extracted dimensions at different timesteps (the timesteps are
  // also included in the data itself)
  int numFeatures = data->at(0).size() - inputDim;

  // First convert the input data into dataset containing filling ratio
  auto fillingRatioData =
      lsSmartPointer<std::vector<std::vector<NumericType>>>::New();
  {
    unsigned numFeaturesRight =
        static_cast<unsigned>(std::ceil(1.0 * numFeatures / 2));
    // unsigned numFeaturesLeft = numFeatures - numFeaturesRight;

    fillingRatioData->reserve(data->size());
    for (auto &d : *data) {
      NumericType aspectRatio = d[0];
      NumericType taperAngle = d[1];
      NumericType stickingProbability = d[2];
      NumericType time = d[3];

      auto tan = std::tan(Utils::deg2rad(taperAngle));

      // The percentage of features that are above the initial trench top
      // surface
      NumericType aboveZeroRatio =
          (1.0 + 2.0 * aspectRatio * tan) /
          (1.0 + 2.0 * aspectRatio * tan + aspectRatio);

      int numRightFeaturesAboveZero =
          static_cast<int>(std::ceil(numFeaturesRight * aboveZeroRatio));

      // Percentage of the trench below the zero line that is filled with the
      // non-conformal coating. Count values smaller than 2% of top opening as
      // closed. Since the trench is symmetrical we only look at the features on
      // the right.
      const NumericType sampleThreshold = 0.0;
      int closedCount = std::count_if(
          std::next(d.begin(), inputDim + numRightFeaturesAboveZero),
          std::next(d.begin(), inputDim + numFeaturesRight),
          [=](NumericType p) { return p <= sampleThreshold; });

      int belowZeroCount = numFeaturesRight - numRightFeaturesAboveZero;

      NumericType fillingRatio = 1.0 * closedCount / belowZeroCount;
      fillingRatioData->emplace_back(std::vector<NumericType>{
          aspectRatio, taperAngle, stickingProbability, time, fillingRatio});
    }
  }

  // Second: extract slice of data at the provided sticking probability and
  // aspect ratio
  auto dataSlice = lsSmartPointer<std::vector<std::vector<NumericType>>>::New();
  bool hasSolution = false;
  {
    SplineGridInterpolation<NumericType> sgi;
    sgi.setDataDimensions(inputDim, 1);
    sgi.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    sgi.setData(fillingRatioData);
    sgi.initialize();
    auto uniqueValues = sgi.getUniqueValues();

    for (auto taperAngle : uniqueValues[1]) {
      for (auto time : uniqueValues[3]) {
        std::vector<NumericType> loc = {aspectRatio, taperAngle,
                                        params.stickingProbability, time};

        // auto grad = sgi.gradient(loc);
        auto valOpt = sgi.estimate(loc);
        if (valOpt) {
          auto [value, isInside] = valOpt.value();
          if (!isInside) {
            lsMessage::getInstance()
                .addError("The provided parameters are not within the "
                          "boundaries of the dataset.")
                .print();
          }
          dataSlice->emplace_back(
              std::vector<NumericType>{taperAngle, time, value[0]});
          if (value[0] >= params.fillingRatioThreshold)
            hasSolution = true;
        }
      }
    }
  }

  if (!hasSolution) {
    std::cout
        << "There exists no solution in the range spanned by the dataset.\n";
    return EXIT_SUCCESS;
  }

  // Now check if the minimum taper angle fulfills the requirement for the
  // maximum top width
  NumericType thresholdAngle;
  NumericType topWidth;
  NumericType normalizedProcessTime;
  {
    SplineGridInterpolation<NumericType> sgi;
    sgi.setDataDimensions(dataSlice->at(0).size() - 1, 1);
    sgi.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    sgi.setData(dataSlice);
    sgi.initialize();
    auto uniqueValues = sgi.getUniqueValues();
    auto angles = uniqueValues.at(0);

    // Do binary search to find the taper angle where the filling ratio has a
    // certain threshold value.
    NumericType eps = 1e-5;
    NumericType minAngle = *angles.begin();
    NumericType maxAngle = *angles.rbegin();
    unsigned maxIter = 100;
    unsigned i = 0;
    for (; i < maxIter; ++i) { // Do a maximum of `maxIter` iterations
      NumericType centerAngle = (maxAngle + minAngle) / 2;
      std::vector<NumericType> loc = {centerAngle, 1.0};
      auto [frEstimate, _] = sgi.estimate(loc).value();
      if (frEstimate.at(0) <= params.fillingRatioThreshold) {
        minAngle = centerAngle;
      } else {
        maxAngle = centerAngle;
      }
      if (maxAngle - minAngle < eps) {
        lsMessage::getInstance()
            .addDebug(std::string("iterations until convergence: ") +
                      std::to_string(i))
            .print();
        break;
      }
    }
    if (i >= maxIter - 1) {
      lsMessage::getInstance()
          .addError(std::string("No convergence in ") +
                    std::to_string(maxIter) + std::string(" iterations"))
          .print();
    }

    // Calculate the resulting top width
    topWidth = params.trenchBottomWidth +
               2.0 * params.trenchDepth * std::tan(Utils::deg2rad(maxAngle));

    // Extend the top width by a few percent,
    NumericType leeway = 0.05;
    thresholdAngle =
        180.0 *
        std::atan((1.0 + leeway) * (topWidth - params.trenchBottomWidth) /
                  (2.0 * params.trenchDepth)) /
        Utils::PI<NumericType>;
    topWidth = (1.0 + leeway) * topWidth;

    // Now check the minimum time required to reach this value
    NumericType stepsize = 1.0 / maxIter;
    for (unsigned i = 0; i < maxIter; ++i) {
      normalizedProcessTime = 1.0 - i * stepsize;
      std::vector<NumericType> loc = {thresholdAngle, normalizedProcessTime};
      auto [frEstimate, _] = sgi.estimate(loc).value();
      if (frEstimate.at(0) < params.fillingRatioThreshold)
        break;
    }
  }

  // Print the result
  bool success = topWidth <= params.topWidthLimit;
  std::cout
      << (success ? "\nSUCCESS!\n" : "\nNO SOLUTION FOUND!\n")
      << "The minimum taper angle required to achieve a filling ratio\nof "
         "at least "
      << 100.0 * params.fillingRatioThreshold << "% is " << thresholdAngle
      << "°.\nThis leads to a trench top width of " << topWidth << " which "
      << (success ? "fulfills" : "does not fulfill")
      << "\nthe specified top width constraint.\nAssuming a rate of 1, the "
         "minimum required process time is "
      << topWidth * normalizedProcessTime << "\n\n";

#ifdef WITH_VIENNAPS
  static constexpr int D = 2;
  if (success) {
    std::cout << "Should I run a physical simulation to verify the "
                 "solution? [y/n]";
    char c;
    std::cin >> c;
    if (c == 'y' || c == 'Y') {
      std::cout << "Running simulation...\n";
      const NumericType gridDelta = std::max(params.trenchDepth / 100, 0.2);

      // The horizontal extent is defined by the trench top width, or the bottom
      // width, whichever is larger. Extend it by 20%.
      NumericType horizontalExtent = topWidth * 1.2;

      // Generate the grid
      std::array<NumericType, 3> origin{0.};
      auto grid = createGrid<NumericType, D>(
          origin, gridDelta, horizontalExtent,
          10.0 /* initial vertical extent */, false /* No periodic BC */);

      // Create the inside of the trench that will be removed from the plane
      auto cutout = createTrenchStamp<NumericType, D>(
          grid, origin, params.trenchDepth, topWidth, thresholdAngle,
          thresholdAngle, true);

      auto trench = createPlane<NumericType, D>(grid, origin);
      lsBooleanOperation<NumericType, D>(
          trench, cutout, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();

      // Since we assume a top rate of 1, the maximum time it takes for the
      // trench to close is equal to the trench top width (if sticking
      // probability 1, otherwise the time until pinchoff should be even less)
      // Additionally ormalize the time scale to the sticking probability, so
      // that we get the same top layer thickness for different sticking
      // probabilities.
      NumericType processDuration =
          normalizedProcessTime * topWidth / params.stickingProbability;

      auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
      geometry->insertNextLevelSet(trench);

      // copy top layer to capture deposition
      auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(trench);
      geometry->insertNextLevelSet(depoLayer);

      // Instantiate the process
      auto processModel =
          SimpleDeposition<NumericType, D>(params.stickingProbability,
                                           1.0 /* particle source power */)
              .getProcessModel();

      psProcess<NumericType, D> process;
      process.setDomain(geometry);
      process.setProcessModel(processModel);
      process.setNumberOfRaysPerPoint(2000);
      process.setProcessDuration(processDuration);

      // Run the process
      process.apply();

      std::cout << "Saving result as 'FilledTrench_volume.vtu'\n";
      lsWriteVisualizationMesh<NumericType, D> visMesh;
      for (auto ls : *geometry->getLevelSets()) {
        visMesh.insertNextLevelSet(ls);
      }
      visMesh.setFileName("FilledTrench");
      visMesh.apply();
    }
  }
#endif

  return EXIT_SUCCESS;
}