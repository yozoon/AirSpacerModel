#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include "CSVReader.hpp"
#include "CubicSplineInterpolation.hpp"
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
    NumericType fillingRatioThreshold = 0.9;

    // Now assign the items
    Utils::AssignItems(
        map, Utils::Item{"stickingProbability", stickingProbability},
        Utils::Item{"trenchBottomWidth", trenchBottomWidth},
        Utils::Item{"trenchDepth", trenchDepth},
        Utils::Item{"topWidthLimit", topWidthLimit},
        Utils::Item{"fillingRatioThreshold", fillingRatioThreshold,
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

      // int numLeftFeaturesAboveZero =
      //     static_cast<int>(std::ceil(numFeaturesLeft * aboveZeroRatio));

      // Percentage of the trench below the zero line that is filled with the
      // non-conformal coating. Count values smaller than 2% of top opening as
      // closed. Since the trench is symmetrical we only look at the features on
      // the right.
      int closedCount = std::count_if(
          std::next(d.begin(), inputDim + numRightFeaturesAboveZero),
          std::next(d.begin(), inputDim + numFeaturesRight),
          [](NumericType p) { return p <= 0.02; });

      int belowZeroCount = numFeaturesRight - numRightFeaturesAboveZero;

      // NumericType cummulativeOpening = 0.;
      // for (unsigned i = numRightFeaturesAboveZero; i < numFeaturesRight; ++i)
      // {
      //   auto diff =
      //       d[inputDim + i] -
      //       d[inputDim + numFeaturesRight + std::min(numFeaturesLeft - 1U,
      //       i)];
      //   cummulativeOpening += std::max(0., diff);
      // }

      // auto filledAreaFraction = 1.0 - (1.0 + 2.0 * aspectRatio * tan) /
      //                                     (1.0 + aspectRatio * tan) *
      //                                     cummulativeOpening /
      //                                     belowZeroCount;

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
          auto [value, _] = valOpt.value();
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
  {
    // Copy the values corresponding to the maximum recorded process time
    std::vector<std::vector<NumericType>> maxTimeBoundary;
    std::copy_if(dataSlice->begin(), dataSlice->end(),
                 std::back_inserter(maxTimeBoundary),
                 [](auto &row) { return row.at(1) == 1.0; });

    // Sort those values based on their taper angle
    std::sort(maxTimeBoundary.begin(), maxTimeBoundary.end(),
              [](auto &a, auto &b) { return a[0] < b[0]; });

    // Extract angles and filling ratios into separate vectors
    std::vector<NumericType> angles;
    std::vector<std::vector<NumericType>> fillingRatios;
    std::transform(maxTimeBoundary.begin(), maxTimeBoundary.end(),
                   std::back_inserter(angles), [](auto &a) { return a[0]; });
    std::transform(maxTimeBoundary.begin(), maxTimeBoundary.end(),
                   std::back_inserter(fillingRatios),
                   [](auto &a) { return std::vector{a[2]}; });

    // Instantiate the cubic spline interpolation
    CubicSplineInterpolation<NumericType> interp(
        angles, fillingRatios, SplineBoundaryConditionType::NOT_A_KNOT);

    // Do binary search to find the taper angle where the filling ratio has a
    // certain threshold value.
    NumericType eps = 1e-5;
    NumericType minAngle = angles.front();
    NumericType maxAngle = angles.back();
    unsigned maxIter = 100;
    unsigned i = 0;
    for (; i < maxIter; ++i) { // Do a maximum of `maxIter` iterations
      NumericType centerAngle = (maxAngle + minAngle) / 2;
      if (interp(centerAngle)[0] <= params.fillingRatioThreshold) {
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
          .addWarning(std::string("No convergence in ") +
                      std::to_string(maxIter) + std::string(" iterations"))
          .print();
    }

    // Use the midpoint as the threshold andle
    thresholdAngle = (maxAngle + minAngle) / 2;

    // Calculate the resulting top width
    topWidth =
        params.trenchBottomWidth +
        2.0 * params.trenchDepth * std::tan(Utils::deg2rad(thresholdAngle));
  }

  // Print the result
  bool success = topWidth <= params.topWidthLimit;
  std::cout << (success ? "SUCCESS: " : "")
            << "The minimum taper angle required to achieve a filling ratio of "
               "at least "
            << 100.0 * params.fillingRatioThreshold << "% is " << thresholdAngle
            << "°.\nThis leads to a trench top width of " << topWidth
            << " which " << (success ? "does" : "does not")
            << " meet the specified top width limit.\n";
  return EXIT_SUCCESS;
}