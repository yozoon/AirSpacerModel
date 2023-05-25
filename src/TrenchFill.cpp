#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#ifdef WITH_VIENNAPS
#include <SimpleDeposition.hpp>
#include <lsWriteVisualizationMesh.hpp>
#include <psProcess.hpp>
#include <psSmartPointer.hpp>

#include "TrenchGeometry.hpp"
#endif

#include "span.hpp"

#include "CSVReader.hpp"
#include "SplineGridInterpolation.hpp"
#include "Utils.hpp"

#define SAMPLE_SLICE 1

template <typename NumericType> class Parameters {
public:
  const NumericType stickingProbability;
  const NumericType trenchBottomWidth;
  const NumericType trenchDepth;
  const NumericType topWidthLimit;
  const NumericType taperAngleTolerance;

  // Factory pattern to construct an instance of the class from a map
  static Parameters
  fromMap(const std::unordered_map<std::string, std::string> &map) {
    // Local variable instances
    NumericType stickingProbability = 0.01;
    NumericType trenchBottomWidth = 5;
    NumericType trenchDepth = 50;
    NumericType topWidthLimit = 10;

    // How much the angle is allowed to increase while following the isoline
    NumericType taperAngleTolerance = 2.0;

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
                       Utils::Item{"taperAngleTolerance", taperAngleTolerance});

    return Parameters(stickingProbability, trenchBottomWidth, trenchDepth,
                      topWidthLimit, taperAngleTolerance);
  }

  void print() const {
    std::cout << "Parameters: " << std::endl;
    std::cout << "- stickingProbability: " << stickingProbability << '\n';
    std::cout << "- trenchBottomWidth: " << trenchBottomWidth << '\n';
    std::cout << "- trenchDepth: " << trenchDepth << '\n';
    std::cout << "- topWidthLimit: " << topWidthLimit << '\n';
    std::cout << "- taperAngleTolerance: " << taperAngleTolerance << '\n';
  }

private:
  // Private constructor, so that the user is forced to use the factory
  Parameters(NumericType passedStickingProbability,
             NumericType passedTrenchBottomWidth, NumericType passedTrenchDepth,
             NumericType passedTopWidthLimit,
             NumericType passedTaperAngleTolerance)
      : stickingProbability(passedStickingProbability),
        trenchBottomWidth(passedTrenchBottomWidth),
        trenchDepth(passedTrenchDepth), topWidthLimit(passedTopWidthLimit),
        taperAngleTolerance(passedTaperAngleTolerance) {}
};

template <typename NumericType>
NumericType toFillingRatio(const NumericType aspectRatio,
                           const NumericType taperAngle,
                           const nonstd::span<const NumericType> rightFeatures,
                           const nonstd::span<const NumericType> leftFeatures) {

  auto tan = std::tan(Utils::deg2rad(taperAngle));

  // The percentage of features that are above the initial trench top
  // surface
  NumericType aboveZeroRatio = (1.0 + 2.0 * aspectRatio * tan) /
                               (1.0 + 2.0 * aspectRatio * tan + aspectRatio);

  int aboveSurfaceCount =
      static_cast<int>(std::ceil(rightFeatures.size() * aboveZeroRatio));

  int closedCount = 0;
  for (int i = 0; i < rightFeatures.size(); ++i) {
    if (rightFeatures[i] - leftFeatures[std::min(
                               static_cast<int>(leftFeatures.size()) - 1, i)] <=
        0) {
      ++closedCount;
    }
  }

  // Percentage of the trench below the zero line that is filled with the
  // non-conformal coating. Count values smaller than 2% of top opening as
  // closed. Since the trench is symmetrical we only look at the features on
  // the right.
  // int closedCount = std::count_if(rightFeatures.begin(), rightFeatures.end(),
  //[=](NumericType p) {return p == 0.0; });
  int closedAboveSurfaceCount =
      std::count_if(rightFeatures.begin(),
                    std::next(rightFeatures.begin(), aboveSurfaceCount),
                    [=](NumericType p) { return p == 0.0; });

  std::cout << rightFeatures.size() << ", " << closedCount << ", "
            << closedAboveSurfaceCount << '\n';

  int closedBelowSurfaceCount = closedCount - closedAboveSurfaceCount;
  int belowSurfaceCount = rightFeatures.size() - aboveSurfaceCount;

  NumericType fillingRatio = 1.0 * closedBelowSurfaceCount / belowSurfaceCount;
  // smoothly extend the value above 1.0, such
  // that gradient based approaches can also be used
  if (fillingRatio >= 1.0)
    fillingRatio += 1.0 * closedAboveSurfaceCount / aboveSurfaceCount;
  return fillingRatio;
  // return 1.0 * closedCount / rightFeatures.size();
}

int main(const int argc, const char *const *const argv) {
  using NumericType = double;

  if (argc != 3) {
    std::cout << "Usage: " << argv[0]
              << " <dataset csv file> <config "
                 "file>\n";
    return EXIT_SUCCESS;
  }

  std::string dataFilename = "trenchfill.csv";
  if (argc > 1)
    dataFilename = argv[1];

  // Load the data captured from previous
  // simulations
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

  // The input dimension is provided by the csv
  // file named parameter `InputDim` In our case
  // it is 4 since we use aspectRatio,
  // leftTaperAngle, stickingProbability and
  // timeStep as input parameters.
  int inputDim{};
  auto namedParams = reader.getNamedParameters();
  if (auto id = namedParams.find("InputDimension"); id != namedParams.end()) {
    inputDim = static_cast<int>(std::round(id->second));
  } else {
    lsMessage::getInstance()
        .addError("`InputDimension` not found "
                  "in provided data CSV file.")
        .print();
  }

  // The dimension of the data that is to be
  // interpolated. In our case this are the
  // extracted dimensions at different timesteps
  // (the timesteps are also included in the
  // data itself)
  int numFeatures = data->at(0).size() - inputDim;

  // First convert the input data into dataset
  // containing filling ratio
  auto fillingRatioData =
      lsSmartPointer<std::vector<std::vector<NumericType>>>::New();
  {
    unsigned numFeaturesRight =
        static_cast<unsigned>(std::ceil(1.0 * numFeatures / 2));

    fillingRatioData->reserve(data->size());
    for (auto &d : *data) {
      NumericType aspectRatio = d[0];
      NumericType taperAngle = d[1];
      NumericType stickingProbability = d[2];
      NumericType time = d[3];

      // Skip all entries that are not at the
      // maximum recorded time
      if (time != 1.0)
        continue;
      const auto rightFeatures =
          nonstd::span(std::next(d.begin(), inputDim),
                       std::next(d.begin(), inputDim + numFeaturesRight));
      const auto leftFeatures = nonstd::span(
          std::next(d.begin(), inputDim + numFeaturesRight), d.end());

      std::cout << aspectRatio << ", " << taperAngle << ", "
                << stickingProbability << ", ";
      auto fillingRatio = toFillingRatio<NumericType>(
          aspectRatio, taperAngle, rightFeatures, leftFeatures);

      auto row = std::vector<NumericType>{aspectRatio, taperAngle,
                                          stickingProbability, fillingRatio};

      fillingRatioData->push_back(row);
    }
  }

  // Data now only contains: aspectRatio,
  // taperAngle, stickingProbability,
  // fillingRatio.
  // aspectRatio is fixed to a value based on
  // the provided depth and width
  // stickingProbability is fixed based on the
  // input Thus all that remains is taper angle
  // vs filling ratio.

  NumericType angleLimit = Utils::rad2deg(
      std::atan((params.topWidthLimit - params.trenchBottomWidth) / 2 /
                params.trenchDepth));

  NumericType thresholdAngle = angleLimit;
  bool hasSolution = false;
  {
    SplineGridInterpolation<NumericType> sgi;
    sgi.setDataDimensions(3, 1);
    sgi.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    sgi.setData(fillingRatioData);
    sgi.initialize();
    auto uniqueValues = sgi.getUniqueValues();
    auto angles = uniqueValues[1];

    NumericType minAngle = *angles.begin();
    NumericType maxAngle = std::min(*angles.rbegin(), angleLimit);
    NumericType aspectRatio = params.trenchDepth / params.trenchBottomWidth;

#if SAMPLE_SLICE
    std::ofstream slice("slice_new.csv");

    size_t resolution = 30;
    NumericType minSticking = *uniqueValues[2].begin();
    NumericType maxSticking = *uniqueValues[2].rbegin();
    for (unsigned i = 0; i < resolution; ++i) {
      for (unsigned j = 0; j < resolution; ++j) {
        NumericType taperAngle =
            minAngle + i * (*angles.rbegin() - minAngle) / (resolution - 1);
        NumericType stickingProbability =
            minSticking + j * (maxSticking - minSticking) / (resolution - 1);
        std::vector<NumericType> loc = {aspectRatio, taperAngle,
                                        stickingProbability};

        auto [frEstimate, isInside] = sgi.estimate(loc).value();
        slice << taperAngle << ',' << stickingProbability << ','
              << frEstimate[0] << '\n';
      }
    }
    // for (auto d : *fillingRatioData) {
    //   slice << Utils::join(d.begin(), d.end()) << '\n';
    // }
    slice.close();
#endif
  }
}

//     // Do binary search to find the taper angle where the filling ratio has
//     a
//     // certain threshold value.
//     NumericType eps = 1e-5;
//     unsigned maxIter = 100;
//     unsigned i = 0;
//     for (; i < maxIter; ++i) { // Do a maximum of `maxIter` iterations
//       NumericType centerAngle = (maxAngle + minAngle) / 2;
//       std::vector<NumericType> loc = {aspectRatio, centerAngle,
//                                       params.stickingProbability, 1.0};
//       auto [frEstimate, isInside] = sgi.estimate(loc).value();
//       if (!isInside) {
//         lsMessage::getInstance()
//             .addError("The provided parameters are not within the "
//                       "boundaries of the dataset.")
//             .print();
//       }
//       if (frEstimate.at(0) <= 1.0) {
//         minAngle = centerAngle;
//       } else {
//         hasSolution = true;
//         maxAngle = centerAngle;
//       }
//       if (maxAngle - minAngle < eps) {
//         lsMessage::getInstance()
//             .addDebug(std::string("binary search iterations: ") +
//                       std::to_string(i))
//             .print();
//         break;
//       }
//     }
//     if (i >= maxIter - 1) {
//       lsMessage::getInstance()
//           .addError(std::string("No convergence in ") +
//                     std::to_string(maxIter) + std::string(" iterations"))
//           .print();
//     }
//     if (!hasSolution) {
//       lsMessage::getInstance()
//           .addError(
//               "There exists no solution in the range spanned by the
//               dataset.")
//           .print();
//       return EXIT_SUCCESS;
//     }
//     thresholdAngle = maxAngle;

//     // Follow the isoline at the provided filling ratio threshold value

//     NumericType eta = 0.01;
//     std::vector<NumericType> loc = {aspectRatio, thresholdAngle,
//                                     params.stickingProbability,
//                                     normalizedProcessTime};

//     maxIter = 5000;
//     for (unsigned i = 0; i < maxIter; ++i) {
//       if (loc[1] >= angleLimit) {
//         lsMessage::getInstance()
//             .addDebug("breaking based on angle limit")
//             .print();
//         break;
//       }

//       auto [frGrad, _] = sgi.gradient(loc).value();
//       NumericType dx = frGrad[3].front();
//       NumericType dy = frGrad[1].front();
//       NumericType ny = dx / std::sqrt(dx * dx + dy * dy);
//       NumericType nx = -std::sqrt(1.0 - ny * ny);
//       loc[3] += eta * nx;
//       loc[1] += eta * ny;
//       normalizedProcessTime = loc[3];
//       thresholdAngle = loc[1];
//       if (loc[1] > maxAngle + params.taperAngleTolerance) {
//         lsMessage::getInstance()
//             .addDebug("breaking based on taper angle tolerance")
//             .print();
//         std::cout << std::get<0>(sgi.estimate(loc).value()).front()
//                   << std::endl;
//         break;
//       }
//     }
//   }

//   // Calculate the resulting top width
//   NumericType topWidth =
//       params.trenchBottomWidth +
//       2.0 * params.trenchDepth * std::tan(Utils::deg2rad(thresholdAngle));
//   bool success = topWidth <= params.topWidthLimit;
//   // Print the result
//   std::cout << (success ? "\nSUCCESS!\n" : "\nNO SOLUTION FOUND!\n")
//             << "The minimum taper angle required to achieve a filling ratio
//             of "
//                "100% is "
//             << thresholdAngle << "°. This leads to a trench top width of "
//             << topWidth << " which "
//             << (success ? "fulfills " : "does not fulfill ")
//             << "the specified top width constraint of " <<
//             params.topWidthLimit
//             << ".\n";
//   if (success)
//     std::cout << "Assuming a top deposition rate of 1, the minimum required
//     "
//                  "process time is "
//               << topWidth * normalizedProcessTime << ".\n";
//   std::cout << std::endl;

// #ifdef WITH_VIENNAPS
//   static constexpr int D = 2;
//   if (success) {
//     std::cout << "Should I run a physical simulation to verify the "
//                  "solution? [y/n]";
//   } else {
//     std::cout << "Should I nevertheless run a physical simulation? [y/n]";
//   }

//   char c;
//   std::cin >> c;
//   if (c == 'y' || c == 'Y') {
//     std::cout << "Running simulation...\n";
//     const NumericType gridDelta = params.trenchDepth / 250;

//     // The horizontal extent is defined by the trench top width, or the
//     bottom
//     // width, whichever is larger. Extend it by 20%.
//     NumericType horizontalExtent = topWidth * 1.2;

//     // Generate the grid
//     std::array<NumericType, 3> origin{0.};
//     auto grid = createGrid<NumericType, D>(origin, gridDelta,
//     horizontalExtent,
//                                            10.0 /* initial vertical extent
//                                            */, false /* No periodic BC */);

//     // Create the inside of the trench that will be removed from the plane
//     auto cutout = createTrenchStamp<NumericType, D>(
//         grid, origin, params.trenchDepth, topWidth, thresholdAngle,
//         thresholdAngle, true);

//     auto trench = createPlane<NumericType, D>(grid, origin);
//     lsBooleanOperation<NumericType, D>(
//         trench, cutout, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
//         .apply();

//     // Since we assume a top rate of 1, the maximum time it takes for the
//     // trench to close is equal to the trench top width (if sticking
//     // probability 1, otherwise the time until pinchoff should be even
//     less)
//     // Additionally ormalize the time scale to the sticking probability, so
//     // that we get the same top layer thickness for different sticking
//     // probabilities.
//     NumericType processDuration =
//         normalizedProcessTime * topWidth / params.stickingProbability;

//     auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
//     geometry->insertNextLevelSet(trench);

//     // copy top layer to capture deposition
//     auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(trench);
//     geometry->insertNextLevelSet(depoLayer);

//     // Instantiate the process
//     auto processModel =
//         SimpleDeposition<NumericType, D>(params.stickingProbability,
//                                          1.0 /* particle source power */)
//             .getProcessModel();

//     psProcess<NumericType, D> process;
//     process.setDomain(geometry);
//     process.setProcessModel(processModel);
//     process.setNumberOfRaysPerPoint(2000);
//     process.setProcessDuration(processDuration);

//     // Run the process
//     process.apply();

//     std::cout << "Saving result as 'FilledTrench_volume.vtu'\n";
//     lsWriteVisualizationMesh<NumericType, D> visMesh;
//     for (auto ls : *geometry->getLevelSets()) {
//       visMesh.insertNextLevelSet(ls);
//     }
//     visMesh.setFileName("FilledTrench");
//     visMesh.apply();
//   }
// #endif

//   return EXIT_SUCCESS;
// }