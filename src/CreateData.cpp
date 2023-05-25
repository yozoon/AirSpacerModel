#include <array>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include <SimpleDeposition.hpp>
#include <psProcess.hpp>

#include "AdvectionCallback.hpp"
#include "FeatureExtraction.hpp"
#include "TrenchGeometry.hpp"

#include "CreateData.hpp"

template <typename NumericType>
std::string createHeader(const int numberOfSamples,
                         const std::vector<NumericType> &sampleLocations,
                         const int inputDimension, bool symmetrical) {
  std::string header =
      "aspectRatio,leftTaperAngle,stickingProbability,normalizedTime";
  for (int i = 0; i < numberOfSamples; ++i) {
    header += ",diameter_" + std::to_string(i);
  }
  header += "\nDiameters are normalized to the trench width.\nSample positions "
            "along vertical axis (normalized to trench depth + trench width):";
  header += "\n!" + Utils::join(sampleLocations.begin(), sampleLocations.end());
  header += "\n!InputDimension=" + std::to_string(inputDimension) +
            ",Symmetrical=" + std::to_string(symmetrical);
  return header;
}

int main(const int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;

  std::array<NumericType, 3> origin{0.};
  const size_t numberOfCrossSectionPoints = 500;
  const NumericType trenchWidth = 5.;

  const int numberOfSamples = 512;

  std::unordered_map<std::string, std::string> config;

  if (argc > 1) {
    config = Utils::readConfigFile(argv[1]);
    if (config.empty()) {
      lsMessage::getInstance().addError("Empty config provided").print();
      return EXIT_FAILURE;
    }
  }
  auto params = Parameters<NumericType>::fromMap(config);
  params.print();

  // The data we are going to store consists of stickingProbability,
  // taperAngle, time and the sampled geometry descriptors as provided by the
  // feature extraction.
  int inputDimension = params.getInputDimension();

  // Instantiate the feature extraction tool
  auto featureExtraction =
      lsSmartPointer<FeatureExtraction<NumericType, D>>::New();
  featureExtraction->setNumberOfSamples(numberOfSamples,
                                        false /* don't include endpoints */);
  // Uniformly spaced sample points
  featureExtraction->setEdgeAffinity(0.0);
  featureExtraction->initializeSampleLocations();

  // The locations at which the features are extracted (normalized to
  // the trench depth at each timestep)
  auto sampleLocations = featureExtraction->getSampleLocations();

  auto writer = lsSmartPointer<psCSVWriter<NumericType>>::New(params.filename);

  // Creation of a header for the csv file
  auto header = createHeader(numberOfSamples, *sampleLocations, inputDimension,
                             params.symmetrical);
  writer->setHeader(header);
  writer->initialize();

  const auto totalNumberOfCombinations = params.getNumberOfCombinations();

  const bool fixTopWidth = params.fixTopWidth;
  const auto maxTime = params.timeSteps;

  unsigned count = 0;
  for (const auto [aspectRatio, leftTaperAngle, rightTaperAngle,
                   stickingProbability] : params) {

    std::vector<NumericType> prefixData = {aspectRatio, leftTaperAngle,
                                           stickingProbability};

    // Calculate the depth of the trench based on the aspect ratio
    NumericType trenchDepth = trenchWidth * aspectRatio;
    NumericType trenchTopWidth = trenchWidth;
    NumericType trenchBottomWidth = trenchWidth;

    // Calculate the offset caused by the tapering
    const NumericType leftOffset =
        trenchDepth * std::tan(Utils::deg2rad(leftTaperAngle));
    const NumericType rightOffset =
        trenchDepth * std::tan(Utils::deg2rad(rightTaperAngle));

    if (fixTopWidth) {
      trenchBottomWidth = trenchBottomWidth - leftOffset - rightOffset;
    } else {
      trenchTopWidth = trenchTopWidth + leftOffset + rightOffset;
    }

    // The horizontal extent is defined by the trench top width, or the bottom
    // width, whichever is larger. Extend it by 10%.
    NumericType horizontalExtent = std::max(
        std::max(trenchTopWidth, trenchBottomWidth) * 1.1, NumericType{10.});

    NumericType crossSectionLength =
        horizontalExtent - trenchTopWidth +
        std::sqrt(leftOffset * leftOffset + trenchDepth * trenchDepth) +
        std::sqrt(rightOffset * rightOffset + trenchDepth * trenchDepth) +
        trenchBottomWidth;

    // Adjustment of grid delta:
    // Calculate a grid delta such that the final geometry contains
    // approximately `numberOfCrossSectionPoints` surface points
    NumericType gridDeltaSpace =
        crossSectionLength / numberOfCrossSectionPoints;

    // Calculate another grid delta that ensures that at least
    // `numberOfStepsPerInterval` steps are taken in each extraction interval.
    int numberOfStepsPerInterval = 1;
    NumericType gridDeltaTime = stickingProbability * (params.timeSteps - 1) /
                                (numberOfStepsPerInterval * 0.4999);

    // Finally chose whichever value of the two is smaller
    const NumericType gridDelta = std::min(gridDeltaSpace, gridDeltaTime);

    std::cout << ++count << '/' << totalNumberOfCombinations
              << " complexity: " << crossSectionLength / gridDelta << std::endl;

    // If one of the widths is negative, the sidewalls intersect. This makes
    // would make the definition of the aspect ratio quite questionable, thus we
    // skip these cases alltogether.
    if (trenchTopWidth <= gridDelta || trenchBottomWidth <= gridDelta) {
      // TODO: create `timeSteps` number of rows in the csv file containing
      // placeholder values or zeros (such that the regular grid is still
      // preserved)
      lsMessage::getInstance()
          .addWarning("The provided taper angle would lead to intersecting "
                      "sidewalls. Skipping configuration.")
          .print();

      std::vector<NumericType> dummyRow(numberOfSamples + prefixData.size() +
                                        1);
      const NumericType placeholderValue = 0.0;
      // Copy the prefix data into the row vector
      std::copy(prefixData.begin(), prefixData.end(), dummyRow.begin());

      // Initialize the sample values with the placeholder value
      std::generate(std::next(dummyRow.begin(), prefixData.size() + 1UL),
                    dummyRow.end(), [=]() { return placeholderValue; });
      for (int i = 0; i < maxTime; ++i) {
        // Update the timestep value in the dummy row
        dummyRow[prefixData.size()] = 1.0 * i / (params.timeSteps - 1);
        // Write the row to the file
        writer->writeRow(dummyRow);
      }
      continue;
    }

    // Set the dimensions, so that the feature extraction can scale the values
    // appropriately
    featureExtraction->setTrenchDimensions(trenchDepth, trenchTopWidth);

    // Generate the grid
    auto grid = createGrid<NumericType, D>(origin, gridDelta, horizontalExtent,
                                           10.0 /* initial vertical extent */,
                                           false /* No periodic BC */);

    // Create the inside of the trench that will be removed from the plane
    auto cutout = createTrenchStamp<NumericType, D>(
        grid, origin, trenchDepth, trenchTopWidth, leftTaperAngle,
        rightTaperAngle, true);

    auto trench = createPlane<NumericType, D>(grid, origin);
    lsBooleanOperation<NumericType, D>(
        trench, cutout, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    // Normalize the time scale to the sticking probability, so that we
    // get the same top layer thickness for different sticking
    // probabilities.
    NumericType timeNormalizationFactor = 1.0 / stickingProbability;

    // Since we assume a top rate of 1, the maximum time it takes for the
    // trench to close is equal to the trench top width (if sticking
    // probability 1, otherwise the time until pinchoff should be even less)
    NumericType processDuration = trenchTopWidth;

    auto advectionCallback = lsSmartPointer<AdvectionCallback<
        NumericType, D, decltype(featureExtraction)::element_type,
        decltype(writer)::element_type>>::New(timeNormalizationFactor,
                                              processDuration /
                                                  (params.timeSteps - 1));

    advectionCallback->setFeatureExtraction(featureExtraction);
    advectionCallback->setPrefixData(prefixData);
    advectionCallback->setWriter(writer);
    advectionCallback->setModifiers(1.0 / (params.timeSteps - 1), 1.0);

    auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
    geometry->insertNextLevelSet(trench);

    // copy top layer to capture deposition
    auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(trench);
    geometry->insertNextLevelSet(depoLayer);

    // Instantiate the process
    auto processModel = psSmartPointer<SimpleDeposition<NumericType, D>>::New(
        stickingProbability, 1.0 /* particle source power */);

    processModel->setAdvectionCallback(advectionCallback);

    psProcess<NumericType, D> process;
    process.setDomain(geometry);
    process.setProcessModel(processModel);
    process.setNumberOfRaysPerPoint(1000);
    process.setProcessDuration(processDuration * timeNormalizationFactor);

    // Run the process
    process.apply();

    // Ensure that values are written to disk
    writer->flush();
  }
}