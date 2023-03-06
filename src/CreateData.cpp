#include <iostream>
#include <vector>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include <SimpleDeposition.hpp>
#include <psProcess.hpp>

#include "AdvectionCallback.hpp"
#include "FeatureExtraction.hpp"
#include "TrenchGeometry.hpp"
#include "Utils.hpp"

int main(int, const char *const *const) {
  using NumericType = double;
  static constexpr int D = 2;

  std::array<NumericType, 3> origin{0.};
  NumericType gridDelta = 0.2;
  NumericType trenchTopWidth = 4.; // -> 20 grid points

  int numberOfSamples = 512;

  // The parameters we are interested in
  std::vector<NumericType> aspectRatios = {10., 30., 50., 70., 100.};
  std::vector<NumericType> leftTaperAngles = {-2.0, -1.0, 0.0, 2.0};
  // std::vector<NumericType> rightTaperAngles = {0.0};
  std::vector<NumericType> stickingProbabilities = {0.1, 0.07, 0.04, 0.01};

  unsigned totalCombinations = aspectRatios.size() * leftTaperAngles.size() *
                               stickingProbabilities.size();

  // The data we are going to store consists of stickingProbability,
  // taperAngle, time and the sampled geometry descriptors as provided by the
  // feature extraction.
  int inputDimension = 4;

  // Instantiate the featureExtraction
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

  const std::string filename = "data.csv";
  auto writer = lsSmartPointer<psCSVWriter<NumericType>>::New(filename);

  // Creation of a descriptive/ detailed header
  std::string header =
      "aspectRatio,leftTaperAngle,stickingProbability,timestep";
  for (int i = 0; i < numberOfSamples; ++i) {
    header += ",diameter_" + std::to_string(i);
  }
  header += "\nDiameters are normalized to the trench width.\nSample positions "
            "along vertical axis (normalized to trench depth + trench width):";
  header += "\n!" + join(sampleLocations->begin(), sampleLocations->end());
  header += "\n!InputDimension=" + std::to_string(inputDimension);

  writer->setHeader(header);
  writer->initialize();

  unsigned count = 1;
  for (auto aspectRatio : aspectRatios) {
    // Calculate the depth of the trench using the aspect ratio
    NumericType trenchDepth = trenchTopWidth * aspectRatio;
    for (auto leftTaperAngle : leftTaperAngles) {
      // for (auto rightTaperAngle : rightTaperAngles) {
      NumericType rightTaperAngle = 0.0;

      NumericType xExtent = 2.0 * trenchTopWidth;

      // Now that we know all geometry parameters generate the geometry
      auto trench = makeTrench<NumericType, D>(
          gridDelta, xExtent, 10., origin, trenchTopWidth, trenchDepth,
          leftTaperAngle, rightTaperAngle, false /* no periodic boundary*/);

      for (auto stickingProbability : stickingProbabilities) {
        std::cout << count << '/' << totalCombinations << std::endl;
        ++count;

        // Normalize the time scale to the sticking probability, so that we
        // get the same top layer thickness for different sticking
        // probabilities.
        NumericType timeScale = 1.0 / stickingProbability;

        // Ensure that we always have 11 samples
        NumericType extractionInterval = trenchTopWidth / 10.;

        featureExtraction->setTrenchDimensions(trenchDepth, trenchTopWidth);

        auto advectionCallback = lsSmartPointer<AdvectionCallback<
            NumericType, D, decltype(featureExtraction)::element_type,
            decltype(writer)::element_type>>::New(timeScale,
                                                  extractionInterval);

        advectionCallback->setFeatureExtraction(featureExtraction);
        advectionCallback->setPrefixData(std::vector<NumericType>{
            aspectRatio, leftTaperAngle, stickingProbability});
        advectionCallback->setWriter(writer);

        auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
        geometry->insertNextLevelSet(trench);

        // copy top layer to capture deposition
        auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(trench);
        geometry->insertNextLevelSet(depoLayer);

        auto processModel =
            SimpleDeposition<NumericType, D>(
                stickingProbability /* particle sticking probability */,
                1.0 /* particle source power */)
                .getProcessModel();

        processModel->setAdvectionCallback(advectionCallback);

        // Since we assume a top rate of 1, the maximum time it takes for the
        // trench to close is equal to the trench top width (if sticking
        // probability 1, otherwise the time until pinchoff will even be less)
        NumericType processDuration = trenchTopWidth / stickingProbability;

        psProcess<NumericType, D> process;
        process.setDomain(geometry);
        process.setProcessModel(processModel);
        process.setNumberOfRaysPerPoint(2000);
        process.setProcessDuration(processDuration);

        // Run the process
        process.apply();

        writer->flush();
        // }
      }
    }
  }

  return EXIT_SUCCESS;
}