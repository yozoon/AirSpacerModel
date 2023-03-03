#include <iostream>
#include <vector>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include <SimpleDeposition.hpp>
#include <psProcess.hpp>

#include "MakeTrench.hpp"
#include "Utils.hpp"
#include "simulation/AdvectionCallback.hpp"
#include "simulation/FeatureExtraction.hpp"

int main(int, const char *const *const) {
  using NumericType = double;
  static constexpr int D = 3;

  std::array<NumericType, 3> origin{0.};
  NumericType gridDelta = 0.2;
  NumericType minXExtent = 10.;
  NumericType trenchDepth = 8.;

  int numberOfSamples = 201;

  // How long to run the process and at which intervals to do the extraction
  NumericType processDuration = 10.0;

  // The parameters we are interested in
  std::vector<NumericType> aspectRatios = {10., 30., 60.};
  std::vector<NumericType> leftTaperAngles = {-2.0, -1.0, 0.0, 1.0, 2.0};
  std::vector<NumericType> stickingProbabilities = {0.1, 0.07, 0.04, 0.01};

  unsigned totalCombinations = aspectRatios.size() * leftTaperAngles.size() *
                               stickingProbabilities.size();

  // The data we are going to store consists of stickingProbability,
  // taperAngle, time and the sampled geometry descriptors as provided by the
  // feature extraction.
  int InputDimension = 3;

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
  std::string header = "aspectRatio,leftTaperAngle,stickingProbability,";
  header += "time,depth";
  for (int i = 0; i < numberOfSamples - 1; ++i)
    header += ",diameter_" + std::to_string(i);

  header += "\nDiameter sample distribution along vertical axis (normalized to "
            "trench minmax):";
  header += "\n!" + join(sampleLocations->begin(), sampleLocations->end());
  header += "\n!InputDimension=" + std::to_string(InputDimension) +
            ",TrenchDepth=" + std::to_string(trenchDepth);

  writer->setHeader(header);
  writer->initialize();

  unsigned count = 1;
  for (auto aspectRatio : aspectRatios) {
    // Calculate the width of the trench using the aspect ratio
    NumericType trenchTopWidth = trenchDepth / aspectRatio;
    for (auto leftTaperAngle : leftTaperAngles) {
      NumericType rightTaperAngle = 0.0;
      const NumericType leftOffset =
          std::tan(leftTaperAngle * M_PI / 180.) * trenchDepth;

      // Check if the tapering would interfere with the depth of the trench
      if (trenchTopWidth / 2. - leftOffset <= gridDelta) {
        std::cout << "Intersection in trench tapering detected. Skipping.\n";
        ++count;
        continue;
      }

      // Make sure that the trench sidewalls stay inside the simulation
      // domain, even if they are tapered.
      NumericType xExtent =
          2. *
          std::max(trenchTopWidth / 2 + leftOffset + gridDelta, minXExtent / 2);

      // Here we already know all geometry parameters, so generate the geometry
      auto trench = MakeTrench<NumericType, D>(
          gridDelta, xExtent, 10., origin, trenchTopWidth, trenchDepth,
          leftTaperAngle, rightTaperAngle, false /* no periodic boundary*/);

      auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
      geometry->insertNextLevelSet(trench);

      for (auto stickingProbability : stickingProbabilities) {
        std::cout << count << '/' << totalCombinations << std::endl;
        ++count;

        // Normalize the time scale to the sticking probability, so that we get
        // the same top layer thickness for different sticking probabilities.
        NumericType timeScale = 1.0 / stickingProbability;

        auto advectionCallback = lsSmartPointer<AdvectionCallback<
            NumericType, D, decltype(featureExtraction)::element_type,
            decltype(writer)::element_type>>::
            New(timeScale, 1.0 /* extraction interval */);

        advectionCallback->setFeatureExtraction(featureExtraction);
        advectionCallback->setPrefixData(std::vector<NumericType>{
            aspectRatio, leftTaperAngle, stickingProbability});
        advectionCallback->setWriter(writer);

        // copy top layer to capture deposition
        auto depoLayer = psSmartPointer<lsDomain<NumericType, D>>::New(
            geometry->getLevelSets()->back());
        geometry->insertNextLevelSet(depoLayer);

        auto processModel =
            SimpleDeposition<NumericType, D>(
                stickingProbability /* particle sticking probability */,
                1.0 /* particle source power */)
                .getProcessModel();

        processModel->setAdvectionCallback(advectionCallback);

        psProcess<NumericType, D> process;
        process.setDomain(geometry);
        process.setProcessModel(processModel);
        process.setNumberOfRaysPerPoint(2000);
        process.setProcessDuration(processDuration / stickingProbability);

        // Run the process
        process.apply();

        writer->flush();
      }
    }
  }

  return EXIT_SUCCESS;
}