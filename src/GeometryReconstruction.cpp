#include <array>
#include <iostream>

#include <lsDomain.hpp>
#include <lsMesh.hpp>
#include <lsSmartPointer.hpp>
#include <lsToDiskMesh.hpp>

#include <SimpleDeposition.hpp>
#include <psProcess.hpp>
#include <psUtils.hpp>

#include "models/MakeTrenchStamp.hpp"
#include "simulation/AdvectionCallback.hpp"
#include "simulation/FeatureExtraction.hpp"

#include "MakeTrench.hpp"

template <typename NumericType, int D> struct Parameters {
  NumericType gridDelta = 0.2;
  NumericType trenchTopWidth = 5.0;
  NumericType trenchDepth = 50.0;
  NumericType leftTaperAngle = 0.0;
  NumericType rightTaperAngle = 0.0;
  NumericType stickingProbability = 0.1;
  NumericType processDuration = 5.0;

  void fromMap(std::unordered_map<std::string, std::string> &m) {
    psUtils::AssignItems(                                          //
        m,                                                         //
        psUtils::Item{"gridDelta", gridDelta},                     //
        psUtils::Item{"trenchDepth", trenchDepth},                 //
        psUtils::Item{"trenchTopWidth", trenchTopWidth},           //
        psUtils::Item{"leftTaperAngle", leftTaperAngle},           //
        psUtils::Item{"rightTaperAngle", rightTaperAngle},         //
        psUtils::Item{"stickingProbability", stickingProbability}, //
        psUtils::Item{"processDuration", processDuration}          //
    );
  }
};

int main(int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;
  static constexpr int numberOfSamples = 512;

  Parameters<NumericType, D> params;
  if (argc > 1) {
    auto config = psUtils::readConfigFile(argv[1]);
    if (config.empty()) {
      std::cerr << "Empty config provided" << std::endl;
      return -1;
    }
    params.fromMap(config);
  }

  const NumericType leftOffset =
      std::tan(params.leftTaperAngle * M_PI / 180.) * params.trenchDepth;
  const NumericType rightOffset =
      std::tan(params.rightTaperAngle * M_PI / 180.) * params.trenchDepth;

  // Check if the tapering would interfere with the depth of the trench
  if (params.trenchTopWidth / 2. - leftOffset + rightOffset <=
      params.gridDelta) {
    std::cout << "Intersection in trench tapering detected. Skipping.\n";
    return EXIT_FAILURE;
  }

  // Make sure that the trench sidewalls stay inside the simulation
  // domain, even if they are tapered.
  NumericType xExtent =
      2. * (params.trenchTopWidth / 2 + std::max(leftOffset, NumericType{0.}) +
            std::max(rightOffset, NumericType{0.}) + 5. * params.gridDelta);

  std::array<NumericType, 3> origin{0.};
  // Generate the initial trench geometry
  auto levelset = MakeTrench<NumericType, D>(
      params.gridDelta, xExtent, 0., origin, params.trenchTopWidth,
      params.trenchDepth, params.leftTaperAngle, params.rightTaperAngle,
      false /* no periodic boundary*/);

  auto geometry = lsSmartPointer<psDomain<NumericType, D>>::New();
  geometry->insertNextLevelSet(levelset);
  geometry->printSurface("GR_initial.vtp");

  // Run a physical deposition simulation
  auto processModel =
      SimpleDeposition<NumericType, D>(
          params.stickingProbability /* particle sticking probability */,
          1.0 /* particle source power */)
          .getProcessModel();

  psProcess<NumericType, D> process;
  process.setDomain(geometry);
  process.setProcessModel(processModel);
  process.setNumberOfRaysPerPoint(2000);
  process.setProcessDuration(params.processDuration /
                             params.stickingProbability);

  // Run the process
  process.apply();

  geometry->printSurface("GR_simulation.vtp");

  // Extract features from the geometry
  FeatureExtraction<NumericType, D> extraction;
  extraction.setDomain(geometry->getLevelSets()->back());
  extraction.setNumberOfSamples(numberOfSamples,
                                false /* don't include endpoints */);
  extraction.setEdgeAffinity(0.0);
  extraction.setOrigin(origin);
  extraction.setTrenchDimensions(params.trenchDepth, params.trenchTopWidth);
  extraction.apply();

  auto sampleLocations = extraction.getSampleLocations();
  auto features = extraction.getFeatures();

#ifndef NDEBUG
  std::cout << "Number of features=" << features->size() << std::endl;
  for (unsigned i = 0; i < sampleLocations->size(); ++i) {
    std::cout << i << ": " << std::setprecision(4)
              << -params.trenchDepth +
                     (params.trenchDepth + params.trenchTopWidth) *
                         sampleLocations->at(i)
              << ", " << features->at(i) << '\n';
  }
#endif
  // Now reconstruct the geometry based on the extracted features
  auto stamp =
      MakeTrenchStamp(levelset->getGrid(), origin, params.trenchDepth,
                      params.trenchTopWidth, *sampleLocations, *features);

  auto reconstructedGeometry =
      lsSmartPointer<lsDomain<NumericType, D>>::New(levelset->getGrid());
  NumericType planeOrigin[D];
  planeOrigin[0] = origin[0];
  planeOrigin[1] = origin[1];
  if constexpr (D == 3)
    planeOrigin[2] = origin[2];
  planeOrigin[D - 1] = params.processDuration;

  NumericType normal[D];
  normal[D - 1] = 1;

  lsMakeGeometry<NumericType, D>(
      reconstructedGeometry,
      lsSmartPointer<lsPlane<NumericType, D>>::New(planeOrigin, normal))
      .apply();

  lsBooleanOperation<NumericType, D>(
      reconstructedGeometry, stamp, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  Utils::printSurface(stamp, "stamp.vtp");
  Utils::printSurface(reconstructedGeometry, "GR_reconstructed.vtp");
}