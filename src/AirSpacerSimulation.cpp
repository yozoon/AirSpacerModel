#include <iostream>
#include <unordered_map>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsMakeGeometry.hpp>
#include <lsSmartPointer.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include <SimpleDeposition.hpp>
#include <psDomain.hpp>
#include <psProcess.hpp>

#include "TrenchGeometry.hpp"
#include "Utils.hpp"

template <typename NumericType> struct Parameters {
  NumericType aspectRatio = 15.0;
  NumericType leftTaperAngle = 0.0;
  NumericType stickingProbability = 0.1;

  void fromMap(std::unordered_map<std::string, std::string> &m) {
    Utils::AssignItems(
        m, Utils::Item{"leftTaperAngle", leftTaperAngle},
        Utils::Item{"aspectRatio", aspectRatio,
                    [](const std::string &s) -> NumericType {
                      auto value = Utils::convert<NumericType>(s);
                      if (value <= 0.0)
                        throw std::invalid_argument(
                            "'aspectRatio' must be strictly positive.");
                      return value;
                    }},
        Utils::Item{
            "stickingProbability", stickingProbability,
            [](const std::string &s) -> NumericType {
              auto value = Utils::convert<NumericType>(s);
              if (value > 1.0 || value <= 0.0)
                throw std::invalid_argument(
                    "'stickingProbability' must be in the range [1,0).");
              return value;
            }});
  }
};

int main(int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;

  const NumericType gridDelta = 0.08;

  // Parameters given by the N7 FinFET specification
  const NumericType initialTrenchTopWidth = 6.0;
  const NumericType initialTrenchDepth = 20.0;

  Parameters<NumericType> params;
  if (argc > 1) {
    auto config = Utils::readConfigFile(argv[1]);
    if (config.empty()) {
      std::cerr << "Empty config provided" << std::endl;
      return -1;
    }
    params.fromMap(config);
  }

  NumericType trenchTopWidth =
      (2 * initialTrenchDepth - initialTrenchTopWidth) /
      (2 * params.aspectRatio - 1);

  NumericType trenchDepth = params.aspectRatio * trenchTopWidth;

  NumericType conformalLayerThickness = initialTrenchDepth - trenchDepth;

  std::cout << "trenchTopWidth=" << trenchTopWidth
            << ", trenchDepth=" << trenchDepth << std::endl;

  std::array<NumericType, 3> origin{0.};

  const NumericType xExtent = initialTrenchTopWidth + 4.0;
  const NumericType yExtent = 5.0;

  const auto [leftBound, rightBound] = calculateHorizontalTrenchBounds(
      initialTrenchTopWidth, initialTrenchDepth, params.leftTaperAngle, 0.0);

  NumericType requiredExtent =
      2.0 * std::max(-leftBound, rightBound) + 2.0 * gridDelta;

  NumericType horizontalExtent = xExtent;
  if (horizontalExtent < requiredExtent) {
    lsMessage::getInstance()
        .addWarning("makeTrench: due to the provided tapering angles, the "
                    "extent of the simulation domain has to be extended.")
        .print();
    horizontalExtent = requiredExtent;
  }

  auto geometry = lsSmartPointer<
      std::vector<lsSmartPointer<lsDomain<NumericType, D>>>>::New();
  geometry->reserve(4);

  std::array<NumericType, 3> baseOrigin = origin;
  baseOrigin[D - 1] -= initialTrenchDepth;

  auto grid = createGrid<NumericType, D>(origin, gridDelta, horizontalExtent,
                                         yExtent, false /* no periodic bc */);
  geometry->emplace_back(createPlane<NumericType, D>(grid, baseOrigin));

  // Create the initial trench geometry
  auto substrate = createPlane<NumericType, D>(grid, origin);
  auto cutout = createTrenchStamp<NumericType, D>(
      grid, origin, initialTrenchDepth, initialTrenchTopWidth,
      params.leftTaperAngle, 0.0);

  lsBooleanOperation<NumericType, D>(
      substrate, cutout, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  geometry->push_back(substrate);

  // Deposit the conformal liner
  auto conformalLayer =
      lsSmartPointer<lsDomain<NumericType, D>>::New(substrate);

  auto dist = lsSmartPointer<lsSphereDistribution<NumericType, D>>::New(
      conformalLayerThickness, gridDelta);

  lsGeometricAdvect<NumericType, D>(conformalLayer, dist).apply();

  // Do CMP to remove top of conformal layer
  {
    auto plane =
        lsSmartPointer<lsDomain<NumericType, D>>::New(substrate->getGrid());

    NumericType normal[D];
    normal[0] = 0.0;
    normal[D - 1] = -1.0;
    lsMakeGeometry<NumericType, D>(
        plane,
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin.data(), normal))
        .apply();

    lsBooleanOperation<NumericType, D>(
        conformalLayer, plane, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
  geometry->push_back(conformalLayer);

  // Create non-conformal layer
  auto nonConformalLayer =
      lsSmartPointer<lsDomain<NumericType, D>>::New(conformalLayer);

  auto psGeom = psSmartPointer<psDomain<NumericType, D>>::New();
  psGeom->insertNextLevelSet(nonConformalLayer);

  auto processModel =
      SimpleDeposition<NumericType, D>(
          params.stickingProbability /* particle sticking probability */,
          1.0 /* particle source power */)
          .getProcessModel();

  // Since we assume a top rate of 1, the maximum time it takes for the
  // trench to close is equal to the trench top width (if sticking
  // probability 1, otherwise the time until pinchoff will even be less)
  NumericType processDuration = trenchTopWidth / params.stickingProbability;

  psProcess<NumericType, D> process;
  process.setDomain(psGeom);
  process.setProcessModel(processModel);
  process.setNumberOfRaysPerPoint(2000);
  process.setProcessDuration(processDuration);

  // Run the process
  process.apply();

  // Do CMP to remove top of non-conformal layer
  {
    auto plane =
        lsSmartPointer<lsDomain<NumericType, D>>::New(substrate->getGrid());

    NumericType normal[D];
    normal[0] = 0.0;
    normal[D - 1] = -1.0;
    lsMakeGeometry<NumericType, D>(
        plane,
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin.data(), normal))
        .apply();

    lsBooleanOperation<NumericType, D>(
        nonConformalLayer, plane, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
  geometry->push_back(nonConformalLayer);

  // Export the resulting geometry as a visualization mesh
  std::cout << "Writing visualization mesh...\n";
  lsWriteVisualizationMesh<NumericType, D> visMesh;
  for (auto ls : *geometry) {
    visMesh.insertNextLevelSet(ls);
  }
  visMesh.setFileName("AirSpacerSimulation");
  visMesh.apply();
}