#include <iostream>
#include <unordered_map>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>

#include <lsSmartPointer.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include "TrenchGeometry.hpp"
#include "Utils.hpp"
#include "GeometricAirSpacerModel.hpp"

template <typename NumericType> struct Parameters {
  NumericType aspectRatio = 15.0;
  NumericType leftTaperAngle = 0.0;
  NumericType stickingProbability = 0.1;
  NumericType offset = 0.0;

  void fromMap(std::unordered_map<std::string, std::string> &m) {
    Utils::AssignItems(
        m, Utils::Item{"leftTaperAngle", leftTaperAngle},
        Utils::Item{"aspectRatio", aspectRatio,
                    [](const std::string &s) -> NumericType {
                      auto value = Utils::convert<NumericType>(s);
                      if (value <= 0.0)
                        throw std::invalid_argument(
                            "`aspectRatio` must be strictly positive.");
                      return value;
                    }},
        Utils::Item{
            "stickingProbability", stickingProbability,
            [](const std::string &s) -> NumericType {
              auto value = Utils::convert<NumericType>(s);
              if (value > 1.0 || value <= 0.0)
                throw std::invalid_argument(
                    "`stickingProbability` must be in the range [1,0).");
              return value;
            }},
        Utils::Item{"offset", offset});
  }
};

int main(int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;

  const NumericType gridDelta = .08;

  // Parameters given by the N7 FinFET specification
  const NumericType initialTrenchTopWidth = 6.0;
  const NumericType initialTrenchDepth = 20.0;

  std::string dataFile = "data.csv";
  if (argc > 1)
    dataFile = argv[1];

  Parameters<NumericType> params;
  if (argc > 2) {
    auto config = Utils::readConfigFile(argv[2]);
    if (config.empty()) {
      lsMessage::getInstance().addError("Empty config provided").print();
      return EXIT_FAILURE;
    }
    params.fromMap(config);
  }

  const NumericType initialAspectRatio =
      initialTrenchDepth / initialTrenchTopWidth;

  if (params.aspectRatio < initialAspectRatio) {
    lsMessage::getInstance()
        .addError(
            std::string("The aspect ratio must not be smaller than the initial "
                        "aspect ratio: ") +
            std::to_string(initialTrenchDepth / initialTrenchTopWidth))
        .print();
    return EXIT_FAILURE;
  }

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
  geometry->reserve(3);

  auto grid = createGrid<NumericType, D>(origin, gridDelta, horizontalExtent,
                                         yExtent, false /* no periodic bc */);
  geometry->emplace_back(createPlane<NumericType, D>(grid, origin));

  // Apply the emulation model to the geometry
  GeometricAirSpacerModel<NumericType, D>(
      geometry, dataFile, origin, initialTrenchDepth, initialTrenchTopWidth,
      params.aspectRatio, params.leftTaperAngle, params.stickingProbability,
      params.offset)
      .apply();

  // Export the resulting geometry as a visualization mesh
  lsMessage::getInstance().addDebug("Writing visualization mesh...").print();
  lsWriteVisualizationMesh<NumericType, D> visMesh;
  for (auto ls : *geometry) {
    visMesh.insertNextLevelSet(ls);
  }
  visMesh.setFileName("AirGapEmulation");
  visMesh.apply();
}