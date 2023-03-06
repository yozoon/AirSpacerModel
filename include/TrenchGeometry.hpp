#ifndef TRENCH_GEOMETRY_HPP
#define TRENCH_GEOMETRY_HPP

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include "Utils.hpp"

template <class NumericType>
std::tuple<NumericType, NumericType> calculateHorizontalTrenchBounds(
    const NumericType trenchTopWidth, const NumericType trenchDepth,
    const NumericType leftTaperAngle, const NumericType rightTaperAngle) {
  const NumericType tanLeft =
      std::tan(-leftTaperAngle * Utils::PI<NumericType> / 180.);
  const NumericType tanRight =
      std::tan(-rightTaperAngle * Utils::PI<NumericType> / 180.);

  const NumericType leftOffset = tanLeft * trenchDepth;
  const NumericType rightOffset = tanRight * trenchDepth;

  return {-trenchTopWidth / 2 - std::max(leftOffset, NumericType{0.}),
          trenchTopWidth / 2 + std::max(rightOffset, NumericType{0.})};
};

template <class NumericType, int D>
hrleGrid<D> createGrid(const std::array<NumericType, 3> &origin,
                       const NumericType gridDelta, const NumericType xExtent,
                       const NumericType yExtent, bool periodicBoundary) {
  NumericType bounds[2 * D];
  bounds[0] = origin[0] - xExtent / 2.;
  bounds[1] = origin[0] + xExtent / 2.;

  if constexpr (D == 3) {
    bounds[2] = origin[1] - yExtent / 2.;
    bounds[3] = origin[1] + yExtent / 2.;
    bounds[4] = origin[2] - gridDelta;
    bounds[5] = origin[2] + gridDelta;
  } else {
    bounds[2] = origin[1] - gridDelta;
    bounds[3] = origin[1] + gridDelta;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];

  for (int i = 0; i < D - 1; i++) {
    if (periodicBoundary) {
      boundaryCons[i] =
          lsDomain<NumericType, D>::BoundaryType::PERIODIC_BOUNDARY;
    } else {
      boundaryCons[i] =
          lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    }
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  hrleIndexType gridMin[D], gridMax[D];
  for (unsigned i = 0; i < D; ++i) {
    gridMin[i] = std::floor(bounds[2 * i] / gridDelta);
    gridMax[i] = std::ceil(bounds[2 * i + 1] / gridDelta);
  }

  return hrleGrid<D>(gridMin, gridMax, gridDelta, boundaryCons);
}

template <class NumericType, int D>
lsSmartPointer<lsDomain<NumericType, D>>
createPlane(const hrleGrid<D> &grid, const std::array<NumericType, 3> &origin) {
  auto plane = lsSmartPointer<lsDomain<NumericType, D>>::New(grid);

  NumericType point[D] = {0.};
  point[0] = origin[0];
  point[1] = origin[1];
  if constexpr (D == 3)
    point[2] = origin[2];

  NumericType normal[D] = {0.};
  normal[D - 1] = 1.;

  lsMakeGeometry<NumericType, D>(
      plane, lsSmartPointer<lsPlane<NumericType, D>>::New(point, normal))
      .apply();
  return plane;
}

template <class NumericType, int D>
lsSmartPointer<lsDomain<NumericType, D>> createTrenchStamp(
    const hrleGrid<D> &grid, const std::array<NumericType, 3> &origin,
    const NumericType passedTrenchDepth, const NumericType trenchTopWidth,
    const NumericType leftTaperAngle, const NumericType rightTaperAngle) {

  const NumericType gridDelta = grid.getGridDelta();

  const NumericType tanLeft =
      std::tan(leftTaperAngle * Utils::PI<NumericType> / 180.);
  const NumericType tanRight =
      std::tan(rightTaperAngle * Utils::PI<NumericType> / 180.);

  NumericType trenchDepth = passedTrenchDepth;
  if (tanLeft + tanRight > 0) {
    NumericType intersectionDepth = trenchTopWidth / (tanLeft + tanRight);
    if (intersectionDepth < trenchDepth) {
      lsMessage::getInstance()
          .addWarning(
              "createTrenchStamp: due to the provided tapering angles, the "
              "trench can only go to a depth of " +
              std::to_string(intersectionDepth))
          .print();
      trenchDepth = intersectionDepth;
    }
  }

  const NumericType leftOffset = tanLeft * trenchDepth;
  const NumericType rightOffset = tanRight * trenchDepth;

  auto stamp = lsSmartPointer<lsDomain<NumericType, D>>::New(grid);

  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();

  if constexpr (D == 2) {
    for (int i = 0; i < 4; i++) {
      std::array<NumericType, 3> node = {0., 0., 0.};
      mesh->insertNextNode(node);
    }
    mesh->nodes[0][0] = origin[0] - trenchTopWidth / 2. + leftOffset;
    mesh->nodes[1][0] = origin[0] + trenchTopWidth / 2. - rightOffset;
    mesh->nodes[2][0] = origin[0] + trenchTopWidth / 2.;
    mesh->nodes[3][0] = origin[0] - trenchTopWidth / 2.;

    mesh->nodes[0][1] = origin[1] - trenchDepth;
    mesh->nodes[1][1] = origin[1] - trenchDepth;
    mesh->nodes[2][1] = origin[1];
    mesh->nodes[3][1] = origin[1];

    mesh->insertNextLine(std::array<unsigned, 2>{0, 3});
    mesh->insertNextLine(std::array<unsigned, 2>{3, 2});
    mesh->insertNextLine(std::array<unsigned, 2>{2, 1});
    mesh->insertNextLine(std::array<unsigned, 2>{1, 0});
    lsFromSurfaceMesh<NumericType, D>(stamp, mesh).apply();
  } else {
    auto bounds = Utils::getBoundsFromGrid<NumericType, D>(grid);
    for (int i = 0; i < 8; i++) {
      std::array<NumericType, 3> node = {0., 0., 0.};
      mesh->insertNextNode(node);
    }
    mesh->nodes[0][0] = origin[0] - trenchTopWidth / 2. + leftOffset;
    mesh->nodes[0][1] = origin[1] + bounds[2] - gridDelta;
    mesh->nodes[0][2] = origin[2] - trenchDepth;

    mesh->nodes[1][0] = origin[0] + trenchTopWidth / 2. - rightOffset;
    mesh->nodes[1][1] = origin[1] + bounds[2] - gridDelta;
    mesh->nodes[1][2] = origin[2] - trenchDepth;

    mesh->nodes[2][0] = origin[0] + trenchTopWidth / 2. - rightOffset;
    mesh->nodes[2][1] = origin[1] + bounds[3] + gridDelta;
    mesh->nodes[2][2] = origin[2] - trenchDepth;

    mesh->nodes[3][0] = origin[0] - trenchTopWidth / 2. + leftOffset;
    mesh->nodes[3][1] = origin[1] + bounds[3] + gridDelta;
    mesh->nodes[3][2] = origin[2] - trenchDepth;

    mesh->nodes[4][0] = origin[0] - trenchTopWidth / 2.;
    mesh->nodes[4][1] = origin[1] + bounds[2] - gridDelta;
    mesh->nodes[4][2] = origin[2];

    mesh->nodes[5][0] = origin[0] + trenchTopWidth / 2.;
    mesh->nodes[5][1] = origin[1] + bounds[2] - gridDelta;
    mesh->nodes[5][2] = origin[2];

    mesh->nodes[6][0] = origin[0] + trenchTopWidth / 2.;
    mesh->nodes[6][1] = origin[1] + bounds[3] + gridDelta;
    mesh->nodes[6][2] = origin[2];

    mesh->nodes[7][0] = origin[0] - trenchTopWidth / 2.;
    mesh->nodes[7][1] = origin[1] + bounds[3] + gridDelta;
    mesh->nodes[7][2] = origin[2];

    mesh->insertNextTriangle(std::array<unsigned, 3>{0, 3, 1});
    mesh->insertNextTriangle(std::array<unsigned, 3>{1, 3, 2});

    mesh->insertNextTriangle(std::array<unsigned, 3>{5, 6, 4});
    mesh->insertNextTriangle(std::array<unsigned, 3>{6, 7, 4});

    mesh->insertNextTriangle(std::array<unsigned, 3>{0, 1, 5});
    mesh->insertNextTriangle(std::array<unsigned, 3>{0, 5, 4});

    mesh->insertNextTriangle(std::array<unsigned, 3>{2, 3, 6});
    mesh->insertNextTriangle(std::array<unsigned, 3>{6, 3, 7});

    mesh->insertNextTriangle(std::array<unsigned, 3>{0, 7, 3});
    mesh->insertNextTriangle(std::array<unsigned, 3>{0, 4, 7});

    mesh->insertNextTriangle(std::array<unsigned, 3>{1, 2, 6});
    mesh->insertNextTriangle(std::array<unsigned, 3>{1, 6, 5});

    lsFromSurfaceMesh<NumericType, D>(stamp, mesh).apply();
  }
  return stamp;
}

template <class NumericType, int D>
lsSmartPointer<lsDomain<NumericType, D>>
makeTrench(const NumericType gridDelta, const NumericType xExtent,
           const NumericType yExtent, const std::array<NumericType, 3> &origin,
           const NumericType trenchTopWidth, const NumericType trenchDepth,
           const NumericType leftTaperAngle, const NumericType rightTaperAngle,
           const bool periodicBoundary) {

  const auto [leftBound, rightBound] = calculateHorizontalTrenchBounds(
      trenchTopWidth, trenchDepth, leftTaperAngle, rightTaperAngle);

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

  auto grid = createGrid<NumericType, D>(origin, gridDelta, horizontalExtent,
                                         yExtent, periodicBoundary);

  auto substrate = createPlane<NumericType, D>(grid, origin);

  auto cutout = createTrenchStamp<NumericType, D>(
      grid, origin, trenchDepth, trenchTopWidth, leftTaperAngle,
      rightTaperAngle);

  lsBooleanOperation<NumericType, D>(
      substrate, cutout, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  return substrate;
}

#endif