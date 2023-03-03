#ifndef MAKE_TRENCH_HPP
#define MAKE_TRENCH_HPP

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsSmartPointer.hpp>

template <class NumericType, int D>
lsSmartPointer<lsDomain<NumericType, D>>
MakeTrench(const NumericType gridDelta, const NumericType xExtent,
           const NumericType yExtent, const std::array<NumericType, 3> &origin,
           const NumericType trenchTopWidth, const NumericType trenchDepth,
           const NumericType leftTaperingAngle,
           const NumericType rightTaperingAngle, const bool periodicBoundary) {
  const static constexpr NumericType PI = std::acos(NumericType{-1});

  NumericType bounds[2 * D];
  bounds[0] = origin[0] - xExtent / 2.;
  bounds[1] = origin[0] + xExtent / 2.;

  if constexpr (D == 3) {
    bounds[2] = origin[1] - yExtent / 2.;
    bounds[3] = origin[1] + yExtent / 2.;
    bounds[4] = origin[2] - trenchDepth - gridDelta;
    bounds[5] = origin[2] + gridDelta;
  } else {
    bounds[2] = origin[1] - trenchDepth - gridDelta;
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

  // Create an empty levelset with the above bounds and extent
  auto levelset = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  // Create the initial substrate surface
  NumericType point[D] = {0.};
  point[0] = origin[0];
  point[1] = origin[1];
  if constexpr (D == 3)
    point[2] = origin[2];

  NumericType normal[D] = {0.};
  normal[D - 1] = 1.;

  lsMakeGeometry<NumericType, D>(
      levelset, lsSmartPointer<lsPlane<NumericType, D>>::New(point, normal))
      .apply();

  // Create the actual trench
  auto cutout =
      lsSmartPointer<lsDomain<NumericType, D>>::New(levelset->getGrid());

  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  const NumericType leftOffset =
      std::tan(leftTaperingAngle * PI / 180.) * trenchDepth;
  const NumericType rightOffset =
      std::tan(rightTaperingAngle * PI / 180.) * trenchDepth;
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
    lsVTKWriter<NumericType>(mesh, "points.vtp").apply();
    lsFromSurfaceMesh<NumericType, D>(cutout, mesh).apply();
  } else {
    for (int i = 0; i < 8; i++) {
      std::array<NumericType, 3> node = {0., 0., 0.};
      mesh->insertNextNode(node);
    }
    mesh->nodes[0][0] = origin[0] - trenchTopWidth / 2. + leftOffset;
    mesh->nodes[0][1] = origin[1] - yExtent / 2. - gridDelta;
    mesh->nodes[0][2] = origin[2] - trenchDepth;

    mesh->nodes[1][0] = origin[0] + trenchTopWidth / 2. - rightOffset;
    mesh->nodes[1][1] = origin[1] - yExtent / 2. - gridDelta;
    mesh->nodes[1][2] = origin[2] - trenchDepth;

    mesh->nodes[2][0] = origin[0] + trenchTopWidth / 2. - rightOffset;
    mesh->nodes[2][1] = origin[1] + yExtent / 2. + gridDelta;
    mesh->nodes[2][2] = origin[2] - trenchDepth;

    mesh->nodes[3][0] = origin[0] - trenchTopWidth / 2. + leftOffset;
    mesh->nodes[3][1] = origin[1] + yExtent / 2. + gridDelta;
    mesh->nodes[3][2] = origin[2] - trenchDepth;

    mesh->nodes[4][0] = origin[0] - trenchTopWidth / 2.;
    mesh->nodes[4][1] = origin[1] - yExtent / 2. - gridDelta;
    mesh->nodes[4][2] = origin[2];

    mesh->nodes[5][0] = origin[0] + trenchTopWidth / 2.;
    mesh->nodes[5][1] = origin[1] - yExtent / 2. - gridDelta;
    mesh->nodes[5][2] = origin[2];

    mesh->nodes[6][0] = origin[0] + trenchTopWidth / 2.;
    mesh->nodes[6][1] = origin[1] + yExtent / 2. + gridDelta;
    mesh->nodes[6][2] = origin[2];

    mesh->nodes[7][0] = origin[0] - trenchTopWidth / 2.;
    mesh->nodes[7][1] = origin[1] + yExtent / 2. + gridDelta;
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

    lsFromSurfaceMesh<NumericType, D>(cutout, mesh).apply();
  }

  lsBooleanOperation<NumericType, D>(
      levelset, cutout, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  return levelset;
}

#endif