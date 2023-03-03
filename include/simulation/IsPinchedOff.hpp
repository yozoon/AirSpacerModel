#ifndef PINCHOFF_DETECTION_HPP
#define PINCHOFF_DETECTION_HPP

#include <algorithm>

#include <lsMarkVoidPoints.hpp>
#include <lsMesh.hpp>
#include <lsSmartPointer.hpp>
#include <lsToSurfaceMesh.hpp>

template <typename NumericType, int D>
bool IsPinchedOff(lsSmartPointer<lsDomain<NumericType, D>> levelset,
                  NumericType threshold = 0.2) {
  lsMarkVoidPoints<NumericType, D>(levelset).apply();
  // Generate a surface mesh
  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  lsToSurfaceMesh<NumericType, D>(levelset, mesh).apply();

  // Check if a sufficient percentage of points belong to a void. If that's
  // the case, mark the trench geometry as closed.
  // If at least one fifth of all points belong to a void, mark the
  // geometry as closed. (The threshold avoids errors which occur when the
  // mark void points algorithm marks single points as voids.)
  auto voidPoints = mesh->getPointData().getScalarData(
      lsMarkVoidPoints<NumericType, D>::voidPointLabel);

  return std::count_if(voidPoints->cbegin(), voidPoints->cend(),
                       [](NumericType p) { return p == 1.; }) >
         threshold * mesh->getNodes().size();
}
#endif