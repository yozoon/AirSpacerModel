#ifndef CREATE_STAMP_FROM_FEATURES_HPP
#define CREATE_STAMP_FROM_FEATURES_HPP

#include <algorithm>
#include <array>
#include <vector>

#include <hrleGrid.hpp>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsSmartPointer.hpp>

#ifndef NDEBUG
#include <lsVTKWriter.hpp>
#endif

#include "Utils.hpp"

template <typename NumericType, int D>
lsSmartPointer<lsDomain<NumericType, D>> createStampFromFeatures(
    const hrleGrid<D> &grid, const std::array<NumericType, 3> &origin,
    const NumericType trenchDepth, const NumericType trenchTopWidth,
    const std::vector<NumericType> &featureLocations,
    const std::vector<NumericType> &features) {
  int verticalDir = D - 1;
  int horizontalDir = 0;
  int trenchDir = D - 2;

  if (features.size() != featureLocations.size() || features.size() < 5) {
    std::cout << "Not enough features provided for reconstruction!\n";
    return nullptr;
  }

  const auto bounds = Utils::getBoundsFromGrid<NumericType, D>(grid);

  // INFO: Mesh point numbering CW: solid is enclosed inside points

  // Now create the stamp for the trench profile (the negative of the trench
  // geometry)
  auto stamp = lsSmartPointer<lsDomain<NumericType, D>>::New(grid);

  NumericType gridDelta = grid.getGridDelta();

  unsigned numSamplesRight =
      static_cast<unsigned>(std::ceil(1.0 * featureLocations.size() / 2));
  unsigned numSamplesLeft = featureLocations.size() - numSamplesRight;

  NumericType reconstructionRange = trenchDepth + trenchTopWidth;
  NumericType trenchBase = -trenchDepth;

#ifndef NDEBUG
  int meshCount = 0;
#endif

  unsigned j = 0;
  while (j < std::min(numSamplesRight, numSamplesLeft)) {
    // Manually create a surface mesh based on the extracted features
    auto mesh = lsSmartPointer<lsMesh<>>::New();

    unsigned k = 0;

    if (j == 0) {
      // Add one point on top of the geometry, in order to avoid potential
      // sharp corners
      std::array<NumericType, 3> point = origin;

      point[horizontalDir] += trenchTopWidth * features.front();
      point[verticalDir] += trenchTopWidth + 2 * gridDelta;
      if constexpr (D == 2) {
        mesh->insertNextNode(point);
      } else if constexpr (D == 3) {
        point[trenchDir] = bounds[2 * trenchDir + 1];
        mesh->insertNextNode(point);

        point[trenchDir] = bounds[2 * trenchDir];
        mesh->insertNextNode(point);
      }
      ++k;
    }

    unsigned nextJ = j;

    // Add points of the right sidewall to the mesh until we encounter a
    // pinch-off (the distance between the left and right sidwall features
    // at that location falling below a certain threshold)
    for (unsigned i = nextJ; i < numSamplesRight; ++i) {
      nextJ = i;

      if (features[i] -
              features[numSamplesRight + std::min(numSamplesLeft - 1, i)] <
          gridDelta / 5) {
#ifndef NDEBUG
        std::cout << "Pinchoff point detected!\n";
#endif
        break;
      }

      std::array<NumericType, 3> point = origin;
      point[horizontalDir] += trenchTopWidth * features[i];
      point[verticalDir] +=
          trenchBase + reconstructionRange * featureLocations[i];

      // std::cout << point[verticalDir] << '\n';
      if constexpr (D == 2) {
        mesh->insertNextNode(point);
      } else if constexpr (D == 3) {
        point[trenchDir] = bounds[2 * trenchDir + 1];
        mesh->insertNextNode(point);

        point[trenchDir] = bounds[2 * trenchDir];
        mesh->insertNextNode(point);
      }
      ++k;
    }

    // Add the points of the left sidewall in reverse order to the mesh,
    // starting from the last, or the one at the previously determined
    // pinch-off location.
    for (int i = std::min(numSamplesLeft - 1, nextJ); i > static_cast<int>(j);
         --i) {
      std::array<NumericType, 3> point = origin;

      point[horizontalDir] += trenchTopWidth * features[numSamplesRight + i];
      point[verticalDir] +=
          trenchBase +
          reconstructionRange * featureLocations[numSamplesRight + i];

      if constexpr (D == 2) {
        mesh->insertNextNode(point);
      } else if constexpr (D == 3) {
        point[trenchDir] = bounds[2 * trenchDir + 1];
        mesh->insertNextNode(point);

        point[trenchDir] = bounds[2 * trenchDir];
        mesh->insertNextNode(point);
      }
      ++k;
    }

    if (j == 0) {
      // Add one point on top of the geometry, in order to avoid potential
      // sharp corners
      std::array<NumericType, 3> point = origin;

      point[horizontalDir] += trenchTopWidth * features.back();
      point[verticalDir] += trenchTopWidth + 2 * gridDelta;
      if constexpr (D == 2) {
        mesh->insertNextNode(point);
      } else if constexpr (D == 3) {
        point[trenchDir] = bounds[2 * trenchDir + 1];
        mesh->insertNextNode(point);

        point[trenchDir] = bounds[2 * trenchDir];
        mesh->insertNextNode(point);
      }
      ++k;
    }

    j = nextJ + 1;

    // If there are not enough points for a triangle skip the boolean
    // operation.
    if (mesh->nodes.size() < 3) {
#ifndef NDEBUG
      std::cout << "Mesh is just a line. Skipping.\n";
#endif
      continue;
    }

    if constexpr (D == 2) {
      // Connect all nodes of the mesh with lines
      for (unsigned i = 0; i < mesh->nodes.size() - 1; ++i)
        mesh->lines.emplace_back(std::array<unsigned, 2>{i, i + 1});
      // Close the hull mesh
      mesh->lines.emplace_back(std::array<unsigned, 2>{k - 1, 0U});
    } else if constexpr (D == 3) {
      // Triangles with two nodes at the front
      for (unsigned i = 0; i < k - 1; ++i) {
        mesh->insertNextTriangle(
            std::array<unsigned, 3>{2 * (i + 1) /* Front */, 2 * i /* Front */,
                                    2 * (i + 1) + 1 /* Back*/});
      }
      mesh->insertNextTriangle(std::array<unsigned, 3>{
          0 /* Front */, 2 * (k - 1) /* Front */, 1 /* Back */});

      // Triangles with two nodes at the back
      for (unsigned i = 0; i < k - 1; ++i) {
        mesh->insertNextTriangle(std::array<unsigned, 3>{
            2 * i + 1 /* Back */, 2 * (i + 1) + 1 /* Back*/,
            2 * i /* Front */});
      }
      mesh->insertNextTriangle(std::array<unsigned, 3>{
          2 * (k - 1) + 1 /* Back */, 1 /* Back*/, 2 * (k - 1) /* Front */});

      // Determine the center point of the cross section
      std::array<NumericType, 3> center{0.};
      std::for_each(mesh->nodes.begin(), mesh->nodes.end(),
                    [&, i = 0](const auto &n) mutable {
                      if (++i % 2 == 1) {
                        center[horizontalDir] += n[horizontalDir];
                        center[verticalDir] += n[verticalDir];
                      }
                    });
      center[horizontalDir] /= k;
      center[verticalDir] /= k;

      // Front cover
      center[trenchDir] = bounds[2 * trenchDir + 1];
      unsigned frontCenterIndex = mesh->insertNextNode(center);
      for (unsigned i = 0; i < k - 1; ++i) {
        mesh->insertNextTriangle(
            std::array<unsigned, 3>{2 * i, 2 * (i + 1), frontCenterIndex});
      }
      mesh->insertNextTriangle(
          std::array<unsigned, 3>{2 * (k - 1), 0, frontCenterIndex});

      // Back cover
      center[trenchDir] = bounds[2 * trenchDir];
      unsigned backCenterIndex = mesh->insertNextNode(center);
      for (unsigned i = 0; i < k - 1; ++i) {
        mesh->insertNextTriangle(std::array<unsigned, 3>{
            2 * (i + 1) + 1, 2 * i + 1, backCenterIndex});
      }
      mesh->insertNextTriangle(
          std::array<unsigned, 3>{1, 2 * (k - 1) + 1, backCenterIndex});
    }

#ifndef NDEBUG
    //  Print the created hull mesh
    lsVTKWriter<NumericType>(mesh,
                             "hullMesh" + std::to_string(meshCount++) + ".vtp")
        .apply();
#endif

    // Create the new levelset based on the mesh and union it with the
    // previously generated levelsets.
    auto hull = lsSmartPointer<lsDomain<NumericType, D>>::New(grid);

    lsFromSurfaceMesh<NumericType, D>(hull, mesh).apply();

    lsBooleanOperation<NumericType, D>(stamp, hull,
                                       lsBooleanOperationEnum::UNION)
        .apply();
  }
  return stamp;
}

#endif