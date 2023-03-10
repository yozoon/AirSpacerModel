#ifndef FEATURE_EXTRACTION_HPP
#define FEATURE_EXTRACTION_HPP

#include <algorithm>
#include <array>
#include <vector>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>

#ifndef NDEBUG
#include <lsVTKWriter.hpp>
#endif

#include <psKDTree.hpp>

template <typename NumericType, int D> class FeatureExtraction {
  using ConstPtr = lsSmartPointer<const std::vector<NumericType>>;
  enum FeatureLabelEnum : unsigned {
    NONE = 0,
    LEFT_SIDEWALL = 1,
    RIGHT_SIDEWALL = 2,
    TRENCH_BOTTOM = 3,
  };

  lsSmartPointer<lsDomain<NumericType, D>> levelset = nullptr;

  std::array<NumericType, 3> origin{0.};
  int numberOfSamples = 32;
  NumericType edgeAffinity = 0.;

  bool includeEndpoints = false;

  int horizontalDir = 0;
  int verticalDir = D - 1;
  int trenchDir = D - 2;

  std::vector<NumericType> sampleLocations;
  std::vector<NumericType> features;
  NumericType trenchDepth;
  NumericType trenchTopWidth;
  unsigned numSamplesRight;

public:
  FeatureExtraction() {}

  FeatureExtraction(lsSmartPointer<lsDomain<NumericType, D>> passedDomain)
      : levelset(passedDomain) {}

  void setDomain(lsSmartPointer<lsDomain<NumericType, D>> passedDomain) {
    levelset = passedDomain;
  }

  void setNumberOfSamples(int passedNumberOfSamples,
                          bool passedIncludeEndpoints = true) {
    assert(numberOfSamples > 1);
    numberOfSamples = passedNumberOfSamples;
    includeEndpoints = passedIncludeEndpoints;
    sampleLocations.clear();
  }

  void setEdgeAffinity(NumericType passedEdgeAffinity) {
    edgeAffinity = passedEdgeAffinity;
    sampleLocations.clear();
  }

  void setOrigin(const std::array<NumericType, 3> &passedOrigin) {
    origin = passedOrigin;
  }

  void setTrenchDimensions(NumericType passedTrenchDepth,
                           NumericType passedTrenchTopWidth) {
    trenchDepth = passedTrenchDepth;
    trenchTopWidth = passedTrenchTopWidth;
  }

  ConstPtr getFeatures() const { return ConstPtr::New(features); }

  ConstPtr getSampleLocations() const { return ConstPtr::New(sampleLocations); }

  void apply() {
    initializeSampleLocations();

    if (!levelset)
      return;

    const NumericType gridDelta = levelset->getGrid().getGridDelta();

    // Re-initialize the feature vector with a value of zero.
    features.clear();
    features.resize(sampleLocations.size(), 0.);

    // Convert the geometry to a surface mesh and extract the nodes as well as
    // the void points
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToDiskMesh<NumericType, D>(levelset, mesh).apply();
    auto nodes = mesh->getNodes();

    auto normalsPtr = mesh->getCellData().getVectorData(
        lsCalculateNormalVectors<NumericType, D>::normalVectorsLabel);
    if (normalsPtr == nullptr) {
      return;
    }

    auto normals = *normalsPtr;

    std::vector<NumericType> featureLabels(nodes.size(), 0.);
    for (unsigned i = 0; i < normals.size(); ++i) {
      const auto &normal = normals[i];
      auto nx = normal[horizontalDir];
      NumericType threshold = 0.05;
      NumericType label = FeatureLabelEnum::NONE;
      if (nx >= threshold) {
        label = FeatureLabelEnum::LEFT_SIDEWALL;
      } else if (nx <= -threshold) {
        label = FeatureLabelEnum::RIGHT_SIDEWALL;
      } else if (nodes[i][verticalDir] < origin[verticalDir] - gridDelta) {
        label = FeatureLabelEnum::TRENCH_BOTTOM;
      }
      featureLabels[i] = label;
    }

#ifndef NDEBUG
    mesh->getPointData().insertNextScalarData(featureLabels, "feature");
    lsVTKWriter<NumericType>(mesh, "FeatureExtractionLabels.vtp").apply();
#endif
    // Project all points onto the hyperplane spanned by all axis except the
    // one of the direction of the trench diameter. (i.e. project on the
    // vertical axis in 2D and in 3D project onto the symmetry plane spanned by
    // the vertical axis and the extrusion axis).
    std::vector<std::vector<NumericType>> vertical;
    vertical.reserve(nodes.size());
    std::transform(
        nodes.begin(), nodes.end(), std::back_inserter(vertical),
        [=](auto &node) {
          if constexpr (D == 3) {
            return std::vector<NumericType>{node[verticalDir], node[trenchDir]};
          } else {
            return std::vector<NumericType>{node[verticalDir]};
          }
        });

    // Initialize a 1D KDTree (~= multiset data structure with convencience
    // functions)
    psKDTree<NumericType> tree;
    tree.setPoints(vertical);
    tree.build();

    NumericType trenchBase = origin[verticalDir] - trenchDepth;
    NumericType extractionRange = trenchDepth + trenchTopWidth;
    NumericType searchRadius = 2.0 * gridDelta;

    // The extract the diameters along its depth at the relative coordinates
    // given by depths
    for (unsigned i = 0; i < sampleLocations.size(); ++i) {
      std::vector<NumericType> loc = {trenchBase +
                                      extractionRange * sampleLocations[i]};
      if constexpr (D == 3) {
        loc.push_back(origin[trenchDir]);
      }

      auto neighborsOpt = tree.findNearestWithinRadius(loc, searchRadius);
      if (!neighborsOpt)
        continue;

      auto neighbors = neighborsOpt.value();

      if (neighbors.empty())
        continue;

      FeatureLabelEnum matchLabel;
      if (i < numSamplesRight) {
        matchLabel = FeatureLabelEnum::RIGHT_SIDEWALL;
      } else {
        matchLabel = FeatureLabelEnum::LEFT_SIDEWALL;
      }

      int upperIdx = -1;
      int lowerIdx = -1;
      NumericType upperDistance = std::numeric_limits<NumericType>::max();
      NumericType lowerDistance = std::numeric_limits<NumericType>::max();
      NumericType upperWidth = 0.;
      NumericType lowerWidth = 0.;
      for (const auto &nb : neighbors) {
        NumericType label = featureLabels[nb.first];
        NumericType nodeZ = nodes[nb.first][verticalDir];
        NumericType nodeX = nodes[nb.first][horizontalDir];
        if (label == matchLabel) {
          if (lowerIdx == -1 && nodeZ < loc[0]) {
            lowerIdx = nb.first;
            lowerDistance = loc[0] - nodeZ;
            lowerWidth = nodeX - origin[horizontalDir];
          } else if (upperIdx == -1 && nodeZ >= loc[0]) {
            upperIdx = nb.first;
            upperDistance = nodeZ - loc[0];
            upperWidth = nodeX - origin[horizontalDir];
          }
        }
        if (upperIdx != -1 && lowerIdx != -1)
          break;
      }

      if (upperIdx != 0 && upperDistance < 1e-4) {
        // If the vertical position of the upper point coincides with the sample
        // position up to a certain epsilon, use its width.
        features[i] = upperWidth / trenchTopWidth;
      } else if (lowerIdx != 0 && lowerDistance < 1e-4) {
        // If the vertical position of the lower point coincides with the sample
        // position up to a certain epsilon, use its width.
        features[i] = lowerWidth / trenchTopWidth;
      } else if (upperIdx != -1 && lowerIdx != -1) {
        // Otherwise linearly interpolate between the two widths based on the
        // offset to to sample position.
        NumericType totalDistance = lowerDistance + upperDistance;
        features[i] =
            (lowerDistance * upperWidth + upperDistance * lowerWidth) /
            totalDistance / trenchTopWidth;
      }
    }
  }

  void initializeSampleLocations() {
    // Sample locations are in the range 0 (bottom) ... 1 (top)
    if (!sampleLocations.empty())
      return;

    // The first sample location is negative - indicating that this is the
    // height measurement
    sampleLocations.clear();
    sampleLocations.resize(numberOfSamples, 0.0);

    // Ensure that the number of samples of the right sidewall is greater than
    // or equal to the number of samples of the left sidewall.
    numSamplesRight =
        static_cast<unsigned>(std::ceil(1.0 * numberOfSamples / 2));

    // The remaining sample locations are distributed in the range 0 to 1.
    // Left sidewall sample locations (top to bottom -> descending)
    distributeSampleLocations(
        sampleLocations.begin(),
        std::next(sampleLocations.begin(), numSamplesRight), edgeAffinity,
        /*not ascending*/ false, includeEndpoints);

    // Right sidewall sample locations (top to bottom -> descending)
    distributeSampleLocations(
        std::next(sampleLocations.begin(), numSamplesRight),
        sampleLocations.end(), edgeAffinity, /*not ascending*/ false,
        includeEndpoints);
  }

private:
  // Populate the given range with a sequence of values in the range from 0
  // to 1 and place them closer to the edges or closer to the center based on
  // the edgeAffinity parameter.
  void
  distributeSampleLocations(typename std::vector<NumericType>::iterator begin,
                            typename std::vector<NumericType>::iterator end,
                            NumericType edgeAffinity = 0.0,
                            bool ascending = true,
                            bool includeEndpoints = true) const {
    auto n = std::distance(begin, end);
    if (n < 1)
      return;

    // Generate n evenly spaced points in the interval [-1, 1] (or (-1, 1) if
    // includeEndpoints==false)
    if (includeEndpoints) {
      std::generate(begin, end, [=, i = 0]() mutable {
        return -1.0 + 2.0 * i++ / (n - 1);
      });
    } else {
      std::generate(begin, end, [=, i = 0]() mutable {
        return -1.0 + 2.0 * (i++ + 0.5) / n;
      });
    }

    NumericType maxVal = 1.0;
    if (edgeAffinity != 0.0) {
      // Spread the points. Increase density around zero (if edgeAffinity < 0)
      // or increase density at the edges (if edgeAffinity > 0)
      std::transform(begin, end, begin, [edgeAffinity](NumericType xi) {
        return xi < 0 ? 1.0 - std::exp(edgeAffinity * xi)
                      : std::expm1(-edgeAffinity * xi);
      });
      maxVal = std::abs(std::expm1(-edgeAffinity));
    }
    // Now transform the points back into the interval [0,1] (or (0,1) if
    // includeEndpoints==False)
    std::transform(begin, end, begin,
                   [=](NumericType xi) { return (xi / maxVal + 1.0) / 2.0; });

    if (ascending) {
      if (edgeAffinity > 0)
        std::transform(begin, end, begin,
                       [](const auto &v) { return 1.0 - v; });
    } else {
      if (edgeAffinity <= 0)
        std::transform(begin, end, begin,
                       [](const auto &v) { return 1.0 - v; });
    }
  }
};
#endif