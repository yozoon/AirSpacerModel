#ifndef GEOMETRIC_AIR_GAP_MODEL_HPP
#define GEOMETRIC_AIR_GAP_MODEL_HPP

#include <array>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include "../Utils.hpp"
#include "../interpolation/SplineGridInterpolation.hpp"
#include "MakeTrenchStamp.hpp"

template <typename NumericType> struct GeometricAirGapModelParameters {
  std::array<NumericType, 3> origin;
  NumericType leftTaperAngle;
  NumericType rightTaperAngle;
};

template <typename NumericType, int D> class GeometricAirGapModel {
  using LSPtrType = lsSmartPointer<lsDomain<NumericType, D>>;

  LSPtrType levelset = nullptr;

  std::vector<NumericType> sampleLocations;
  std::array<NumericType, 3> origin;

  NumericType initialHeight;
  NumericType initialWidth;
  NumericType initialTaperAngle;

  SplineGridInterpolation<NumericType> gridInterpolation;

  int horizontalDir = 0;
  int verticalDir = D - 1;
  int trenchDir = D - 2;
  bool initialized = false;

public:
  GeometricAirGapModel(LSPtrType passedLevelset,
                       const std::vector<NumericType> &passedSampleLocations,
                       const std::array<NumericType, 3> &passedOrigin,
                       const NumericType passedInitialHeight,
                       const NumericType passedInitialWidth,
                       const NumericType passedInitialTaperAngle = 0.0)
      : levelset(passedLevelset), sampleLocations(passedSampleLocations),
        origin(passedOrigin), initialHeight(passedInitialHeight),
        initialWidth(passedInitialWidth),
        initialTaperAngle(passedInitialTaperAngle) {}

  void setLevelset(LSPtrType passedLevelset) { levelset = passedLevelset; }

  void
  setSampleLocations(const std::vector<NumericType> &passedSampleLocations) {
    sampleLocations = passedSampleLocations;
  }

  void setOrigin(const std::array<NumericType, 3> &passedOrigin) {
    origin = passedOrigin;
  }

  void setData(lsSmartPointer<const std::vector<std::vector<NumericType>>> data,
               int inputDim) {
    if (!data) {
      lsMessage::getInstance().addError("No data provided!").print();
      return;
    }
    if (data->empty()) {
      lsMessage::getInstance().addError("Provided data is empty!").print();
      return;
    }

    auto numFeatures = data->at(0).size() - inputDim;

    if (numFeatures < 1) {
      lsMessage::getInstance()
          .addError("Invalid input data dimension!")
          .print();
      return;
    }

    gridInterpolation.setDataDimensions(inputDim, numFeatures);
    gridInterpolation.setData(data);
    gridInterpolation.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    gridInterpolation.initialize();

    initialized = true;
  }

  void apply(const NumericType aspectRatio,
             const NumericType stickingProbability,
             const NumericType taperAngle = 0.0) {
    if (!levelset) {
      lsMessage::getInstance().addError("No levelset provided.").print();
      return;
    }

    if (!initialized) {
      lsMessage::getInstance()
          .addError("Interpolation does not have up-to-date data. Call "
                    "`setData` first.")
          .print();
    }

    std::vector<NumericType> evaluationPoint = {aspectRatio,
                                                stickingProbability};

    auto estimationOpt = gridInterpolation.estimate(evaluationPoint);
    if (!estimationOpt) {
      lsMessage::getInstance().addError("Value estimation failed.").print();
      return;
    }

    auto [estimatedFeatures, isInside] = estimationOpt.value();

    lsMessage::getInstance()
        .addDebug(
            std::string("Evaluation point within data grid boundaries: ") +
            (isInside ? std::string("true") : std::string("false")))
        .print();

    if (sampleLocations.size() != estimatedFeatures.size()) {
      lsMessage::getInstance()
          .addError("Mismatch of feature dimensions!")
          .print();
      return;
    }

    std::array<NumericType, 3> adjustedOrigin;
    std::copy(origin.begin(), origin.end(), adjustedOrigin.begin());
    adjustedOrigin[verticalDir] += processTime;

    // First: Fill the volume of the box such that we can later subtract the
    // trench stamp from it.
    NumericType minCorner[D] = {0.};
    minCorner[verticalDir] = origin[verticalDir] - initialHeight;
    minCorner[horizontalDir] = origin[horizontalDir] - initialWidth / 2;

    NumericType maxCorner[D] = {0.};
    maxCorner[verticalDir] = origin[verticalDir];
    maxCorner[horizontalDir] = origin[horizontalDir] + initialWidth / 2;

    if constexpr (D == 3) {
      const auto bounds =
          psUtils::getBoundsFromGrid<NumericType, D>(levelset->getGrid());
      minCorner[trenchDir] = bounds[2 * trenchDir];
      maxCorner[trenchDir] = bounds[2 * trenchDir + 1];
    }

    auto boxFill =
        lsSmartPointer<lsDomain<NumericType, D>>::New(levelset->getGrid());
    auto box = lsSmartPointer<lsBox<NumericType, D>>::New(minCorner, maxCorner);
    lsMakeGeometry<NumericType, D>(boxFill, box).apply();
    lsBooleanOperation<NumericType, D>(levelset, boxFill,
                                       lsBooleanOperationEnum::UNION)
        .apply();

    // Second: Generate the trench stamp based on the interpolated features
    auto stamp = MakeTrenchStamp(levelset->getGrid(), adjustedOrigin,
                                 sampleLocations, estimatedFeatures);
    if (!stamp) {
      lsMessage::getInstance()
          .addError("Error while generating the trench stamp.")
          .print();
      return;
    }

    // Only use the part of the stamp that is within the filled area (don't need
    // any features above the origin, since we are going to remove them by CMP
    // anyhow).
    lsBooleanOperation<NumericType, D>(stamp, boxFill,
                                       lsBooleanOperationEnum::INTERSECT)
        .appy();

    // Third: stamp out the trench
    lsBooleanOperation<NumericType, D>(
        levelset, stamp, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
};
#endif