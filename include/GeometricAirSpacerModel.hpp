#ifndef GEOMETRIC_AIR_SPACER_MODEL_HPP
#define GEOMETRIC_AIR_SPACER_MODEL_HPP

#include <array>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include "CSVReader.hpp"
#include "CreateStampFromFeatures.hpp"
#include "SplineGridInterpolation.hpp"
#include "TrenchGeometry.hpp"
#include "Utils.hpp"

template <typename NumericType, int D> class GeometricAirSpacerModel {
  using LSPtrType = lsSmartPointer<lsDomain<NumericType, D>>;
  using GeometryPtrType = lsSmartPointer<std::vector<LSPtrType>>;

  GeometryPtrType geometry = nullptr;
  std::string dataFilename;
  std::array<NumericType, 3> origin;
  NumericType initialTrenchDepth;
  NumericType initialTrenchTopWidth;

  NumericType aspectRatio;
  NumericType leftTaperAngle;
  NumericType stickingProbability;

public:
  GeometricAirSpacerModel(GeometryPtrType passedGeometry,
                          const std::string &passedDataFilename,
                          const std::array<NumericType, 3> &passedOrigin,
                          const NumericType passedInitialTrenchDepth,
                          const NumericType passedInitialTrenchTopWidth,
                          const NumericType passedAspectRatio,
                          const NumericType passedLeftTaperAngle,
                          const NumericType passedStickingProbability)
      : geometry(passedGeometry), dataFilename(passedDataFilename),
        origin(passedOrigin), initialTrenchDepth(passedInitialTrenchDepth),
        initialTrenchTopWidth(passedInitialTrenchTopWidth),
        aspectRatio(passedAspectRatio), leftTaperAngle(passedLeftTaperAngle),
        stickingProbability(passedStickingProbability) {}

  void apply() {
    if (!geometry) {
      lsMessage::getInstance().addError("No geometry provided").print();
      return;
    }
    if (geometry->size() < 1) {
      lsMessage::getInstance().addError("Empty geometry provided").print();
      return;
    }

    const auto &grid = geometry->back()->getGrid();
    const auto gridDelta = grid.getGridDelta();

    // Load the data captured from previous simulations
    CSVReader<NumericType> reader;
    reader.setFilename(dataFilename);
    auto data = reader.getData();

    // Also read the data that was stored alongside the actual data describing
    // where it was sampled from (at which relative height)
    auto sampleLocations = reader.getPositionalParameters();

    // The input dimension is provided by the csv file named parameter
    // `InputDim` In our case it is 4 since we use aspectRatio, leftTaperAngle,
    // stickintProbability and timeStep as input parameters.
    int inputDim;
    auto namedParams = reader.getNamedParameters();
    if (auto id = namedParams.find("InputDimension"); id != namedParams.end()) {
      inputDim = static_cast<int>(std::round(id->second));
    } else {
      lsMessage::getInstance()
          .addError("`InputDimension` not found in provided data CSV file.")
          .print();
      return;
    }

    // The dimension of the data that is to be interpolated. In our case this
    // are the extracted dimensions at different timesteps (the timesteps are
    // also included in the data itself)
    int numFeatures = data->at(0).size() - inputDim;

    // Instantiate the Spline grid interpolation class
    SplineGridInterpolation<NumericType> gridInterpolation;
    gridInterpolation.setDataDimensions(inputDim, numFeatures);
    gridInterpolation.setData(data);
    gridInterpolation.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
    if (!gridInterpolation.initialize()) {
      lsMessage::getInstance()
          .addError("Initialization of the spline grid interpolator failed.")
          .print();
      return;
    }

    // Use maximum timestep, since this ensures that the trench is closed (but
    // it does not ensure that the trench is closed below the origin)
    const auto uniqueValues = gridInterpolation.getUniqueValues();
    auto maxTime = *(uniqueValues.back().rbegin());
    std::vector<NumericType> evaluationPoint = {aspectRatio, leftTaperAngle,
                                                stickingProbability, maxTime};

    auto estimationOpt = gridInterpolation.estimate(evaluationPoint);
    if (!estimationOpt) {
      lsMessage::getInstance().addError("Value estimation failed.").print();
      return;
    }

    auto [estimatedFeatures, isInside] = estimationOpt.value();

#ifndef NDEBUG
    lsMessage::getInstance()
        .addDebug(
            std::string("Evaluation point within data grid boundaries: ") +
            (isInside ? std::string("true") : std::string("false")))
        .print();
#endif

    // Calculate the actual trench dimensions using the assumption that there is
    // a conformal coating step before the non-conformal one.
    NumericType trenchTopWidth =
        (2 * initialTrenchDepth - initialTrenchTopWidth) /
        (2 * aspectRatio - 1);
    NumericType trenchDepth = aspectRatio * trenchTopWidth;

#ifndef NDEBUG
    NumericType conformalLayerThickness = initialTrenchDepth - trenchDepth;
    lsMessage::getInstance()
        .addDebug("trenchTopWidth=" + std::to_string(trenchTopWidth) +
                  ", trenchDepth=" + std::to_string(trenchDepth) +
                  ", conformalLayerThickness=" +
                  std::to_string(conformalLayerThickness))
        .print();
#endif

    const bool fixTopWidth = true;

    auto outerTrench = createTrenchStamp<NumericType, D>(
        grid, origin, initialTrenchDepth, initialTrenchTopWidth, leftTaperAngle,
        0.0, fixTopWidth);

    auto innerTrench = createTrenchStamp<NumericType, D>(
        grid, origin, trenchDepth, trenchTopWidth, leftTaperAngle, 0.0,
        fixTopWidth);

    if (!outerTrench || !innerTrench)
      return;

    // Copy the substrate into two new layers
    auto conformalLayer =
        lsSmartPointer<lsDomain<NumericType, D>>::New(geometry->back());

    auto nonConformalLayer =
        lsSmartPointer<lsDomain<NumericType, D>>::New(geometry->back());

    // Cut the outer trench out of the substrate
    lsBooleanOperation<NumericType, D>(
        geometry->back(), outerTrench,
        lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

#ifndef NDEBUG
    Utils::printSurface(geometry->back(), "trench.vtp");
#endif

    // Cut the inner trench out of the conformal layer
    lsBooleanOperation<NumericType, D>(
        conformalLayer, innerTrench,
        lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    geometry->push_back(conformalLayer);

    trenchTopWidth = trenchDepth / aspectRatio;

    auto stamp = createStampFromFeatures<NumericType, D>(
        grid, origin, trenchDepth, trenchTopWidth, sampleLocations,
        estimatedFeatures, true /* moveOffset */);

#ifndef NDEBUG
    Utils::printSurface(stamp, "stamp.vtp");
#endif

    // And finally cut the interpolated air spacer out of the non-conformal
    // layer
    lsBooleanOperation<NumericType, D>(
        nonConformalLayer, stamp, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    geometry->push_back(nonConformalLayer);
  }
};
#endif