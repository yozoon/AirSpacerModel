#include <numeric>
#include <string>

#include <fmt/core.h>

#include <lsMessage.hpp>

#include "CSVReader.hpp"
#include "SplineGridInterpolation.hpp"

int main(const int argc, const char *const *const argv) {
  // Input: sticking probability, initial_bottom_width, desired depth, maximum
  // top width Output: is there a solution at all. If there is a solution:
  // minimum taper angle that ensures trench filling.
  using NumericType = double;
  static constexpr int D = 2;

  std::string dataFilename = "trenchfill.csv";
  if (argc > 1)
    dataFilename = argv[1];

  // Load the data captured from previous simulations
  CSVReader<NumericType> reader;
  reader.setFilename(dataFilename);
  auto data = reader.getData();

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
  }

  // The dimension of the data that is to be interpolated. In our case this
  // are the extracted dimensions at different timesteps (the timesteps are
  // also included in the data itself)
  int numFeatures = data->at(0).size() - inputDim;

  unsigned numFeaturesRight =
      static_cast<unsigned>(std::ceil(1.0 * numFeatures / 2));
  unsigned numFeaturesLeft = numFeatures - numFeaturesRight;

  auto fillingRatioData =
      lsSmartPointer<std::vector<std::vector<NumericType>>>::New();
  fillingRatioData->reserve(data->size());
  for (auto &d : *data) {
    NumericType aspectRatio = d[0];
    NumericType taperAngle = d[1];
    NumericType stickingProbability = d[2];
    NumericType time = d[3];

    auto tan = std::tan(Utils::deg2rad(taperAngle));

    // The percentage of features that are above the initial trench top surface
    NumericType aboveZeroRatio = (1.0 + 2.0 * aspectRatio * tan) /
                                 (1.0 + 2.0 * aspectRatio * tan + aspectRatio);

    int numRightFeaturesAboveZero =
        static_cast<int>(std::ceil(numFeaturesRight * aboveZeroRatio));

    // int numLeftFeaturesAboveZero =
    //     static_cast<int>(std::ceil(numFeaturesLeft * aboveZeroRatio));

    // Percentage of the trench below the zero line that is filled with the
    // non-conformal coating. Count values smaller than 2% of top opening as
    // closed. Since the trench is symmetrical we only look at the features on
    // the right.
    int closedCount = std::count_if(
        std::next(d.begin(), inputDim + numRightFeaturesAboveZero),
        std::next(d.begin(), inputDim + numFeaturesRight),
        [](NumericType p) { return p <= 0.02; });

    int belowZeroCount = numFeaturesRight - numRightFeaturesAboveZero;

    // NumericType cummulativeOpening = 0.;
    // for (unsigned i = numRightFeaturesAboveZero; i < numFeaturesRight; ++i) {
    //   auto diff =
    //       d[inputDim + i] -
    //       d[inputDim + numFeaturesRight + std::min(numFeaturesLeft - 1U, i)];
    //   cummulativeOpening += std::max(0., diff);
    // }

    // auto filledAreaFraction = 1.0 - (1.0 + 2.0 * aspectRatio * tan) /
    //                                     (1.0 + aspectRatio * tan) *
    //                                     cummulativeOpening / belowZeroCount;

    NumericType fillingRatio = 1.0 * closedCount / belowZeroCount;
    fillingRatioData->emplace_back(std::vector<NumericType>{
        aspectRatio, taperAngle, stickingProbability, time, fillingRatio});

    fmt::print("{},{},{:.2f},{},{:.4f}\n", static_cast<int>(aspectRatio),
               static_cast<int>(taperAngle), stickingProbability, time,
               fillingRatio);
  }

  auto first = fillingRatioData->at(0);

  SplineGridInterpolation<NumericType> sgi;
  sgi.setDataDimensions(inputDim, 1);
  sgi.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);
  sgi.setData(fillingRatioData);
  sgi.initialize();

  std::vector<NumericType> loc;
  for (int i = 0; i < inputDim; ++i) {
    loc.push_back(first[i]);
  }
  auto grad = sgi.gradient(loc);
  fmt::print("Done\n");
}