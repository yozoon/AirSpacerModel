#include <iostream>
#include <unordered_map>
#include <vector>

#include "CSVReader.hpp"
#include "Utils.hpp"
#include "interpolation/SplineGridInterpolation.hpp"
#include "models/MakeTrenchStamp.hpp"

template <typename NumericType> struct Parameters {
  NumericType gridDelta = 0.2;
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
  static constexpr int numberOfSamples = 512;

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
      std::cerr << "Empty config provided" << std::endl;
      return -1;
    }
    params.fromMap(config);
  }

  NumericType trenchTopWidth =
      (2 * initialTrenchDepth - initialTrenchTopWidth) /
      (2 * params.aspectRatio - 1);

  NumericType trenchDepth = params.aspectRatio * trenchTopWidth;

  std::cout << "trenchTopWidth=" << trenchTopWidth
            << ", trenchDepth=" << trenchDepth << std::endl;

  CSVReader<NumericType> reader;
  reader.setFilename(dataFile);

  // Get a copy of the data from the data source
  auto data = reader.getData();

  // Also read the data that was stored alongside the actual data describing
  // where it was sampled from (at which relative height)
  auto sampleLocations = reader.getPositionalParameters();

  // The input dimension is provided by the csv file named parameter
  // 'InputDim' In our case it is 2 since we use taper angle and sticking
  // probability as input parameters.
  int inputDim;
  auto namedParams = reader.getNamedParameters();
  if (auto id = namedParams.find("InputDimension"); id != namedParams.end()) {
    inputDim = static_cast<int>(std::round(id->second));
  } else {
    std::cout << "'InputDimension' not found in provided data CSV file.\n";
    return EXIT_FAILURE;
  }

  std::cout << "InputDimension: " << inputDim << std::endl;

  // The dimension of the data that is to be interpolated. In our case this
  // are the extracted dimensions at different timesteps (the timesteps are
  // also included in the data itself)
  int numFeatures = data->at(0).size() - inputDim;

  SplineGridInterpolation<NumericType> gridInterpolation;
  gridInterpolation.setDataDimensions(inputDim, numFeatures);
  gridInterpolation.setData(data);
  gridInterpolation.setBCType(SplineBoundaryConditionType::NOT_A_KNOT);

  std::vector<NumericType> evaluationPoint = {params.aspectRatio,
                                              params.leftTaperAngle,
                                              params.stickingProbability, 1.0};

  auto estimationOpt = gridInterpolation.estimate(evaluationPoint);
  if (!estimationOpt) {
    lsMessage::getInstance().addError("Value estimation failed.").print();
    return EXIT_FAILURE;
  }

  auto [estimatedFeatures, isInside] = estimationOpt.value();

  lsMessage::getInstance()
      .addDebug(std::string("Evaluation point within data grid boundaries: ") +
                (isInside ? std::string("true") : std::string("false")))
      .print();
}