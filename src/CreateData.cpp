#include <algorithm>
#include <iostream>
#include <vector>

#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>

#include <SimpleDeposition.hpp>
#include <psProcess.hpp>

#include "AdvectionCallback.hpp"
#include "FeatureExtraction.hpp"
#include "TrenchGeometry.hpp"
#include "Utils.hpp"

template <typename NumericType> struct Parameters {
  std::string filename = "data.csv";

  std::vector<NumericType> aspectRatios;
  std::vector<NumericType> taperAngles;
  std::vector<NumericType> stickingProbabilities;

  bool taperBothSidewalls = true;

  void fromMap(std::unordered_map<std::string, std::string> &m) {
    static auto strictlyPositive = [](const std::string &s) -> NumericType {
      auto value = Utils::convert<NumericType>(s);
      if (value <= 0.0)
        throw std::invalid_argument("`aspectRatio` must be strictly positive.");
      return value;
    };

    static auto unitRange = [](const std::string &s) -> NumericType {
      auto value = Utils::convert<NumericType>(s);
      if (value > 1.0 || value <= 0.0)
        throw std::invalid_argument(
            "`stickingProbability` must be in the range [1,0).");
      return value;
    };

    static auto toBool = [](const std::string &s) -> bool {
      auto lower = s;
      std::transform(lower.begin(), lower.end(), lower.begin(),
                     [](unsigned char c) { return std::tolower(c); });
      if (lower == "true") {
        return true;
      }
      return false;
    };

    Utils::AssignItems(
        m,
        Utils::Item{"filename", filename,
                    [](const std::string &s) { return s; }},
        Utils::Item{"taperAngles", taperAngles, &Utils::toVector<NumericType>},
        Utils::Item{
            "aspectRatios", aspectRatios,
            [](const std::string &s) -> std::vector<NumericType> {
              return Utils::toVector<NumericType, decltype(strictlyPositive)>(
                  s, strictlyPositive);
            }},
        Utils::Item{"stickingProbabilities", stickingProbabilities,
                    [](const std::string &s) -> std::vector<NumericType> {
                      return Utils::toVector<NumericType, decltype(unitRange)>(
                          s, unitRange);
                    }},
        Utils::Item{
            "taperBothSidewalls",
            taperBothSidewalls,
            toBool,
        });
  }

  std::size_t getTotalCombinations() const {
    return aspectRatios.size() * taperAngles.size() *
           stickingProbabilities.size();
  }

  void print() const {
    std::cout << "Parameters: " << std::endl;
    std::cout << "- filename: '" << filename << "'\n";
    std::cout << "- aspectRatios: ["
              << Utils::join(aspectRatios.begin(), aspectRatios.end()) << "]\n";
    std::cout << "- taperAngles: ["
              << join(taperAngles.begin(), taperAngles.end()) << "]\n";
    std::cout << "- stickingProbabilities: ["
              << Utils::join(stickingProbabilities.begin(),
                             stickingProbabilities.end())
              << "]\n";
    std::cout << "- taperBothSidewalls: "
              << ((taperBothSidewalls) ? "true\n" : "false\n");
    std::cout << " -> Total number of unique configurations: "
              << getTotalCombinations() << std::endl;
  }
};

int main(const int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;

  std::array<NumericType, 3> origin{0.};
  NumericType gridDelta = 0.2;
  NumericType trenchTopWidth = 4.; // -> 20 grid points

  int numberOfSamples = 512;

  Parameters<NumericType> params;
  if (argc > 1) {
    auto config = Utils::readConfigFile(argv[1]);
    if (config.empty()) {
      lsMessage::getInstance().addError("Empty config provided").print();
      return EXIT_FAILURE;
    }
    params.fromMap(config);
  }
  params.print();
  return 0;

  // // The data we are going to store consists of stickingProbability,
  // // taperAngle, time and the sampled geometry descriptors as provided by the
  // // feature extraction.
  // int inputDimension = 4;

  // // Instantiate the featureExtraction
  // auto featureExtraction =
  //     lsSmartPointer<FeatureExtraction<NumericType, D>>::New();
  // featureExtraction->setNumberOfSamples(numberOfSamples,
  //                                       false /* don't include endpoints */);
  // // Uniformly spaced sample points
  // featureExtraction->setEdgeAffinity(0.0);
  // featureExtraction->initializeSampleLocations();

  // // The locations at which the features are extracted (normalized to
  // // the trench depth at each timestep)
  // auto sampleLocations = featureExtraction->getSampleLocations();

  // auto writer =
  // lsSmartPointer<psCSVWriter<NumericType>>::New(params.filename);

  // // Creation of a descriptive header
  // std::string header =
  //     "aspectRatio,leftTaperAngle,stickingProbability,normalizedTime";
  // for (int i = 0; i < numberOfSamples; ++i) {
  //   header += ",diameter_" + std::to_string(i);
  // }
  // header += "\nDiameters are normalized to the trench width.\nSample
  // positions "
  //           "along vertical axis (normalized to trench depth + trench
  //           width):";
  // header +=
  //     "\n!" + Utils::join(sampleLocations->begin(), sampleLocations->end());
  // header += "\n!InputDimension=" + std::to_string(inputDimension);

  // writer->setHeader(header);
  // writer->initialize();

  // unsigned count = 1;
  // for (auto aspectRatio : params.aspectRatios) {
  //   // Calculate the depth of the trench using the aspect ratio
  //   NumericType trenchDepth = trenchTopWidth * aspectRatio;
  //   for (auto leftTaperAngle : params.leftTaperAngles) {
  //     // for (auto rightTaperAngle : rightTaperAngles) {
  //     NumericType rightTaperAngle = 0.0;

  //     NumericType xExtent = 2.0 * trenchTopWidth;

  //     // Now that we know all geometry parameters generate the geometry
  //     auto trench = makeTrench<NumericType, D>(
  //         gridDelta, xExtent, 10., origin, trenchTopWidth, trenchDepth,
  //         leftTaperAngle, rightTaperAngle, false /* no periodic boundary*/);

  //     for (auto stickingProbability : params.stickingProbabilities) {
  //       std::cout << count << '/' << params.getTotalCombinations() <<
  //       std::endl;
  //       ++count;

  //       // Normalize the time scale to the sticking probability, so that we
  //       // get the same top layer thickness for different sticking
  //       // probabilities.
  //       NumericType timeScale = 1.0 / stickingProbability;

  //       // Ensure that we always have 11 samples
  //       NumericType maxTime = 10.0;
  //       NumericType extractionInterval = trenchTopWidth / maxTime;

  //       featureExtraction->setTrenchDimensions(trenchDepth, trenchTopWidth);

  //       auto advectionCallback = lsSmartPointer<AdvectionCallback<
  //           NumericType, D, decltype(featureExtraction)::element_type,
  //           decltype(writer)::element_type>>::New(timeScale,
  //                                                 extractionInterval);

  //       advectionCallback->setFeatureExtraction(featureExtraction);
  //       advectionCallback->setPrefixData(std::vector<NumericType>{
  //           aspectRatio, leftTaperAngle, stickingProbability});
  //       advectionCallback->setWriter(writer);
  //       advectionCallback->setModifiers(1.0 / maxTime, 1.0);

  //       auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
  //       geometry->insertNextLevelSet(trench);

  //       // copy top layer to capture deposition
  //       auto depoLayer = lsSmartPointer<lsDomain<NumericType,
  //       D>>::New(trench); geometry->insertNextLevelSet(depoLayer);

  //       auto processModel =
  //           SimpleDeposition<NumericType, D>(
  //               stickingProbability /* particle sticking probability */,
  //               1.0 /* particle source power */)
  //               .getProcessModel();

  //       processModel->setAdvectionCallback(advectionCallback);

  //       // Since we assume a top rate of 1, the maximum time it takes for the
  //       // trench to close is equal to the trench top width (if sticking
  //       // probability 1, otherwise the time until pinchoff will even be
  //       less) NumericType processDuration = trenchTopWidth /
  //       stickingProbability;

  //       psProcess<NumericType, D> process;
  //       process.setDomain(geometry);
  //       process.setProcessModel(processModel);
  //       process.setNumberOfRaysPerPoint(2000);
  //       process.setProcessDuration(processDuration);

  //       // Run the process
  //       process.apply();

  //       writer->flush();
  //       // }
  //     }
  //   }
  // }

  // return EXIT_SUCCESS;
}