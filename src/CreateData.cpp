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

template <typename NumericType> class Parameters {
  const std::string filename = "data.csv";
  const std::vector<NumericType> aspectRatios;
  const std::vector<NumericType> taperAngles;
  const std::vector<NumericType> stickingProbabilities;
  const bool symmetrical = true;

  // Private constructor, such that the parameters object can only be
  // constructed using the `fromMap` function.
  Parameters(std::string &&passedFilename,
             std::vector<NumericType> &&passedAspectRatios,
             std::vector<NumericType> &&passedTaperAngles,
             std::vector<NumericType> &&passedStickingProbabilities,
             bool passedSymmetrical)
      : filename(std::move(passedFilename)),
        aspectRatios(std::move(passedAspectRatios)),
        taperAngles(std::move(passedTaperAngles)),
        stickingProbabilities(std::move(passedStickingProbabilities)),
        symmetrical(passedSymmetrical) {}

public:
  // Factory pattern to construct an instance of the class from a map
  static Parameters
  fromMap(const std::unordered_map<std::string, std::string> &map) {
    // Local variable instances
    std::string filename = "data.csv";
    std::vector<NumericType> aspectRatios;
    std::vector<NumericType> taperAngles;
    std::vector<NumericType> stickingProbabilities;
    bool symmetrical = true;

    // Now assign the items
    Utils::AssignItems(
        map,
        Utils::Item{"filename", filename,
                    [](const std::string &s) { return s; }},
        Utils::Item{"taperAngles", taperAngles, &Utils::toVector<NumericType>},
        Utils::Item{"aspectRatios", aspectRatios,
                    [](const std::string &s) -> std::vector<NumericType> {
                      return Utils::toVector<
                          NumericType,
                          decltype(&Utils::toStrictlyPositive<NumericType>)>(
                          s, &Utils::toStrictlyPositive<NumericType>);
                    }},
        Utils::Item{
            "stickingProbabilities", stickingProbabilities,
            [](const std::string &s) -> std::vector<NumericType> {
              return Utils::toVector<
                  NumericType, decltype(&Utils::toUnitRange<NumericType>)>(
                  s, &Utils::toUnitRange<NumericType>);
            }},
        Utils::Item{
            "symmetrical",
            symmetrical,
            Utils::toBool,
        });

    return Parameters(std::move(filename), std::move(aspectRatios),
                      std::move(taperAngles), std::move(stickingProbabilities),
                      symmetrical);
  }

  std::size_t getNumberOfCombinations() const {
    return aspectRatios.size() * taperAngles.size() *
           stickingProbabilities.size();
  }

  std::string getFilename() const { return filename; }

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
    std::cout << "- symmetrical: " << ((symmetrical) ? "true\n" : "false\n");
    std::cout << "  -> total number of combinations: "
              << getNumberOfCombinations() << '\n';
  }

  struct const_iterator {
    using iter_t = typename std::vector<NumericType>::const_iterator;
    iter_t aspectRatioIterator;
    iter_t taperAngleIterator;
    iter_t stickingProbabilityIterator;

    const std::vector<NumericType> &aspectRatios;
    const std::vector<NumericType> &taperAngles;
    const std::vector<NumericType> &stickingProbabilities;
    const bool symmetrical;

    const_iterator(iter_t passedAspectRatioIterator,
                   iter_t passedTaperAngleIterator,
                   iter_t passedStickingProbabilityIterator,
                   const std::vector<NumericType> &passedAspectRatios,
                   const std::vector<NumericType> &passedTaperAngles,
                   const std::vector<NumericType> &passedStickingProbabilities,
                   const bool passedSymmetrical)
        : aspectRatioIterator(passedAspectRatioIterator),
          taperAngleIterator(passedTaperAngleIterator),
          stickingProbabilityIterator(passedStickingProbabilityIterator),
          aspectRatios(passedAspectRatios), taperAngles(passedTaperAngles),
          stickingProbabilities(passedStickingProbabilities),
          symmetrical(passedSymmetrical) {}

    // Return a tuple so that we can use structural bindings to unpack
    std::tuple<NumericType, NumericType, NumericType, NumericType>
    operator*() const {
      return {*aspectRatioIterator, *taperAngleIterator,
              (symmetrical ? *taperAngleIterator : 0.0),
              *stickingProbabilityIterator};
    };

    const_iterator &operator++() {
      if (++stickingProbabilityIterator == stickingProbabilities.cend()) {
        if (++taperAngleIterator == taperAngles.cend()) {
          if (++aspectRatioIterator == aspectRatios.cend()) {
            return *this;
          }
          taperAngleIterator = taperAngles.cbegin();
        }
        stickingProbabilityIterator = stickingProbabilities.cbegin();
      }
      return *this;
    };

    bool operator!=(const const_iterator &other) const {
      return (aspectRatioIterator != other.aspectRatioIterator) ||
             (taperAngleIterator != other.taperAngleIterator) ||
             (stickingProbabilityIterator != other.stickingProbabilityIterator);
    }
  };

  const_iterator begin() const {
    return const_iterator(aspectRatios.cbegin(), taperAngles.cbegin(),
                          stickingProbabilities.cbegin(), aspectRatios,
                          taperAngles, stickingProbabilities, symmetrical);
  }

  const_iterator end() const {
    return const_iterator(aspectRatios.cend(), taperAngles.cend(),
                          stickingProbabilities.cend(), aspectRatios,
                          taperAngles, stickingProbabilities, symmetrical);
  }
};

template <typename NumericType>
std::string createHeader(const int numberOfSamples,
                         const std::vector<NumericType> &sampleLocations,
                         const int inputDimension) {
  std::string header =
      "aspectRatio,leftTaperAngle,stickingProbability,normalizedTime";
  for (int i = 0; i < numberOfSamples; ++i) {
    header += ",diameter_" + std::to_string(i);
  }
  header += "\nDiameters are normalized to the trench width.\nSample positions "
            "along vertical axis (normalized to trench depth + trench width):";
  header += "\n!" + Utils::join(sampleLocations.begin(), sampleLocations.end());
  header += "\n!InputDimension=" + std::to_string(inputDimension);
  return header;
}

int main(const int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;

  std::array<NumericType, 3> origin{0.};
  NumericType gridDelta = 0.2;
  NumericType trenchTopWidth = 4.; // -> 20 grid points

  int numberOfSamples = 512;

  std::unordered_map<std::string, std::string> config;

  if (argc > 1) {
    config = Utils::readConfigFile(argv[1]);
    if (config.empty()) {
      lsMessage::getInstance().addError("Empty config provided").print();
      return EXIT_FAILURE;
    }
  }
  auto params = Parameters<NumericType>::fromMap(config);
  params.print();

  for (const auto [aspectRatio, leftTaperAngle, rightTaperAngle,
                   stickingProbability] : params) {
    std::cout << aspectRatio << ", " << leftTaperAngle << ", "
              << rightTaperAngle << ", " << stickingProbability << '\n';
  }
  return 0;

  // The data we are going to store consists of stickingProbability,
  // taperAngle, time and the sampled geometry descriptors as provided by the
  // feature extraction.
  int inputDimension = 4;

  // Instantiate the featureExtraction
  auto featureExtraction =
      lsSmartPointer<FeatureExtraction<NumericType, D>>::New();
  featureExtraction->setNumberOfSamples(numberOfSamples,
                                        false /* don't include endpoints */);
  // Uniformly spaced sample points
  featureExtraction->setEdgeAffinity(0.0);
  featureExtraction->initializeSampleLocations();

  // The locations at which the features are extracted (normalized to
  // the trench depth at each timestep)
  auto sampleLocations = featureExtraction->getSampleLocations();

  auto writer =
      lsSmartPointer<psCSVWriter<NumericType>>::New(params.getFilename());

  // Creation of a descriptive header
  auto header = createHeader(numberOfSamples, *sampleLocations, inputDimension);
  writer->setHeader(header);
  writer->initialize();

  unsigned count = 1;
  for (const auto [aspectRatio, leftTaperAngle, rightTaperAngle,
                   stickingProbability] : params) {
    // Calculate the depth of the trench using the aspect ratio
    NumericType trenchDepth = trenchTopWidth * aspectRatio;

    NumericType xExtent = 2.0 * trenchTopWidth;
    // Now that we know all geometry parameters generate the geometry
    auto trench = makeTrench<NumericType, D>(
        gridDelta, xExtent, 10., origin, trenchTopWidth, trenchDepth,
        leftTaperAngle, rightTaperAngle, false /* no periodic boundary*/);

    std::cout << count++ << '/' << params.getNumberOfCombinations()
              << std::endl;

    // Normalize the time scale to the sticking probability, so that we
    // get the same top layer thickness for different sticking
    // probabilities.
    NumericType timeScale = 1.0 / stickingProbability;

    // Ensure that we always have 11 samples
    NumericType maxTime = 10.0;
    NumericType extractionInterval = trenchTopWidth / maxTime;

    featureExtraction->setTrenchDimensions(trenchDepth, trenchTopWidth);

    auto advectionCallback = lsSmartPointer<AdvectionCallback<
        NumericType, D, decltype(featureExtraction)::element_type,
        decltype(writer)::element_type>>::New(timeScale, extractionInterval);

    advectionCallback->setFeatureExtraction(featureExtraction);
    advectionCallback->setPrefixData(std::vector<NumericType>{
        aspectRatio, leftTaperAngle, stickingProbability});
    advectionCallback->setWriter(writer);
    advectionCallback->setModifiers(1.0 / maxTime, 1.0);

    auto geometry = psSmartPointer<psDomain<NumericType, D>>::New();
    geometry->insertNextLevelSet(trench);

    // copy top layer to capture deposition
    auto depoLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(trench);
    geometry->insertNextLevelSet(depoLayer);

    auto processModel =
        SimpleDeposition<NumericType, D>(
            stickingProbability /* particle sticking probability */,
            1.0 /* particle source power */)
            .getProcessModel();

    processModel->setAdvectionCallback(advectionCallback);

    // Since we assume a top rate of 1, the maximum time it takes for the
    // trench to close is equal to the trench top width (if sticking
    // probability 1, otherwise the time until pinchoff will be even less)
    NumericType processDuration = trenchTopWidth / stickingProbability;

    psProcess<NumericType, D> process;
    process.setDomain(geometry);
    process.setProcessModel(processModel);
    process.setNumberOfRaysPerPoint(2000);
    process.setProcessDuration(processDuration);

    // Run the process
    process.apply();

    writer->flush();
  }

  // return EXIT_SUCCESS;
}