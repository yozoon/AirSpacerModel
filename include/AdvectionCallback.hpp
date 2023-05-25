#ifndef ADVECTION_CALLBACK_HPP
#define ADVECTION_CALLBACK_HPP

#include <algorithm>
#include <filesystem>
#include <vector>

#include <psAdvectionCallback.hpp>
#include <psCSVWriter.hpp>
#include <psSmartPointer.hpp>

#include <lsWriter.hpp>

namespace fs = std::filesystem;

template <typename NumericType, int D, typename FeatureExtractionType,
          typename WriterType>
class AdvectionCallback : public psAdvectionCallback<NumericType, D> {
private:
  NumericType timeScale = 1.0;
  NumericType extractionInterval = 1.0;

  NumericType processTime = 0.0;
  NumericType lastUpdateTime = 0.0;
  size_t counter = 0;
  std::vector<NumericType> prefixData;
  NumericType timeModifier = 1.0;
  NumericType lengthModifier = 1.0;

  psSmartPointer<FeatureExtractionType> featureExtraction = nullptr;
  psSmartPointer<WriterType> writer = nullptr;
  psSmartPointer<std::vector<NumericType>> dataPtr = nullptr;

protected:
  using psAdvectionCallback<NumericType, D>::domain;

public:
  AdvectionCallback() {}

  AdvectionCallback(NumericType passedTimeScale,
                    NumericType passedExtractionInterval)
      : timeScale(passedTimeScale),
        extractionInterval(passedExtractionInterval) {}

  void setFeatureExtraction(
      psSmartPointer<FeatureExtractionType> passedFeatureExtraction) {
    featureExtraction = passedFeatureExtraction;
  }

  void setModifiers(NumericType passedTimeModifier,
                    NumericType passedLengthModifier) {
    timeModifier = passedTimeModifier;
    lengthModifier = passedLengthModifier;
  }

  void setWriter(psSmartPointer<WriterType> passedWriter) {
    writer = passedWriter;
  }

  void setPrefixData(const std::vector<NumericType> &passedPrefixData) {
    prefixData = passedPrefixData;
  }

  void setDataPtr(psSmartPointer<std::vector<NumericType>> passedDataPtr) {
    dataPtr = passedDataPtr;
  }

  void apply() {
    std::ostringstream name;
    for (unsigned i = 0; i < prefixData.size(); ++i)
      name << prefixData[i] << "_";
    name << counter;
    // Utils::printSurface(domain->getLevelSets()->back(), name.str() + ".vtp");
    lsWriter<NumericType, D>(
        domain->getLevelSets()->back(),
        (fs::path("lvst") / (name.str() + ".lvst")).string())
        .apply();
    std::cout << "-- " << counter << '\n';

    featureExtraction->setDomain(domain->getLevelSets()->back());
    featureExtraction->apply();

    auto features = featureExtraction->getFeatures();
    if (features) {
      std::vector<NumericType> row(prefixData.begin(), prefixData.end());
      row.push_back(timeModifier * counter);
      std::transform(features->begin(), features->end(),
                     std::back_inserter(row),
                     [=](auto &v) { return lengthModifier * v; });

      if (writer)
        writer->writeRow(row);
      if (dataPtr)
        std::copy(row.begin(), row.end(), std::back_inserter(*dataPtr));
    }
  }

  bool applyPreAdvect(const NumericType passedProcessTime) override {
    if (passedProcessTime == 0.) {
      if (!featureExtraction)
        return false;
      if (!fs::exists("lvst"))
        fs::create_directories("lvst");
      apply();
      ++counter;
      lastUpdateTime = 0.;
    }
    processTime = passedProcessTime;
    return true;
  }

  bool applyPostAdvect(const NumericType advectionTime) override {
    processTime += advectionTime;
    if (processTime - lastUpdateTime >= extractionInterval * timeScale) {
      apply();
      lastUpdateTime = counter * extractionInterval * timeScale;
      ++counter;
    }
    return true;
  }
};
#endif