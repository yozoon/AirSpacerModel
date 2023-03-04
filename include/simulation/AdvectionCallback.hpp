#ifndef ADVECTION_CALLBACK_HPP
#define ADVECTION_CALLBACK_HPP

#include <algorithm>
#include <vector>

#include <psAdvectionCallback.hpp>
#include <psCSVWriter.hpp>
#include <psSmartPointer.hpp>

#include <lsWriter.hpp>

// #include "Utils.hpp"

template <typename NumericType, int D, typename FeatureExtractionType,
          typename WriterType>
class AdvectionCallback : public psAdvectionCalback<NumericType, D> {
private:
  NumericType timeScale = 1.0;
  NumericType extractionInterval = 1.0;

  NumericType processTime = 0.0;
  NumericType lastUpdateTime = 0.0;
  size_t counter = 0;
  std::vector<NumericType> prefixData;

  psSmartPointer<FeatureExtractionType> featureExtraction = nullptr;
  psSmartPointer<WriterType> writer = nullptr;
  psSmartPointer<std::vector<NumericType>> dataPtr = nullptr;

protected:
  using psAdvectionCalback<NumericType, D>::domain;

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
    lsWriter<NumericType, D>(domain->getLevelSets()->back(),
                             name.str() + ".lvst")
        .apply();
    std::cout << "-- " << counter << '\n';

    featureExtraction->setDomain(domain->getLevelSets()->back());
    featureExtraction->apply();

    auto features = featureExtraction->getFeatures();
    if (features) {
      std::vector<NumericType> row(prefixData.begin(), prefixData.end());
      row.push_back(counter);
      std::copy(features->begin(), features->end(), std::back_inserter(row));

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