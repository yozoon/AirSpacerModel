#ifndef CREATE_DATA_HPP
#define CREATE_DATA_HPP

#include <algorithm>
#include <iostream>
#include <vector>

#include "Utils.hpp"

// Data creation parameters with const iterator interface
template <typename NumericType> class Parameters {
public:
  const std::string filename = "data.csv";
  const std::vector<NumericType> aspectRatios;
  const std::vector<NumericType> taperAngles;
  const std::vector<NumericType> stickingProbabilities;
  const std::size_t timeSteps;
  const bool symmetrical;
  const bool fixTopWidth;

private:
  // Private constructor, such that the parameters object can only be
  // constructed using the `fromMap` function.
  Parameters(std::string &&passedFilename,
             std::vector<NumericType> &&passedAspectRatios,
             std::vector<NumericType> &&passedTaperAngles,
             std::vector<NumericType> &&passedStickingProbabilities,
             std::size_t passedTimeSteps, bool passedSymmetrical,
             bool passedFixTopWidth)
      : filename(std::move(passedFilename)),
        aspectRatios(std::move(passedAspectRatios)),
        taperAngles(std::move(passedTaperAngles)),
        stickingProbabilities(std::move(passedStickingProbabilities)),
        timeSteps(passedTimeSteps), symmetrical(passedSymmetrical),
        fixTopWidth(passedFixTopWidth) {}

public:
  // Factory pattern to construct an instance of the class from a map
  static Parameters
  fromMap(const std::unordered_map<std::string, std::string> &map) {
    // Local variable instances
    std::string filename = "data.csv";
    std::vector<NumericType> aspectRatios;
    std::vector<NumericType> taperAngles;
    std::vector<NumericType> stickingProbabilities;
    std::size_t timeSteps = 10;
    bool symmetrical = true;
    bool fixTopWidth = true;

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
        Utils::Item{"timeSteps", timeSteps, &Utils::toStrictlyPositive<int>},
        Utils::Item{
            "symmetrical",
            symmetrical,
            Utils::toBool,
        },
        Utils::Item{
            "fixTopWidth",
            fixTopWidth,
            Utils::toBool,
        });

    return Parameters(std::move(filename), std::move(aspectRatios),
                      std::move(taperAngles), std::move(stickingProbabilities),
                      timeSteps, symmetrical, fixTopWidth);
  }

  std::size_t getNumberOfCombinations() const {
    return aspectRatios.size() * taperAngles.size() *
           stickingProbabilities.size();
  }

  int getInputDimension() const { return 4; }

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
    std::cout << "- timeSteps: " << timeSteps << "\n";
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

#endif
