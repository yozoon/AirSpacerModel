#ifndef UTILS_HPP
#define UTILS_HPP

#include <array>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <lsDomain.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>
#include <lsSmartPointer.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace Utils {
template <typename NumericType>
inline constexpr NumericType PI = std::acos(NumericType{-1});

template <typename NumericType>
inline constexpr NumericType deg2rad(NumericType deg) {
  return deg * PI<NumericType> / NumericType{180.};
}

template <typename NumericType>
inline constexpr NumericType rad2deg(NumericType rad) {
  return rad * NumericType{180.} / PI<NumericType>;
}

template <typename NumericType, int D>
void printSurface(lsSmartPointer<lsDomain<NumericType, D>> levelset,
                  const std::string &filename) {
  auto mesh = lsSmartPointer<lsMesh<>>::New();
  lsToSurfaceMesh<NumericType, D>(levelset, mesh).apply();
  lsVTKWriter<NumericType>(mesh, filename).apply();
}

template <typename NumericType, int D>
std::array<NumericType, 2 * D> getBoundsFromGrid(const hrleGrid<D> &grid) {
  std::array<NumericType, 2 * D> bounds{0.};
  for (unsigned i = 0; i < D; ++i) {
    // Retrieve min and max bounds for each direction that does not have
    // infinite boundary condition
    if (grid.getBoundaryConditions(i) !=
        lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY) {
      auto minIndex = grid.getMinGridPoint(i);
      auto maxIndex = grid.getMaxGridPoint(i);
      bounds[2 * i] = grid.index2Coordinate(minIndex);
      bounds[2 * i + 1] = grid.index2Coordinate(maxIndex);
    }
  }
  return bounds;
}

// Checks if a string starts with a - or not
bool isSigned(const std::string &s) {
  auto pos = s.find_first_not_of(' ');
  if (pos == s.npos)
    return false;
  if (s[pos] == '-')
    return true;
  return false;
}

// Converts string to the given numeric datatype
template <typename T> T convert(const std::string &s) {
  if constexpr (std::is_same_v<T, int>) {
    return std::stoi(s);
  } else if constexpr (std::is_same_v<T, unsigned int>) {
    if (isSigned(s))
      throw std::invalid_argument("The value must be unsigned");
    unsigned long int val = std::stoul(s);
    unsigned int num = static_cast<unsigned int>(val);
    if (val != num)
      throw std::out_of_range("The value is larger than the largest value "
                              "representable by `unsigned int`.");
    return num;
  } else if constexpr (std::is_same_v<T, long int>) {
    return std::stol(s);
  } else if constexpr (std::is_same_v<T, unsigned long int>) {
    if (isSigned(s))
      throw std::invalid_argument("The value must be unsigned");
    return std::stoul(s);
  } else if constexpr (std::is_same_v<T, long long int>) {
    return std::stoll(s);
  } else if constexpr (std::is_same_v<T, unsigned long long int>) {
    if (isSigned(s))
      throw std::invalid_argument("The value must be unsigned");
    return std::stoull(s);
  } else if constexpr (std::is_same_v<T, float>) {
    return std::stof(s);
  } else if constexpr (std::is_same_v<T, double>) {
    return std::stod(s);
  } else if constexpr (std::is_same_v<T, long double>) {
    return std::stold(s);
  } else if constexpr (std::is_same_v<T, std::string>) {
    return s;
  } else {
    // Throws a compile time error for all types but void
    return;
  }
}

template <typename V, typename C>
std::vector<V> toVector(const std::string &s, C conv) {
  const char delimiter = ',';
  std::istringstream istr(s);
  std::vector<V> vec;
  std::string item;
  while (std::getline(istr, item, delimiter)) {
    try {
      vec.push_back(conv(item));
    } catch (std::exception &e) {
      lsMessage::getInstance()
          .addError(
              std::string(
                  "Value couldn't be converted to requested type. Reason: ") +
              std::string(e.what()))
          .print();
    }
  }
  return vec;
}

template <typename V> inline std::vector<V> toVector(const std::string &s) {
  return toVector<V, decltype(&Utils::convert<V>)>(s, &Utils::convert<V>);
}

// Special conversion functions
template <typename NumericType>
inline NumericType toStrictlyPositive(const std::string &s) {
  auto value = Utils::convert<NumericType>(s);
  if (value <= 0.0)
    throw std::invalid_argument("Value must be strictly positive.");
  return value;
};

template <typename NumericType>
inline NumericType toUnitRange(const std::string &s) {
  auto value = Utils::convert<NumericType>(s);
  if (value > 1.0 || value <= 0.0)
    throw std::invalid_argument("Value must be in the range [1,0).");
  return value;
};

bool toBool(const std::string &s) {
  auto lower = s;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  bool value = false;
  if (lower == "true" || lower == "1") {
    value = true;
  } else if (lower == "false" || lower == "0") {
    value = false;
  } else {
    throw std::invalid_argument("Failed to convert value to boolean.");
  }
  return value;
};

template <class Iterator>
std::string join(Iterator begin, Iterator end,
                 const std::string &separator = ",") {
  std::ostringstream ostr;
  if (begin != end)
    ostr << *begin++;
  while (begin != end)
    ostr << separator << *begin++;
  return ostr.str();
}

// safeConvert wraps the convert function to catch exceptions. If an error
// occurs the default initialized value is returned.
template <typename T> std::optional<T> safeConvert(const std::string &s) {
  T value;
  try {
    value = convert<T>(s);
  } catch (std::exception &e) {
    lsMessage::getInstance()
        .addWarning(std::string("`") + s +
                    std::string("` couldn't be converted to type  `") +
                    typeid(value).name() + std::string("`"))
        .print();
    return std::nullopt;
  }
  return {value};
}

std::string trim(const std::string &str,
                 const std::string &whitespace = " \t") {
  const auto strBegin = str.find_first_not_of(whitespace);
  if (strBegin == std::string::npos)
    return std::string("");

  const auto strEnd = str.find_last_not_of(whitespace);
  const auto strRange = strEnd - strBegin + 1;

  return str.substr(strBegin, strRange);
}

std::unordered_map<std::string, std::string>
parseConfigStream(std::istream &input) {
  // Reads a simple config file containing a single <key>=<value> pair per line
  // and returns the content as an unordered map
  std::unordered_map<std::string, std::string> paramMap;
  std::string line;
  while (std::getline(input, line)) {
    // Remove trailing and leading whitespaces
    line = trim(line);

    // Skip this line if it is marked as a comment
    if (line.find('#') == 0 || line.empty())
      continue;

    auto splitPos = line.find("=");
    if (splitPos == std::string::npos)
      continue;

    auto keyStr = trim(line.substr(0, splitPos));
    auto valStr = trim(line.substr(splitPos + 1, line.length()));

    if (keyStr.empty() || valStr.empty())
      continue;

    paramMap.insert({keyStr, valStr});
  }
  return paramMap;
}

// Opens a file and forwards its stream to the config stream parser.
std::unordered_map<std::string, std::string>
readConfigFile(const std::string &filename) {
  std::ifstream f(filename);
  if (!f.is_open()) {
    lsMessage::getInstance()
        .addWarning(std::string("Failed to open config file `") + filename +
                    std::string("`"))
        .print();
    return {};
  }
  return parseConfigStream(f);
}

// Class that can be used during the assigning process of a param map to the
// param struct
template <typename K, typename V, typename C = decltype(&convert<V>)>
class Item {
private:
  C conv;

public:
  K key;
  V &value;

  Item(K key_, V &value_) : conv(&convert<V>), key(key_), value(value_) {}

  Item(K key_, V &value_, C conv_) : conv(conv_), key(key_), value(value_) {}

  void operator()(const std::string &k) {
    try {
      value = conv(k);
    } catch (std::exception &e) {
      lsMessage::getInstance()
          .addError(
              std::string("`") + k +
              std::string("` couldn't be converted to type of parameter `") +
              key + std::string("`. Reason: ") + std::string(e.what()))
          .print();
    }
  }
};

// If the key is found inthe unordered_map, then the
template <typename K, typename V, typename C>
void AssignItems(const std::unordered_map<std::string, std::string> &map,
                 Item<K, V, C> &&item) {
  if (auto it = map.find(item.key); it != map.end()) {
    item(it->second);
  } else {
    lsMessage::getInstance()
        .addWarning(
            std::string("Couldn't find `") + item.key +
            std::string("` in parameter file. Using default value instead."))
        .print();
  }
}

// Peels off items from parameter pack
template <typename K, typename V, typename C, typename... ARGS>
void AssignItems(const std::unordered_map<std::string, std::string> &map,
                 Item<K, V, C> &&item, ARGS &&...args) {
  AssignItems(map, std::forward<Item<K, V, C>>(item));
  AssignItems(map, std::forward<ARGS>(args)...);
}

}; // namespace Utils
#endif