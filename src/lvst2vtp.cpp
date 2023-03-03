#include <filesystem>
#include <iostream>

#include <lsMesh.hpp>
#include <lsReader.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace fs = std::filesystem;

int main(int argc, const char *const *const argv) {
  using NumericType = double;
  static constexpr int D = 2;

  std::string lvstFilename;
  std::string vtpFilename;
  if (argc > 2) {
    lvstFilename = argv[1];
    vtpFilename = argv[2];

    if (!fs::is_regular_file(lvstFilename)) {
      std::cout << "Input file '" << lvstFilename << "' couldn't be found!"
                << std::endl;
      return 0;
    }
  } else {
    std::cout << "Usage:\n"
              << argv[0] << " <lvst_filename> <output_filename>" << std::endl;
    return 0;
  }

  auto levelset = lsSmartPointer<lsDomain<NumericType, D>>::New();

  lsReader<NumericType, D>(levelset, lvstFilename).apply();

  auto mesh = lsSmartPointer<lsMesh<>>::New();
  lsToSurfaceMesh<NumericType, D>(levelset, mesh).apply();
  lsVTKWriter<NumericType>(mesh, vtpFilename).apply();
}