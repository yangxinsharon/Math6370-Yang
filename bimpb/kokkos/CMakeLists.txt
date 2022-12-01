cmake_minimum_required (VERSION 3.10)
project (bimpbKokkos)

set(SPACK_CXX g++)
if(SPACK_CXX)
  message("found spack compiler g++")
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "the C++ compiler" FORCE)
  set(ENV{CXX} g++)
endif()

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

add_executable(bimpb_kokkos.exe main.cpp matvec.cpp readin.cpp)
target_link_libraries(bimpb_kokkos.exe Kokkos::kokkos)
