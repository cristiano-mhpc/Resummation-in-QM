cmake_minimum_required(VERSION 3.16)
project(lu_auto LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless overridden
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# MPI
find_package(MPI REQUIRED)

# Boost components to link (others are header-only)
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
find_package(Threads REQUIRED)

# Executable name must be 'auto' to match your run script
add_executable(auto LU.cpp)

# Needed for fmemopen on glibc
target_compile_definitions(auto PRIVATE _GNU_SOURCE)

# Include dirs
target_include_directories(auto PRIVATE ${Boost_INCLUDE_DIRS})

# Try to locate MPFR/GMP (optional). Fall back to bare names if not found.
find_library(MPFR_LIB mpfr)
find_library(GMP_LIB  gmp)

# Link exactly what your manual command links
target_link_libraries(auto
  PRIVATE
    MPI::MPI_CXX
    Boost::mpi
    Boost::serialization
    Boost::thread
    Boost::system
    Boost::chrono
    Threads::Threads
    ${MPFR_LIB}
    ${GMP_LIB}
    mpfr gmp              # fallback if find_library fails
)

# Warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(auto PRIVATE -O3 -Wall -Wextra -Wno-unknown-pragmas)
else()
  target_compile_options(auto PRIVATE -O3)
endif()

