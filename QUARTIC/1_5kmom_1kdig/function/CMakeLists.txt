cmake_minimum_required(VERSION 3.16)
project(function_mpi_mpfr LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless user overrides
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# MPI
find_package(MPI REQUIRED)

# Boost libs we actually link (others are header-only)
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
find_package(Threads REQUIRED)

add_executable(function function.cpp)

# Needed for fmemopen on glibc
target_compile_definitions(function PRIVATE _GNU_SOURCE)

# Include dirs
target_include_directories(function PRIVATE ${Boost_INCLUDE_DIRS})

# Try to locate libs (optional; if not found we still link by name below)
find_library(MPC_LIB mpc)
find_library(MPFR_LIB mpfr)
find_library(GMP_LIB gmp)

# Core links (mirror your working command)
target_link_libraries(function
  PRIVATE
    MPI::MPI_CXX
    Boost::mpi
    Boost::serialization
    Boost::thread
    Boost::system
    Boost::chrono
    Threads::Threads
    # Prefer found paths, else fall back to bare names
    ${MPC_LIB} ${MPFR_LIB} ${GMP_LIB}
    mpc mpfr gmp
)

# Warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(function PRIVATE -Wall -Wextra -Wno-unknown-pragmas)
endif()

