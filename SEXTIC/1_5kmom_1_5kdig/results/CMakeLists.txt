cmake_minimum_required(VERSION 3.16)
project(result_mpi_mpfr LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless the user overrides
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# ---- MPI ----
find_package(MPI REQUIRED)

# ---- Boost ----
# Required linked components:
#   mpi, serialization, thread, system, chrono
#   (ublas, math, multiprecision are header-only)
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
include_directories(${Boost_INCLUDE_DIRS})

# ---- MPFR & GMP ----
# We can just link by name (like your manual command) without pkg-config
find_library(MPFR_LIB mpfr)
find_library(GMP_LIB gmp)

# Threads
find_package(Threads REQUIRED)

# ---- Target ----
add_executable(result result.cpp)

# Needed for fmemopen on glibc
target_compile_definitions(result PRIVATE _GNU_SOURCE)

# Link libraries
target_link_libraries(result
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
    mpfr   # fallback if not found with find_library
    gmp    # fallback
)

# Optional warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(result PRIVATE -Wall -Wextra -Wno-unknown-pragmas)
endif()

