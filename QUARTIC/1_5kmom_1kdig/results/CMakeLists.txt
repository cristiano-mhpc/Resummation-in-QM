cmake_minimum_required(VERSION 3.16)
project(result LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Prefer mpicxx if MPI is required, else g++ works too ----
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# ---- MPI ----
find_package(MPI REQUIRED)

# ---- Boost ----
# Required linked components:
#   mpi, serialization, thread, system, chrono
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
find_package(Threads REQUIRED)

# ---- MPFR & GMP ----
find_library(MPFR_LIB mpfr)
find_library(GMP_LIB  gmp)

# ---- Target ----
add_executable(result result.cpp)

# Needed for fmemopen (GNU extension)
target_compile_definitions(result PRIVATE _GNU_SOURCE)

# Include directories
target_include_directories(result PRIVATE ${Boost_INCLUDE_DIRS})

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
    mpfr
    gmp
)

# Optional: add optimization and warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(result PRIVATE -O3 -Wall -Wextra -Wno-unknown-pragmas)
else()
  target_compile_options(result PRIVATE -O3)
endif()

