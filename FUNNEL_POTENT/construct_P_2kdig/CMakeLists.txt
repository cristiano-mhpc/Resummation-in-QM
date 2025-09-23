cmake_minimum_required(VERSION 3.16)
project(lu_mpi_mpfr LANGUAGES CXX)

# Use C++17 (adjust if you prefer 14/20)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless user overrides
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# MPI
find_package(MPI REQUIRED)

# Boost (headers + libs used by the code)
# Components: mpi, serialization, thread, system, chrono
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
include_directories(${Boost_INCLUDE_DIRS})

# MPFR & GMP via pkg-config (common on Ubuntu)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MPFR REQUIRED mpfr)
pkg_check_modules(GMP  REQUIRED gmp)

# Threads (Boost::thread will pull it, but be explicit)
find_package(Threads REQUIRED)

add_executable(lu_mpi LU.cpp)

# fmemopen needs _GNU_SOURCE on glibc
target_compile_definitions(lu_mpi PRIVATE _GNU_SOURCE)

# Include dirs from pkg-config
target_include_directories(lu_mpi PRIVATE
  ${MPFR_INCLUDE_DIRS}
  ${GMP_INCLUDE_DIRS}
)

# Link everything
target_link_libraries(lu_mpi
  PRIVATE
    MPI::MPI_CXX
    Boost::mpi
    Boost::serialization
    Boost::thread
    Boost::system
    Boost::chrono
    Threads::Threads
    ${MPFR_LINK_LIBRARIES}
    ${GMP_LINK_LIBRARIES}
)

# If your pkg-config variables are empty on your system,
# you can fall back explicitly to mpfr/gmp:
if (NOT MPFR_LINK_LIBRARIES)
  target_link_libraries(lu_mpi PRIVATE mpfr)
endif()
if (NOT GMP_LINK_LIBRARIES)
  target_link_libraries(lu_mpi PRIVATE gmp)
endif()

# Optional: treat warnings nicely
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(lu_mpi PRIVATE -Wall -Wextra -Wno-unknown-pragmas)
endif()

