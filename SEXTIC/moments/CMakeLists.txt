cmake_minimum_required(VERSION 3.16)
project(moments_mpi_mpfr LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless user overrides
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# MPI
find_package(MPI REQUIRED)

# Boost: headers + these components are used/linked per your script
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
include_directories(${Boost_INCLUDE_DIRS})

# MPFR & GMP (via pkg-config)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MPFR REQUIRED mpfr)
pkg_check_modules(GMP  REQUIRED gmp)

# Threads (Boost::thread will need it)
find_package(Threads REQUIRED)

add_executable(moments moments.cpp)

# fmemopen requires this on glibc systems
target_compile_definitions(moments PRIVATE _GNU_SOURCE)

# Include dirs from pkg-config
target_include_directories(moments PRIVATE
  ${MPFR_INCLUDE_DIRS}
  ${GMP_INCLUDE_DIRS}
)

# Link everything
target_link_libraries(moments
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

# Fallback if pkg-config didn't set link libs (some distros)
if (NOT MPFR_LINK_LIBRARIES)
  target_link_libraries(moments PRIVATE mpfr)
endif()
if (NOT GMP_LINK_LIBRARIES)
  target_link_libraries(moments PRIVATE gmp)
endif()

# Warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(moments PRIVATE -Wall -Wextra -Wno-unknown-pragmas)
endif()

