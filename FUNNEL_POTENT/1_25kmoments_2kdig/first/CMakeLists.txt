cmake_minimum_required(VERSION 3.16)
project(first_mpi_mpfr LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless user overrides
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# ---- MPI ----
find_package(MPI REQUIRED)

# ---- Boost ----
# Components used by your source:
#   mpi, serialization (need linking)
#   thread, system, chrono (needed by Boost.Thread and timing)
# Ublas, lexical_cast, math are header-only => no link target required.
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
include_directories(${Boost_INCLUDE_DIRS})

# ---- MPFR & GMP (for Boost.Multiprecision mpfr backend) ----
find_package(PkgConfig REQUIRED)
pkg_check_modules(MPFR REQUIRED mpfr)
pkg_check_modules(GMP  REQUIRED gmp)

# Threads (pulled by Boost::thread but be explicit)
find_package(Threads REQUIRED)

# ---- Target ----
add_executable(first first.cpp)

# fmemopen (glibc) needs this
target_compile_definitions(first PRIVATE _GNU_SOURCE)

# Include dirs from pkg-config
target_include_directories(first PRIVATE
  ${MPFR_INCLUDE_DIRS}
  ${GMP_INCLUDE_DIRS}
)

# Link everything
target_link_libraries(first
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

# Fallback if pkg-config didn't report link libs
if (NOT MPFR_LINK_LIBRARIES)
  target_link_libraries(first PRIVATE mpfr)
endif()
if (NOT GMP_LINK_LIBRARIES)
  target_link_libraries(first PRIVATE gmp)
endif()

# Warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(first PRIVATE -Wall -Wextra -Wno-unknown-pragmas)
endif()

