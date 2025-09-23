cmake_minimum_required(VERSION 3.16)
project(second_mpi_mpfr LANGUAGES CXX)

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
# Components to link: mpi, serialization, thread, system, chrono
#   (ublas and multiprecision are header-only)
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
include_directories(${Boost_INCLUDE_DIRS})

# ---- MPFR & GMP (for Boost.Multiprecision mpfr backend) ----
find_package(PkgConfig REQUIRED)
pkg_check_modules(MPFR REQUIRED mpfr)
pkg_check_modules(GMP  REQUIRED gmp)

# Threads (Boost::thread needs it)
find_package(Threads REQUIRED)

# ---- Target ----
add_executable(second second.cpp)

# fmemopen (glibc) needs this
target_compile_definitions(second PRIVATE _GNU_SOURCE)

# Include dirs from pkg-config
target_include_directories(second PRIVATE
  ${MPFR_INCLUDE_DIRS}
  ${GMP_INCLUDE_DIRS}
)

# Link everything
target_link_libraries(second
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
  target_link_libraries(second PRIVATE mpfr)
endif()
if (NOT GMP_LINK_LIBRARIES)
  target_link_libraries(second PRIVATE gmp)
endif()

# Warnings (optional)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(second PRIVATE -Wall -Wextra -Wno-unknown-pragmas)
endif()

