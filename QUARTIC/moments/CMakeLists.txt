cmake_minimum_required(VERSION 3.16)
project(moments_mpi_mpfr LANGUAGES CXX)

# --- C++ standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prefer mpicxx unless user overrides
if (NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()

# --- MPI ---
find_package(MPI REQUIRED)

# --- Boost ---
# Required linked components: mpi, serialization, thread, system, chrono
# (ublas and multiprecision headers are header-only)
find_package(Boost 1.71 REQUIRED COMPONENTS mpi serialization thread system chrono)
find_package(Threads REQUIRED)

# --- MPFR / GMP / MPC ---
# Try to locate the libs; fall back to bare names if not found.
find_library(MPFR_LIB mpfr)
find_library(GMP_LIB  gmp)
find_library(MPC_LIB  mpc)

# --- Target ---
add_executable(moments moments.cpp)

# Needed for fmemopen on glibc systems
target_compile_definitions(moments PRIVATE _GNU_SOURCE)

# Include Boost headers
target_include_directories(moments PRIVATE ${Boost_INCLUDE_DIRS})

# Link libraries (using found full paths if available, otherwise names)
target_link_libraries(moments
  PRIVATE
    MPI::MPI_CXX
    Boost::mpi
    Boost::serialization
    Boost::thread
    Boost::system
    Boost::chrono
    Threads::Threads
    ${MPC_LIB} ${MPFR_LIB} ${GMP_LIB}
    mpc mpfr gmp       # fallback if find_library fails
)

# Optional: add optimization and warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(moments PRIVATE -O3 -Wall -Wextra -Wno-unknown-pragmas)
else()
  target_compile_options(moments PRIVATE -O3)
endif()

