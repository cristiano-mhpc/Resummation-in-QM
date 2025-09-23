cmake_minimum_required(VERSION 3.16)
project(constant LANGUAGES CXX)

# ---- C++ standard ----
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Target ----
add_executable(constant constant.cpp)

# ---- MPFR & GMP ----
# Try to locate the libraries (optional: CMake will fall back to bare names)
find_library(MPFR_LIB mpfr)
find_library(GMP_LIB  gmp)

# Add include directories if needed (usually /usr/include is enough)
# target_include_directories(constant PRIVATE /usr/include)

# Link libraries (prefer found paths, else fall back to names)
target_link_libraries(constant
    PRIVATE
        ${MPFR_LIB}
        ${GMP_LIB}
        mpfr
        gmp
)

# Optimization and warnings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(constant PRIVATE -O3 -Wall -Wextra)
else()
    target_compile_options(constant PRIVATE -O3)
endif()

