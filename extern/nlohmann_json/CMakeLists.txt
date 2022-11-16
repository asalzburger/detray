cmake_minimum_required( VERSION 3.11 )
include( FetchContent )

# Silence FetchContent warnings with CMake >=3.24.
if( POLICY CMP0135 )
   cmake_policy( SET CMP0135 NEW )
endif()

# Tell the user what's happening.
message( STATUS "Building nlohmann_json as part of the detray project" )

# Declare where to get VecMem from.
set( DETRAY_NLOHMANN_JSON_GIT_TAG "v3.10.5" CACHE STRING "Version of nlohmann_json to build" )
set( DETRAY_NLOHMANN_JSON_SHA1 "8969f5ad1a422e01f040ff48dcae9c0e6ad0811d" CACHE STRING "SHA1 hash of the downloaded zip" )
mark_as_advanced( DETRAY_NLOHMANN_JSON_GIT_REPOSITORY DETRAT_NLOHMANN_JSON_GIT_TAG )
FetchContent_Declare( nlohmann_json
   URL "https://github.com/nlohmann/json/archive/refs/tags/${DETRAY_NLOHMANN_JSON_GIT_TAG}.tar.gz"
   URL_HASH SHA1=${DETRAY_NLOHMANN_JSON_SHA1})

# Now set up its build.
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install ON CACHE INTERNAL "")
FetchContent_MakeAvailable( nlohmann_json )