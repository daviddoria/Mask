cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# This is needed for shared_ptr and the trick using enable_if and if_fundamental to allow scalars
# to be treated as the 0th component of a vector.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x")

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

SET(Mask_BuildTests ON CACHE STRING "Build Mask tests?" FORCE)
if(Mask_BuildTests)
  add_subdirectory(Tests)
endif(Mask_BuildTests)

#add_library(Mask Mask.cpp MaskOperations.cpp MaskQt.cpp)
add_library(Mask Mask.cpp MaskOperations.cpp)
target_link_libraries(Mask ITKHelpers)
