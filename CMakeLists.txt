cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(Qt4 REQUIRED)
INCLUDE(${QT_USE_FILE})

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
  include_directories(${ITKHelpers_includes})
endif()

if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(Mask_includes ${ITKHelpers_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
endif()

# This is needed for c++0x features used in Helpers.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase
ITKImageIntensity ITKImageFeature ITKMathematicalMorphology ITKBinaryMathematicalMorphology ITKDistanceMap)
INCLUDE(${ITK_USE_FILE})

SET(Mask_BuildTests OFF CACHE BOOL "Build Mask tests?")
if(Mask_BuildTests)
  add_subdirectory(Tests)
endif(Mask_BuildTests)

add_library(Mask Mask.cpp MaskOperations.cpp)
target_link_libraries(Mask ITKHelpers)

# Only build the MaskQt file if Qt is available
FIND_PACKAGE(Qt4)
if(DEFINED QT4_FOUND)
  INCLUDE(${QT_USE_FILE})
  add_library(MaskQt Mask.cpp MaskQt.cpp)
endif()
