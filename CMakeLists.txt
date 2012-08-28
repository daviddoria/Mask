project(Mask)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

cmake_minimum_required(VERSION 2.8)

if(NOT TARGET ITKHelpers)
  add_subdirectory(ITKHelpers)
  include_directories(${ITKHelpers_includes})
  set(Mask_includes ${Mask_includes} ${ITKHelpers_includes})
  set(Mask_libraries ${Mask_libraries} ITKHelpers)
endif()

if(NOT ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(Mask_includes ${Mask_includes} ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
  set(Mask_libraries ${Mask_libraries} Mask PARENT_SCOPE)
  message("Mask submodule used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# This is needed for c++0x features used in Helpers.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology ITKDistanceMap)
  INCLUDE(${ITK_USE_FILE})
endif()

SET(Mask_BuildTests OFF CACHE BOOL "Build Mask tests?")
if(Mask_BuildTests)
  add_subdirectory(Tests)
endif(Mask_BuildTests)

add_library(Mask Mask.cpp MaskOperations.cpp)
target_link_libraries(Mask ITKHelpers)

add_custom_target(MaskSources SOURCES Mask.h MaskQt.h Mask.hpp MaskOperations.hpp)

# Only build the MaskQt library if the user wants it
SET(Mask_BuildMaskQt OFF CACHE BOOL "Build MaskQt?")
if(Mask_BuildMaskQt)
  FIND_PACKAGE(Qt4)
  INCLUDE(${QT_USE_FILE})
  add_library(MaskQt Mask.cpp MaskQt.cpp)
endif(Mask_BuildMaskQt)


GET_DIRECTORY_PROPERTY(AllMaskIncludes INCLUDE_DIRECTORIES)
message("AllMaskIncludes include directories: ${AllMaskIncludes}")
