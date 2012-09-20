cmake_minimum_required(VERSION 2.8)

project(Mask)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

# Only include Helpers if it is not already available
get_property(FoundITKHelpers GLOBAL PROPERTY ITKHelpersIncludeDirs SET)
if(NOT FoundITKHelpers)
  add_subdirectory(ITKHelpers)
endif()

get_property(ITKHelpersIncludeDirs GLOBAL PROPERTY ITKHelpersIncludeDirs)
set(Mask_include_dirs ${Mask_include_dirs} ${ITKHelpersIncludeDirs})
get_property(ITKHelpersLibs GLOBAL PROPERTY ITKHelpersLibs)
set(Mask_libraries ${Mask_libraries} ${ITKHelpersLibs})

# This is needed for c++0x features used in Helpers.
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
endif(UNIX)

# VTK
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
set(Mask_libraries ${Mask_libraries} ${VTK_LIBRARIES})

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology ITKDistanceMap)
  INCLUDE(${ITK_USE_FILE})
endif()
set(Mask_libraries ${Mask_libraries} ${ITK_LIBRARIES})

# Give the compiler all of the required include directories
include_directories(${Mask_include_dirs})

# Create the library
add_library(Mask Mask.cpp MaskOperations.cpp)
target_link_libraries(Mask ${Mask_libraries})
set(Mask_libraries ${Mask_libraries} Mask)

# Add non-compiled sources to the project
add_custom_target(MaskSources SOURCES Mask.h MaskQt.h Mask.hpp MaskOperations.hpp)

# Only build the MaskQt library if the user wants it
SET(Mask_BuildMaskQt OFF CACHE BOOL "Build MaskQt?")
if(Mask_BuildMaskQt)
  FIND_PACKAGE(Qt4)
  INCLUDE(${QT_USE_FILE})
  add_library(MaskQt Mask.cpp MaskQt.cpp)
endif(Mask_BuildMaskQt)

# Allow this project to be detected and used as a submodule
set(Mask_include_dirs ${Mask_include_dirs} ${CMAKE_CURRENT_SOURCE_DIR})
set_property(GLOBAL PROPERTY MaskIncludeDirs ${Mask_include_dirs})
set_property(GLOBAL PROPERTY MaskLibs ${Mask_libraries})

# Build the tests if requested
SET(Mask_BuildTests OFF CACHE BOOL "Build Mask tests?")
if(Mask_BuildTests)
  add_subdirectory(Tests)
endif(Mask_BuildTests)
