cmake_minimum_required(VERSION 2.8)

project(Mask)
SET(CMAKE_INCLUDE_CURRENT_DIR ON)
enable_testing()

add_subdirectory(CMakeHelpers)

UseSubmodule(ITKHelpers Mask)

# Enable c++11
if(UNIX)
  message("System is UNIX - enabling C++11")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAS_CPP11") # Can't use HAS_C++11 - it doesn't like the + characters in the define
endif(UNIX)

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology ITKDistanceMap)
  INCLUDE(${ITK_USE_FILE})

  if( "${ITK_VERSION_MAJOR}" LESS 4 )
    MESSAGE(FATAL_ERROR "You must build Mask with ITK >= 4.0!")
  endif( "${ITK_VERSION_MAJOR}" LESS 4 )
endif()
set(Mask_libraries ${Mask_libraries} ${ITK_LIBRARIES})

# Give the compiler all of the required include directories
include_directories(${Mask_include_dirs})

# Create the library
add_library(Mask Mask.cpp MaskOperations.cpp
ForegroundBackgroundSegmentMask.cpp)
target_link_libraries(Mask ${Mask_libraries})
set(Mask_libraries ${Mask_libraries} Mask)

# Add non-compiled sources to the project
add_custom_target(MaskSources SOURCES
ForegroundBackgroundSegmentMask.h
Mask.h
Mask.hpp
MaskQt.h
MaskVTK.h
MaskVTK.hpp
MaskOperations.hpp
SegmentMask.h)

# Only build the MaskQt library if the user wants it
option(Mask_BuildMaskQt "Build MaskQt?")
if(Mask_BuildMaskQt)
  FIND_PACKAGE(Qt4 REQUIRED)
  INCLUDE(${QT_USE_FILE})

  add_library(MaskQt MaskQt.cpp)
  target_link_libraries(MaskQt Mask)
  set(Mask_libraries ${Mask_libraries} MaskQt)
endif(Mask_BuildMaskQt)

# Only build the MaskVTK library if the user wants it
option(Mask_BuildMaskVTK "Build MaskQt?")
if(Mask_BuildMaskVTK)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
  set(Mask_libraries ${Mask_libraries} ${VTK_LIBRARIES})

  add_library(MaskVTK MaskVTK.cpp)
  target_link_libraries(MaskVTK Mask)
  set(Mask_libraries ${Mask_libraries} MaskVTK)
endif(Mask_BuildMaskVTK)

# Allow this project to be detected and used as a submodule
CreateSubmodule(Mask)

# Build the tests if requested
option(Mask_BuildTests "Build Mask tests?" OFF)
if(Mask_BuildTests)
  add_subdirectory(Tests)
endif(Mask_BuildTests)
