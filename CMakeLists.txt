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

# Enable c++11
if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif(UNIX)

# VTK
option(MaskUseVTK "Enable VTK functions in MaskOperations?" OFF)
if(MaskUseVTK)
  FIND_PACKAGE(VTK)
  if(VTK_FOUND)
    INCLUDE(${VTK_USE_FILE})
    set(Mask_libraries ${Mask_libraries} ${VTK_LIBRARIES})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMASK_USE_VTK")
  endif()
endif()

# ITK
if(NOT ITK_FOUND)
  FIND_PACKAGE(ITK REQUIRED ITKCommon ITKIOImageBase ITKIOPNG ITKIOMeta
                            ITKImageIntensity ITKImageFeature ITKMathematicalMorphology
                            ITKBinaryMathematicalMorphology ITKDistanceMap)
  INCLUDE(${ITK_USE_FILE})

  if( "${ITK_VERSION_MAJOR}" LESS 4 )
    MESSAGE(FATAL_ERROR "You must build ImageGraphCutSegmentation with ITK >= 4.0!")
  endif( "${ITK_VERSION_MAJOR}" LESS 4 )
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
option(Mask_BuildTests "Build Mask tests?" OFF)
if(Mask_BuildTests)
  add_subdirectory(Tests)
endif(Mask_BuildTests)

# Display where this code was used from (if it is used as a submodule, there may be multiple instances of this submodule in the project, only the first of which is used)
option(Mask_ShowSubmoduleLocation "Show the path from which Mask was used?" OFF)
if(Mask_ShowSubmoduleLocation)
  message("Mask used from ${CMAKE_CURRENT_SOURCE_DIR}")
endif(Mask_ShowSubmoduleLocation)
