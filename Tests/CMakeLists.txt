# Allow headers in tests to be included like
# #include "Mask.h" rather than needing
# #include "Mask/Mask.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestMask TestMask.cpp)
target_link_libraries(TestMask Mask)

add_executable(TestMaskOperations TestMaskOperations.cpp)
target_link_libraries(TestMaskOperations Mask)
