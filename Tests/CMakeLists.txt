# Allow headers in tests to be included like
# #include "Mask.h" rather than needing
# #include "Mask/Mask.h"
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(TestMask TestMask.cpp)
target_link_libraries(TestMask ${Mask_libraries})
add_test(TestMask TestMask)

add_executable(TestMaskOperations TestMaskOperations.cpp)
target_link_libraries(TestMaskOperations ${Mask_libraries})
add_test(TestMaskOperations TestMaskOperations)

add_executable(TestMaskRead TestMaskRead.cpp)
target_link_libraries(TestMaskRead ${Mask_libraries})
add_test(TestMaskRead TestMaskRead ${CMAKE_SOURCE_DIR}/Tests/data/TestMask.mask)

add_executable(TestForegroundBackgroundSegmentMaskRead TestForegroundBackgroundSegmentMaskRead.cpp)
target_link_libraries(TestForegroundBackgroundSegmentMaskRead ${Mask_libraries})
add_test(TestForegroundBackgroundSegmentMaskRead TestForegroundBackgroundSegmentMaskRead
         ${CMAKE_SOURCE_DIR}/Tests/data/TestMask.fbmask)
