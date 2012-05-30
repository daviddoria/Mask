add_executable(TestMask TestMask.cpp)
target_link_libraries(TestMask Mask)

add_executable(TestMaskOperations TestMaskOperations.cpp)
target_link_libraries(TestMaskOperations Mask)
