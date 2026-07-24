add_executable(common_tests test_main.cpp test_inputfile.cpp test_typses.cpp)

target_link_libraries(common_tests PRIVATE sirun_utils)

target_include_directories(common_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_test(NAME common_tests COMMAND common_tests)
