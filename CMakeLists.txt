cmake_minimum_required (VERSION 3.10)

project(GSEXP)

add_executable(GSEXP src/utilities.h src/ann_data.h src/input_param.h src/var_list.h src/pop_data.h src/main_analysis_module.h src/thread_analysis_module.h src/PFA-0.99.cpp)

target_compile_features(GSEXP PUBLIC cxx_std_17)

set_target_properties(GSEXP PROPERTIES CXX_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED True)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)

target_link_libraries(GSEXP PRIVATE Threads::Threads)

project(test_program)

add_executable(Postprocessing src/Postprocessing.cpp)

project(utility_program)

add_executable(Get_trimming_parameter_len utilities/Get_trimming_parameter_len.cpp)

add_executable(Get_trimming_parameter_qual utilities/Get_trimming_parameter_qual.cpp)