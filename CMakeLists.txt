cmake_minimum_required (VERSION 3.10)

project(GSP)

# Auto detect SDK path on macOS
if(APPLE)
  if(NOT DEFINED CMAKE_OSX_SYSROOT)
    execute_process(COMMAND xcrun --sdk macosx --show-sdk-path
                    OUTPUT_VARIABLE CMAKE_OSX_SYSROOT
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
endif()

add_executable(GSP src/utilities.h src/ann_data.h src/input_param.h src/var_list.h src/pop_data.h src/main_analysis_module.h src/thread_analysis_module.h src/segment_analysis_module.h src/segment_output_module.h src/segment_integrated_analysis.h src/GSP.cpp)

target_compile_features(GSP PUBLIC cxx_std_17)

set_target_properties(GSP PROPERTIES CXX_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED True)

# Remove the stdlib flag to use default
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)

target_link_libraries(GSP PRIVATE Threads::Threads)

# Add system include directories for macOS
if(APPLE)
  target_include_directories(GSP SYSTEM PRIVATE "${CMAKE_OSX_SYSROOT}/usr/include/c++/v1")
endif()

project(test_program)

add_executable(Postprocessing src/Postprocessing.cpp)
if(APPLE)
  target_include_directories(Postprocessing SYSTEM PRIVATE "${CMAKE_OSX_SYSROOT}/usr/include/c++/v1")
endif()

project(utility_program)

add_executable(Get_trimming_parameter_len utilities/Get_trimming_parameter_len.cpp)
if(APPLE)
  target_include_directories(Get_trimming_parameter_len SYSTEM PRIVATE "${CMAKE_OSX_SYSROOT}/usr/include/c++/v1")
endif()

add_executable(Get_trimming_parameter_qual utilities/Get_trimming_parameter_qual.cpp)
if(APPLE)
  target_include_directories(Get_trimming_parameter_qual SYSTEM PRIVATE "${CMAKE_OSX_SYSROOT}/usr/include/c++/v1")
endif()

add_executable(Get_signature_clusters utilities/Get_signature_clusters.cpp)
if(APPLE)
  target_include_directories(Get_signature_clusters SYSTEM PRIVATE "${CMAKE_OSX_SYSROOT}/usr/include/c++/v1")
endif()

# Add test executable for segment module
add_executable(test_segment_module test/test_segment_module.cpp)
target_compile_features(test_segment_module PUBLIC cxx_std_17)
if(APPLE)
  target_include_directories(test_segment_module SYSTEM PRIVATE "${CMAKE_OSX_SYSROOT}/usr/include/c++/v1")
endif()