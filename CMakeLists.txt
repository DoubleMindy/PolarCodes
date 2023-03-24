cmake_minimum_required(VERSION 3.10.0)

project(PolarCodes CXX)
set(CMAKE_CXX_STANDARD 17)
include(cmake/StandardProjectSettings.cmake)
# Link this 'library' to set the c++ standard / compile-time options requested
add_library(default INTERFACE)
target_compile_features(default INTERFACE cxx_std_17)

include(cmake/CompilerWarnings.cmake)
set_project_warnings(default)

option(DOMAIN_LLR "Use LLR domain, not P1 domain" ON)
target_compile_definitions(default INTERFACE $<IF:$<BOOL:${DOMAIN_LLR}>,DOMAIN_LLR=1,DOMAIN_P1=1>)

find_package(Matlab COMPONENTS MX_LIBRARY)
if(Matlab_FOUND)
    add_subdirectory(lib/mexbind0x)
endif()


add_subdirectory(src)
