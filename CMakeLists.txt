cmake_minimum_required(VERSION 2.8.4)

if (POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)
endif ()

if (POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif ()

PROJECT(GATest)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/mex-it/cmake_modules")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DMX_COMPAT_32 -DMATLAB_MEX_FILE -lm -lpython2.7")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

message(STATUS "Looking for MATLAB")

find_package(MATLAB REQUIRED COMPONENTS MX_LIBRARY ENG_LIBRARY)
if (MATLAB_FOUND)
    message(STATUS "Found Matlab at dir:     ${Matlab_ROOT_DIR}")
    message(STATUS "Matlab include dir:      ${Matlab_INCLUDE_DIRS}")
else()
    message( STATUS "No Matlab Found!" )
endif ()

find_package(PythonLibs REQUIRED)

include_directories(${Matlab_INCLUDE_DIRS} "mex-it")

set(LINK_LIBRARIES ${Matlab_LIBRARIES} ${PYTHON_LIBRARIES})
aux_source_directory("./ga" GA_SRC)

# message( STATUS "GASRC:" ${GA_SRC})

set(MEXES locator_optimize)

# This comes from FindMATLAB.cmake
FOREACH (t ${MEXES})
    add_matlab_mex(NAME ${t} SRC ${t}.cpp ${GA_SRC} LINK_TO ${LINK_LIBRARIES})
    add_custom_command(
            TARGET ${t}
            PRE_BUILD
            COMMAND sed '9d' ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp > ${CMAKE_CURRENT_SOURCE_DIR}/mex_function.cpp
    )
ENDFOREACH (t ${MEXES})

# Standalone example with external mex function - Linux/Mac Only
if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_executable(main main.cpp mex_function.cpp ${GA_SRC})

    target_link_libraries(main ${LINK_LIBRARIES})
    set_property(TARGET main PROPERTY POSITION_INDEPENDENT_CODE TRUE)
    set_property(TARGET main PROPERTY CXX_STANDARD 11)
endif ()
