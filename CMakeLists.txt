cmake_minimum_required(VERSION 3.12.1)
project(Regenerator)

macro(configure_files srcDir destDir)
    message(STATUS "Configuring directory ${destDir}")
    make_directory(${destDir})

    file(GLOB templateFiles RELATIVE ${srcDir} ${srcDir}/*)
    foreach(templateFile ${templateFiles})
        set(srcTemplatePath ${srcDir}/${templateFile})
        if(NOT IS_DIRECTORY ${srcTemplatePath})
            message(STATUS "Configuring file ${templateFile}")
            configure_file(
                    ${srcTemplatePath}
                    ${destDir}/${templateFile}
                    @ONLY)
        endif(NOT IS_DIRECTORY ${srcTemplatePath})
    endforeach(templateFile)
endmacro(configure_files)

set(BUILD_PATH \"${PROJECT_BINARY_DIR}\")
configure_file("${PROJECT_SOURCE_DIR}/build_config.h.in" "${PROJECT_SOURCE_DIR}/src/build_config.h")
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/spdlog)
set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

file(GLOB_RECURSE REGEN_SOURCES "src/*.cpp" "src/*.h" "spdlog/spdlog.h")
configure_files(${PROJECT_SOURCE_DIR}/PropertyFiles ${PROJECT_BINARY_DIR}/PropertyFiles)


	
add_executable(Regenerator ${REGEN_SOURCES})

