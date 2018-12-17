cmake_minimum_required(VERSION 2.8.12.2)
project(Regenerator CXX)

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

include_directories(.)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 -pthread -O3 -flto")

file(GLOB_RECURSE REGEN_SOURCES "src/*.cpp" "src/*.h")
configure_files(${PROJECT_SOURCE_DIR}/PropertyFiles ${PROJECT_BINARY_DIR}/PropertyFiles)

add_executable(Regenerator ${REGEN_SOURCES})
