# SPDX-FileCopyrightText: Copyright (c) 2022 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

if (NOT OSD_LITE_FOLDER_NAME)
    set(OSD_LITE_FOLDER_NAME "${PROJECT_NAME}")
endif()

if (NOT CMAKE_DEBUG_POSTFIX)
    set(CMAKE_DEBUG_POSTFIX d)
endif()

add_subdirectory(sdc)
add_subdirectory(vtr)
add_subdirectory(far)
add_subdirectory(bfr)
add_subdirectory(tmr)

set(src_files version.cpp)

if (MSVC)
    # add natvis templates to introspect various array data structures in the debugger
    set(OPENSUBDIV_VERSION "v${PROJECT_VERSION_MAJOR}_${PROJECT_VERSION_MINOR}_${PROJECT_VERSION_PATCH}")
    configure_file(msvc_osd.natvis.in msvc_osd.natvis)
    if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/msvc_osd.natvis")
        list(APPEND src_files "${CMAKE_CURRENT_BINARY_DIR}/msvc_osd.natvis")
        message(STATUS "MSVC natvis templates added")
    else()
        message(WARNING "MSVC generator detected, but unable to configure natvis templates")
    endif()
endif()

configure_file("version.h.in" "${CMAKE_CURRENT_BINARY_DIR}/version.h")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/version.h" DESTINATION "include/opensubdiv")

include(GNUInstallDirs)

#
# Static
#

add_library(osd_lite_static STATIC ${src_files}
    $<TARGET_OBJECTS:osd_lite_sdc_obj>
    $<TARGET_OBJECTS:osd_lite_vtr_obj>
    $<TARGET_OBJECTS:osd_lite_far_obj>  
    $<TARGET_OBJECTS:osd_lite_bfr_obj>  
    $<TARGET_OBJECTS:osd_lite_tmr_obj>)

target_compile_options(osd_lite_static PUBLIC
    "$<TARGET_PROPERTY:osd_lite_tmr_obj,INTERFACE_COMPILE_OPTIONS>"
    "$<TARGET_PROPERTY:osd_lite_bfr_obj,INTERFACE_COMPILE_OPTIONS>"
    "$<TARGET_PROPERTY:osd_lite_far_obj,INTERFACE_COMPILE_OPTIONS>"
)

target_include_directories(osd_lite_static PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_include_directories(osd_lite_static PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..;${CMAKE_CURRENT_BINARY_DIR}/sdc;${CMAKE_CURRENT_BINARY_DIR}/..>" $<INSTALL_INTERFACE:include/>)

set_target_properties(osd_lite_static PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX} FOLDER ${OSD_LITE_FOLDER_NAME})

install(TARGETS osd_lite_static DESTINATION "${CMAKE_INSTALL_LIBDIR}" EXPORT "osd_lite_targets")

#
# Shared
#
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

add_library(osd_lite_shared SHARED ${src_files}
    $<TARGET_OBJECTS:osd_lite_sdc_obj>
    $<TARGET_OBJECTS:osd_lite_vtr_obj>
    $<TARGET_OBJECTS:osd_lite_far_obj>
    $<TARGET_OBJECTS:osd_lite_bfr_obj>
    $<TARGET_OBJECTS:osd_lite_tmr_obj>)

target_compile_options(osd_lite_shared PUBLIC
    "$<TARGET_PROPERTY:osd_lite_tmr_obj,INTERFACE_COMPILE_OPTIONS>"
    "$<TARGET_PROPERTY:osd_lite_bfr_obj,INTERFACE_COMPILE_OPTIONS>"
    "$<TARGET_PROPERTY:osd_lite_far_obj,INTERFACE_COMPILE_OPTIONS>"
)

target_include_directories(osd_lite_shared PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_include_directories(osd_lite_shared PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..;${CMAKE_CURRENT_BINARY_DIR}/sdc>" $<INSTALL_INTERFACE:include/>)

set_target_properties(osd_lite_shared PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX} FOLDER ${OSD_LITE_FOLDER_NAME})

install(TARGETS osd_lite_shared DESTINATION "${CMAKE_INSTALL_BINDIR}" EXPORT "osd_lite_targets")

#
# Export
#

include(CMakePackageConfigHelpers)

set(osd_lite_config_path "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${osd_lite_config_path}"
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  VERSION "${CMAKE_PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion)

install(EXPORT "osd_lite_targets" 
    FILE "${PROJECT_NAME}Targets.cmake" 
    NAMESPACE "osd::"
    DESTINATION "${osd_lite_config_path}")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "${osd_lite_config_path}")

