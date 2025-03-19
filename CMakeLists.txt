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

cmake_minimum_required(VERSION 3.22)

project(OSD_lite VERSION 3.6.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(CMAKE_INSTALL_MESSAGE LAZY)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
    set(CMAKE_COMPILER_IS_CLANGCC 1)
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
    set(CMAKE_COMPILER_IS_GNUCC 1)
elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
    set(CMAKE_COMPILER_IS_ICC 1)
endif()

if (MSVC)

    # Turn on all warnings
    add_definitions(/W3 /WX)

    # Make sure WinDef.h does not define min and max macros which
    # will conflict with std::min() and std::max().
    add_definitions(/DNOMINMAX /D_USE_MATH_DEFINES)
    
    # disable _snprintf warning
    add_definitions(/wd4996)

    # Eliminate deprecation warnings for the older, less secure
    # CRT functions
    add_definitions(/D_CRT_SECURE_NO_WARNINGS)

else()

    # Turn on all warnings
    if(CMAKE_COMPILER_IS_ICC)
        add_definitions(-w2 -wd1572 -wd1418 -wd981 -wd383 -wd193 -wd444)
    elseif(CMAKE_COMPILER_IS_CLANGCC OR CMAKE_COMPILER_IS_GNUCC)
        add_definitions(-Wall -Wextra)
    endif()

    add_definitions(-Werror)

endif()

add_subdirectory(opensubdiv)

install(FILES "LICENSE.txt" "NOTICE.txt" DESTINATION ".")

export(EXPORT "osd_lite_targets" FILE "${CMAKE_CURRENT_BINARY_DIR}/osd_liteTargets.cmake")

#
# Package
#

if (WIN32)
    set(os_suffix "win64")
elseif (UNIX AND NOT APPLE)
    set(os_suffix "linux64")
endif()

set(CPACK_GENERATOR "ZIP")

set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}.${CMAKE_PROJECT_VERSION}.${os_suffix}")

include(CPack)
