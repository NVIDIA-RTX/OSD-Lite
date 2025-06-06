#
#   Copyright 2013 Pixar
#
#   Licensed under the Apache License, Version 2.0 (the "Apache License")
#   with the following modification; you may not use this file except in
#   compliance with the Apache License and the following modification to it:
#   Section 6. Trademarks. is deleted and replaced with:
#
#   6. Trademarks. This License does not grant permission to use the trade
#      names, trademarks, service marks, or product names of the Licensor
#      and its affiliates, except as required to comply with Section 4(c) of
#      the License and to reproduce the content of the NOTICE file.
#
#   You may obtain a copy of the Apache License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the Apache License with the above modification is
#   distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#   KIND, either express or implied. See the Apache License for the specific
#   language governing permissions and limitations under the Apache License.
#

set(PUBLIC_HEADER_FILES
    bilinearScheme.h
    catmarkScheme.h
    crease.h
    loopScheme.h
    options.h
    scheme.h
    types.h
)

set(PRIVATE_HEADER_FILES )

# source & headers
set(SOURCE_FILES
    crease.cpp
    typeTraits.cpp
)

add_library(osd_lite_sdc_obj
    OBJECT
        ${SOURCE_FILES}
        ${PRIVATE_HEADER_FILES}
        ${PUBLIC_HEADER_FILES}
)

target_include_directories(osd_lite_sdc_obj PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

set_target_properties(osd_lite_sdc_obj PROPERTIES POSITION_INDEPENDENT_CODE TRUE FOLDER "${OSD_LITE_FOLDER_NAME}")

install(
    FILES
        ${PUBLIC_HEADER_FILES}
    DESTINATION
        "include/opensubdiv/sdc"
    PERMISSIONS
        OWNER_READ
        GROUP_READ
        WORLD_READ )

