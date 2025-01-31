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

# source & headers
set(SOURCE_FILES
    bilinearPatchBuilder.cpp
    catmarkPatchBuilder.cpp
    error.cpp
    loopPatchBuilder.cpp
    patchBasis.cpp
    patchBuilder.cpp
    patchDescriptor.cpp
    patchMap.cpp
    patchTable.cpp
    patchTableFactory.cpp
    ptexIndices.cpp
    stencilTable.cpp
    stencilTableFactory.cpp
    stencilBuilder.cpp
    topologyDescriptor.cpp
    topologyRefiner.cpp
    topologyRefinerFactory.cpp
)

set(PRIVATE_HEADER_FILES
    bilinearPatchBuilder.h
    catmarkPatchBuilder.h
    loopPatchBuilder.h
    patchBuilder.h
    stencilBuilder.h
)

set(PUBLIC_HEADER_FILES
    error.h
    patchBasis.h    
    patchDescriptor.h
    patchParam.h
    patchMap.h
    patchTable.h
    patchTableFactory.h
    primvarRefiner.h
    ptexIndices.h
    sparseMatrix.h
    stencilTable.h
    stencilTableFactory.h
    topologyDescriptor.h
    topologyLevel.h
    topologyRefiner.h
    topologyRefinerFactory.h
    types.h
)


# Compile objs first for both the CPU and GPU libs -----
add_library(osd_lite_far_obj
    OBJECT
        ${SOURCE_FILES}
        ${PRIVATE_HEADER_FILES}
        ${PUBLIC_HEADER_FILES}
)

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_CLANGCC)
    target_compile_options(osd_lite_far_obj PUBLIC -Wno-unused-but-set-variable)
endif()

target_include_directories(osd_lite_far_obj PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

set_target_properties(osd_lite_far_obj PROPERTIES POSITION_INDEPENDENT_CODE TRUE FOLDER "${OSD_LITE_FOLDER_NAME}")

install(TARGETS osd_lite_far_obj EXPORT "osd_lite_targets")

install(
    FILES
        ${PUBLIC_HEADER_FILES}
    DESTINATION
        "include/opensubdiv/far"
    PERMISSIONS
        OWNER_READ
        GROUP_READ
        WORLD_READ )

#-------------------------------------------------------------------------------
