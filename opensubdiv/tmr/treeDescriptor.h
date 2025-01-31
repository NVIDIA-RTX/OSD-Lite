//
//   Copyright 2016 Nvidia
//
//   Licensed under the Apache License, Version 2.0 (the "Apache License")
//   with the following modification; you may not use this file except in
//   compliance with the Apache License and the following modification to it:
//   Section 6. Trademarks. is deleted and replaced with:
//
//   6. Trademarks. This License does not grant permission to use the trade
//      names, trademarks, service marks, or product names of the Licensor
//      and its affiliates, except as required to comply with Section 4(c) of
//      the License and to reproduce the content of the NOTICE file.
//
//   You may obtain a copy of the Apache License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the Apache License with the above modification is
//   distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//   KIND, either express or implied. See the Apache License for the specific
//   language governing permissions and limitations under the Apache License.
//


#ifndef OPENSUBDIV3_TMR_TREE_DESCRIPTOR_H
#define OPENSUBDIV3_TMR_TREE_DESCRIPTOR_H

#ifdef __cplusplus
#include "../tmr/types.h"
#include <cstdint>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {
namespace Tmr {
#endif

/// \brief Set of bitfields for SubdivisionPlan quadtrees (aka 'root node')
/// 
///  Field              | Bits | Content
///  -------------------|:----:|---------------------------------------------------
///  FIELD0             |      |
///  irreg face         | 1    | true if the base face has an irregular size
///  (ununsed)          | 15   |
///  face size          | 16   | number of vertices in the base face
///  -------------------|:----:|---------------------------------------------------
///  FIELD1             |      |
///  subface index      | 16   | relative sub-face index for non-quads
///  num control points | 16   | number of control vertices in the 1-ring
///  -------------------|:----:|---------------------------------------------------
///  numPatchPoints     |      | number of patch points required for a given level 
///                     |      | of isolation of the surface
///  -------------------|:----:|---------------------------------------------------
///  padding            |      | bring size of the descriptor to 12 bytes 
/// 
struct TreeDescriptor {

    void Set(bool regularFace, uint16_t faceSize, uint16_t subfaceIndex, uint16_t numControlPoints);

    constexpr bool IsRegularFace() const { return unpack(field0, 1, 0) != 0; }

    constexpr uint32_t GetFaceSize() const { return unpack(field0, 16, 16); }

    constexpr uint32_t GetSubfaceIndex() const { return unpack(field1, 16, 0); }

    constexpr uint32_t GetNumControlPoints() const { return unpack(field1, 16, 16); }

    constexpr uint32_t GetNumPatchPoints(int level) const { return numPatchPoints[level]; }

    uint32_t field0 = 0;
    uint32_t field1 = 0;
    uint32_t numPatchPoints[kMaxIsolationLevel + 1] = { 0 };
    uint32_t padding = 0;
};

inline void TreeDescriptor::Set(bool regularFace, uint16_t faceSize, uint16_t subfaceIndex, uint16_t numControlPoints) {    
    field0 = pack(faceSize, 16, 16) |
             // pack(padding, 15, 1) |
             pack(regularFace, 1, 0);
    field1 = pack(numControlPoints, 16, 16) |
             pack(subfaceIndex, 16, 0);
}

#ifdef __cplusplus
} // end namespace Tmr
} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
#endif

#endif // OPENSUBDIV3_TMR_NODE_BASE_H
