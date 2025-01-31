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

#ifndef OPENSUBDIV3_TMR_NODE_DESCRIPTOR_H
#define OPENSUBDIV3_TMR_NODE_DESCRIPTOR_H

#ifdef __cplusplus
#include "../tmr/types.h"
#include <cstdint>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {
namespace Tmr {
#endif

enum class NodeType : uint16_t {
    NODE_REGULAR = 0,
    NODE_RECURSIVE = 1,
    NODE_TERMINAL = 2,
    NODE_END = 3,
};

enum class EndCapType : uint16_t {
    ENDCAP_NONE = 0,             ///< no endcap
    ENDCAP_BILINEAR_BASIS,       ///< use bilinear quads (4 cp) as end-caps
    ENDCAP_BSPLINE_BASIS,        ///< use BSpline basis patches (16 cp) as end-caps
    ENDCAP_GREGORY_BASIS,        ///< use Gregory basis patches (20 cp) as end-caps
};

struct NodeDescriptor {

    /// \brief Set bitfields for REGULAR nodes
    ///
    ///  Field         | Bits | Content
    ///  --------------|:----:|---------------------------------------------------
    ///  type          | 2    | NodeType
    ///  single crease | 1    | Whether the patch is of "single crease" type
    ///  depth         | 4    | level of isolation of the patch
    ///  boundary      | 5    | boundary edge mask encoding
    ///  v             | 10   | log2 value of u parameter at first patch corner
    ///  u             | 10   | log2 value of v parameter at first patch corner
    void SetRegular(bool singleCrease, uint8_t depth, uint8_t boundary, short u, short v);

    /// \brief Set bitfields for END nodes
    ///
    ///  Field         | Bits | Content
    ///  --------------|:----:|---------------------------------------------------
    ///  type          | 2    | NodeType
    ///  (unused)      | 1    | 
    ///  depth         | 4    | level of isolation of the patch
    ///  boundary      | 5    | boundary edge mask encoding
    ///  v             | 10   | log2 value of u parameter at first patch corner
    ///  u             | 10   | log2 value of v parameter at first patch corner
    void SetEnd(uint8_t depth, uint8_t boundary, short u, short v);

    /// \brief Set bitfields for RECURSIVE nodes
    ///
    ///  Field         | Bits | Content
    ///  --------------|:----:|---------------------------------------------------
    ///  type          | 2    | NodeType
    ///  has end-cap   | 1    | whether the patch has an end-cap
    ///  depth         | 4    | level of isolation of the patches
    ///  (unused)      | 5    | 
    ///  v             | 10   | log2 value of u parameter at first patch corner
    ///  u             | 10   | log2 value of v parameter at first patch corner
    void SetRecursive(uint8_t depth, short u, short v, bool hasEndcap);

    /// \brief Set bitfields for TERMINAL nodes
    ///
    ///  Field         | Bits | Content
    ///  --------------|:----:|---------------------------------------------------
    ///  type          | 2    | NodeType
    ///  has end-cap   | 1    | whether the patch has an end-cap
    ///  depth         | 4    | level of isolation of the patches
    ///  (unused)      | 1    |  
    ///  evIndex       | 4    | local index of the extraordinary vertex
    ///  v             | 10   | log2 value of u parameter at first patch corner
    ///  u             | 10   | log2 value of v parameter at first patch corner
    void SetTerminal(uint8_t depth, uint8_t evIndex, short u, short v, bool hasEndcap);

    void Clear() { field0 = 0; }

    //
    // Generic accessors
    //

    // The following accessors decode bitfields shared by all the nodes.

    constexpr NodeType GetType() const { return NodeType(unpack(field0, 2, 0)); }

    /// \brief Returns the depth of the node in the tree, which corresponds to the
    /// isolation level of a sub-patch
    constexpr uint32_t GetDepth() const { return unpack(field0, 4, 3); }

    /// \brief Returns the log2 value of the u parameter at the top left corner of
    /// the patch
    constexpr uint32_t GetU() const { return unpack(field0, 10, 12); }

    /// \brief Returns the log2 value of the v parameter at the top left corner of
    /// the patch
    constexpr uint32_t GetV() const { return unpack(field0, 10, 22); }

    /// \brief Returns the fraction of normalized parametric space covered by the
    /// sub-patch.
    float GetParamFraction(bool regularFace) const;

    /// \brief Maps the (u,v) parameterization from coarse to refined
    /// The (u,v) pair is mapped from the coarse face parameterization to
    /// the refined face parameterization
    void MapCoarseToRefined(float& u, float& v, bool regularFace) const;

    /// \brief Maps the (u,v) parameterization from refined to coarse
    /// The (u,v) pair is mapped from the refined face parameterization to
    /// the coarse face parameterization
    void MapRefinedToCoarse(float& u, float& v, bool regularFace) const;

    //
    // Type-specific accessors
    // 

    // The following accessors decode bitfields that are specific to certain types
    // of nodes only. Behavior is otherwise 'undefined': proceed with care !

    /// \brief Returns the boundary edge encoding mask (see Far::PatchParam)
    /// (REGULAR and END node only)
    constexpr uint32_t GetBoundaryMask() const { return unpack(field0, 5, 7); }

    /// \brief Returns the number of boundary edges in the sub-patch (-1 for invalid mask)
    /// (REGULAR and END node only)
    uint32_t GetBoundaryCount() const;

    /// \brief Returns local index of the extraordinary vertex 
    /// (TERMINAL node only)
    constexpr uint32_t GetEvIndex() const { return unpack(field0, 4, 8); }

    /// \brief True if the node has a fall-back irregular patch that can be used for
    /// dynamic isolation 
    /// (TERMINAL or RECURSIVE node only)
    constexpr bool HasEndcap() const { return unpack(field0, 1, 2) != 0; }

    /// \brief Returns true if the node has a 'sharpness' value
    /// (REGULAR node only)
    constexpr bool HasSharpness() const { return unpack(field0, 1, 2) != 0; }

    uint32_t field0 = 0;
};

static_assert(sizeof(NodeDescriptor) == sizeof(int));

inline void
NodeDescriptor::SetRegular(bool singleCrease, uint8_t depth, uint8_t boundary, short u, short v) {
    field0 = pack(v, 10, 22) |
             pack(u, 10, 12) |
             pack(boundary, 5, 7) |
             pack(depth, 4, 3) |
             pack(singleCrease, 1, 2) |
             pack(uint8_t(NodeType::NODE_REGULAR), 2, 0);
}
inline void
NodeDescriptor::SetEnd(uint8_t depth, uint8_t boundary, short u, short v) {
    field0 = pack(v, 10, 22) |
             pack(u, 10, 12) |
             pack(boundary, 5, 7) |
             pack(depth, 4, 3) |
             // pack(unused, 1, 3);
             pack(uint8_t(NodeType::NODE_END), 2, 0);
}
inline void
NodeDescriptor::SetRecursive(uint8_t depth, short u, short v, bool hasEndcap) {
    field0 = pack(v, 10, 22) |
             pack(u, 10, 12) |
             // pack(unused, 5, 7);
             pack(depth, 4, 3) |
             pack(hasEndcap, 1, 2) |
             pack(uint8_t(NodeType::NODE_RECURSIVE), 2, 0);
}
inline void
NodeDescriptor::SetTerminal(uint8_t depth, uint8_t evIndex, short u, short v, bool hasEndcap) {
    field0 = pack(v, 10, 22) |
             pack(u, 10, 12) |
             pack(evIndex, 4, 8) |
             // pack(unused, 1, 7);
             pack(depth, 4, 3) |
             pack(hasEndcap, 1, 2) |
             pack(uint8_t(NodeType::NODE_TERMINAL), 2, 0);
}

inline uint32_t
NodeDescriptor::GetBoundaryCount() const {
    return countbits(GetBoundaryMask());
}

inline float
NodeDescriptor::GetParamFraction(bool regularFace) const {
    uint32_t depth = regularFace ? GetDepth() : GetDepth() - 1;
    return 1.0f / float(1 << depth);
}

inline void
NodeDescriptor::MapCoarseToRefined(float& u, float& v, bool regularFace) const {
    float frac = GetParamFraction(regularFace);
    float pu = (float)GetU() * frac;
    float pv = (float)GetV() * frac;
    u = (u - pu) / frac;
    v = (v - pv) / frac;
}

inline void
NodeDescriptor::MapRefinedToCoarse(float& u, float& v, bool regularFace) const {
    float frac = GetParamFraction(regularFace);
    float pu = (float)GetU() * frac;
    float pv = (float)GetV() * frac;
    u = u * frac + pu;
    v = v * frac + pv;
}

#ifdef __cplusplus
} // end namespace Tmr
} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
#endif

#endif // OPENSUBDIV3_TMR_NODE_DESCRIPTOR_H