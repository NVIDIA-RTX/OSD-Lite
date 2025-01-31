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


#ifndef OPENSUBDIV3_TMR_SURFACE_H
#define OPENSUBDIV3_TMR_SURFACE_H

#ifdef __cplusplus

#include "../version.h"

#include "../tmr/types.h"

#include <cstdint>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {
namespace Tmr {
#endif

enum class Domain : uint8_t {
    Tri = 0,
    Quad,
    Quad_Subface,
};

///
///  \brief Linear Surface descriptor
///
/// Specialized descriptor for linearly interpolated surfaces. Linear surfaces
/// do not require subdivision plans or other external data and can be evaluated
/// directly from the SurfaceTable.
/// 
/// (held by SurfaceTable)
/// 
/// Encoding:
/// 
///  field0        | Bits | Content
///  --------------|:----:|---------------------------------------------------
///  face size     | 16   | number of control points in the face
///  subface index | 16   | index of the quad sub-face (or invalid index = 0xFF)
/// 

struct LinearSurfaceDescriptor {

    void Set(unsigned int firstPoint, uint16_t faceSize, uint16_t quadSubface = ~uint16_t(0));

    void SetNoLimit() { field0 = 0; firstControlPoint = ~uint32_t(0); };
    bool HasLimit() const { return GetFaceSize() != 0; }

    constexpr uint16_t GetFaceSize() const { return uint16_t(unpack(field0, 16, 0)); }
    constexpr LocalIndex GetQuadSubfaceIndex() const { return LocalIndex(unpack(field0, 16, 16)); }

    Index GetPatchPoint(int pointIndex);
    static Index GetPatchPoint(int pointInex, uint16_t faceSize, LocalIndex subfaceIndex);

    static Domain getDomain(uint16_t faceSize, LocalIndex subfaceIndex);
    Domain GetDomain() const { return getDomain(GetFaceSize(), GetQuadSubfaceIndex()); }

    uint32_t field0 = 0;
    uint32_t firstControlPoint = 0;
};

inline void 
LinearSurfaceDescriptor::Set(
    unsigned int firstPoint, uint16_t faceSize, LocalIndex quadSubface) {

    field0 = pack(faceSize, 16, 0) |
             pack(quadSubface, 16, 16);
    firstControlPoint = firstPoint;
}

inline Index
LinearSurfaceDescriptor::GetPatchPoint(int pointIndex, uint16_t faceSize, LocalIndex subfaceIndex) {
    if (subfaceIndex == LOCAL_INDEX_INVALID) {
        assert(pointIndex < faceSize);
        return pointIndex;
    } else {
        assert(pointIndex < 4);
        // patch point indices layout (N = faceSize) :
        // [ N control points ] 
        // [ 1 face-point ] 
        // [ N edge-points ]
        int N = faceSize;
        switch (pointIndex) {
            case 0: return subfaceIndex;
            case 1: return N + 1 + subfaceIndex; // edge-point after
            case 2: return N; // central face-point
            case 3: return N + (subfaceIndex > 0 ? subfaceIndex : N);
        }
    }
    return INDEX_INVALID;
}

inline Index
LinearSurfaceDescriptor::GetPatchPoint(int pointIndex) {
    return GetPatchPoint(pointIndex, GetFaceSize(), GetQuadSubfaceIndex());
}

inline Domain 
LinearSurfaceDescriptor::getDomain(uint16_t faceSize, LocalIndex subfaceIndex) {
    if (subfaceIndex == LOCAL_INDEX_INVALID)
        return faceSize == 4 ? Domain::Quad : Domain::Tri;
    return Domain::Quad_Subface;
}

///
///  \brief Surface descriptor
///
/// Aggregates pointers into multiple sets of data that need to be assembled in
/// order to evaluate the limit surface for the face of a mesh:
/// 
///   - the indices of the 1-ring set of control points around the face
///   - a pointer to the SubdivisionPlan with all the topological information
///     (composed of an index to a TopologyMap, and the index of the Plan itself
///     within that map)
///   - a subset of flags affecting the evaluation of the surface.
/// 
/// (held by SurfaceTable)
/// 
/// Encoding:
/// 
///  field0              | Bits | Content
///  --------------------|:----:|---------------------------------------------------
///  has limit           | 1    | limit surface cannot be evaluated if false (implies
///                      |      | other fields are expected to be set to 0)
///  param rotation      | 2    | parametric rotation of the subdivision plan
///  edges adjacency     | 4    | per-edge bits set: true if one or more surfaces
///                      |      | adjacent to that edge are irregular (the edge is a
///                      |      | T-junction) ; always false if the surface is irregular
///  topology map        | 5    | index of topology map (optional)
///  plan index          | 20   | index of the plan within the topology map selected
///  

struct SurfaceDescriptor {

    static constexpr uint32_t const kMaxMapIndex = (1 << 5) - 1;
    static constexpr uint32_t const kMaxPlanIndex = (1 << 20) - 1;

    void SetNoLimit() { field0 = 0; firstControlPoint = ~uint32_t(0); };
    constexpr bool HasLimit() const { return unpack(field0, 1, 0); }

    void Set(unsigned int firstPoint, unsigned int planIndex, uint8_t rotation, uint8_t adjacency, unsigned int mapIndex = 0);

    constexpr uint8_t GetParametricRotation() const { return unpack(field0, 2, 1); }

    constexpr uint8_t GetEdgeAdjacencyBits() const { return unpack(field0, 4, 3); }
    constexpr bool GetEdgeAdjacencyBit(uint8_t edgeIndex) const { uint8_t edgebits = unpack(field0, 4, 3); return (edgebits >> edgeIndex) & 0x1; }

    unsigned int GetTopologyMapIndex() const { return unpack(field0, 5, 7); }
    constexpr unsigned int GetSubdivisionPlanIndex() const { return unpack(field0, 20, 12); }

    uint32_t field0 = 0;
    uint32_t firstControlPoint = 0;
};

inline void SurfaceDescriptor::Set(
    unsigned int firstPoint, unsigned int planIndex, uint8_t rotation, uint8_t adjacency, unsigned int mapIndex) {

    assert(planIndex < kMaxPlanIndex && rotation < 4 && mapIndex < kMaxMapIndex);

    field0 = pack(true, 1, 0) |
             pack(rotation, 2, 1) |
             pack(adjacency, 4, 3) |
             pack(mapIndex, 5, 7) |
             pack(planIndex, 20, 12);

    firstControlPoint = firstPoint;
}

#ifdef __cplusplus
} // end namespace Tmr
} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
#endif

#endif // OPENSUBDIV3_TMR_SURFACE_H
