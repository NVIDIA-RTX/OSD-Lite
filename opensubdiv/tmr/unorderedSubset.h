//
//   Copyright 2015 Nvidia
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

#ifndef OPENSUBDIV3_TMR_UNOREDERED_SUBSET_H
#define OPENSUBDIV3_TMR_UNOREDERED_SUBSET_H

#include "../version.h"
#include "../tmr/types.h"
#include "../vtr/stackBuffer.h"
#include "../sdc/crease.h"

#include <cassert>
#include <cstdint>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Vtr::internal {
    class Level;
}

namespace Tmr {


unsigned int getFaceVertexIncidentFaceVertexIndices(
    Vtr::internal::Level const& level, Index vIndex, Index indices[], int fvarChannel = -1);

// \brief Unordered topology traversal helper
// 
// Neighborhood equivalence relies on the ordered traversal of the components within
// the 1-ring yielding bit-wise identical data-sets. However, Vtr cannot guarantee
// component ordering around non-manifold configurations, which results in false
// neighborhood equivalences and hashing collisions, leading to incorrect limit surface
// evaluation.
// 
// UnorderedSubset is a work-around for when a non-manifold vertex is detected in the
// 0-ring. Given a vertex in the 0-ring base-face, it attempts to build an ordered
// subset of the faces incident to the vertex. The subset is composed only of faces
// in the contiguous edge-cycle that contains the base face. Other discontinuous subsets
// are removed from the neighborhood, since they do not contribute to the limit surface.
// 
// This helper exposes a basic interface for an orderly traversal of the resulting subset.
// Assuming the subset is not empty, accessors to faces previous & next faces contiguous
// to the base face.
//
// note : this interface forces a bi-directional traversal
//   - first starting from the base face & going forward
//   - then rewinding to the base face & going backward
// We may need to be revisited if we still observe neighborhood collisions when hashing

class UnorderedSubset {

public:

    UnorderedSubset(Vtr::internal::Level const& level,
        Index faceIndex, Index vertIndex, int vertInFace, int regFaceSize);

    bool IsValid() const { return _vdesc.isValid && (_faceInRing != ~uint16_t(0)); }
    bool IsBoundary() const { return _tag.boundaryVerts; }

    int GetFaceStart() const { return _faceInRing; }

    inline int GetFaceNext(int face) const;
    inline int GetFacePrev(int face) const;

    inline int GetNumFacesAfter() const;
    inline int GetNumFacesBefore() const;

    inline int GetFaceAfter(int step) const;
    inline int GetFaceBefore(int step) const;

    inline bool FaceIndicesMatchAtCorner(int facePrev, int faceNext, Index const* indices) const;
    inline bool FaceIndicesMatchAtEdgeEnd(int facePrev, int faceNext, Index const* indices) const;
    inline bool FaceIndicesMatchAcrossEdge(int facePrev, int faceNext, Index const* indices) const;

    inline int GetNumFaceVerts() const;

private:

    union VertexTag {
        uint16_t bits = 0;
        struct {
            uint16_t boundaryVerts : 1;
            uint16_t infSharpVerts : 1;
            uint16_t infSharpEdges : 1;
            uint16_t infSharpDarts : 1;
            uint16_t semiSharpVerts : 1;
            uint16_t semiSharpEdges : 1;
            uint16_t unCommonFaceSizes : 1;
            uint16_t irregularFaceSizes : 1;
            uint16_t nonManifoldVerts : 1;
            uint16_t boundaryNonSharp : 1;
        };
        inline bool hasSharpVertices() const { return semiSharpVerts || infSharpVerts; }
        inline bool hasSharpEdges() const { return semiSharpEdges || infSharpEdges; }
    };

    struct VertexDescriptor {

        float vertSharpness = 0.f;

        uint16_t numFaces = 0;
        uint16_t commonFaceSize = 0;

        uint16_t isValid : 1 = 0;
        uint16_t isInitialized : 1 = 0;
        uint16_t isFinallized : 1 = 0;
        uint16_t hasFaceSizes : 1 = 0;
        uint16_t hasEdgeSharpness : 1 = 0;

        Vtr::internal::StackBuffer<int, 16, true> faceSizeOffsets;
        Vtr::internal::StackBuffer<float, 32, true> faceEdgeSharpness;

        bool initialize(int numIncidentFaces);
        bool finalize();

        void setIncidentFaceSize(int incFaceIndex, int faceSize);
        void setIncidentFaceEdgeSharpness(int faceIndex, float leadingEdgeSharpness, float trailingEdgeSharpness);

        int getFaceSize(int face) const;
        int getFaceIndexOffset(int face) const;

        Index getFaceIndexAtCorner(int face, Index const indices[]) const;
        Index getFaceIndexLeading(int face, Index const indices[]) const;
        Index getFaceIndexTrailing(int face, Index const indices[]) const;
    } _vdesc;

    Index _faceIndex = INDEX_INVALID;
    Index _vertIndex = INDEX_INVALID;
    int _vertInFace = -1;

    uint16_t _faceInRing = ~uint16_t(0);
    uint16_t _regFaceSize = 0;

    uint16_t _isExpInfSharp : 1 = 0;
    uint16_t _isExpSemiSharp : 1 = 0;
    uint16_t _isImpInfSharp : 1 = 0;
    uint16_t _isImpSemiSharp : 1 = 0;

    VertexTag _tag;

    Vtr::internal::StackBuffer<short, 32, true> _faceEdgeNeighbors;

    struct Edge;
    int createUnOrderedEdges(Edge* edges, short* faceEdgeIndices, Index const* faceVertIndices) const;
    void markDuplicateEdges(Edge* edges, short* faceEdgeIndices, Index const* faceVertIndices) const;
    void assignUnOrderedFaceNeighbors(Edge* edges, short* faceEdgeIndices);
    void finalizeUnOrderedTags(Edge const* edges, int numEdges);

    int getConnectedFaceNext(int face) const;
    int getConnectedFacePrev(int face) const;

    int populateFaceVertexDescriptor(Vtr::internal::Level const& level, Index faceIndex, Index vertIndex, int vertInface);
    float getImplicitVertexSharpness() const;
    int findConnectedSubsetExtent() const;
    int finalizeVertexTag(int faceInVertex);
    void connectUnOrderedFaces(Index const* fvIndices);
};

inline int UnorderedSubset::GetFaceNext(int face) const {
    return getConnectedFaceNext(face);
}
inline int UnorderedSubset::GetFacePrev(int face) const {
    return getConnectedFacePrev(face);
}

inline int UnorderedSubset::GetNumFacesAfter() const {
    assert(_vdesc.isFinallized);
    int fstart = GetFaceStart(), numFacesAfter = 0;
    for (int f = GetFaceNext(fstart); f >= 0; f = GetFaceNext(f), ++numFacesAfter)
        if (f == fstart)
            break;
    return numFacesAfter;
}

inline int UnorderedSubset::GetNumFacesBefore() const {
    assert(_vdesc.isFinallized);
    int fstart = GetFaceStart(), numFacesBefore = 0;
    for (int f = GetFacePrev(fstart); f >= 0; f = GetFacePrev(f), ++numFacesBefore)
        if (f == fstart)
            break;
    return numFacesBefore;
}

inline int UnorderedSubset::GetNumFaceVerts() const {
    assert(_vdesc.isFinallized);
    return _vdesc.hasFaceSizes ?
        _vdesc.faceSizeOffsets[_vdesc.numFaces] : _vdesc.numFaces * _vdesc.commonFaceSize;
}

inline bool UnorderedSubset::FaceIndicesMatchAtCorner(
    int facePrev, int faceNext, Index const* indices) const {
    assert(_vdesc.isFinallized);
    return _vdesc.getFaceIndexAtCorner(facePrev, indices) ==
        _vdesc.getFaceIndexAtCorner(faceNext, indices);
}

inline bool UnorderedSubset::FaceIndicesMatchAtEdgeEnd(
    int facePrev, int faceNext, Index const indices[]) const {
    assert(_vdesc.isFinallized);
    return _vdesc.getFaceIndexTrailing(facePrev, indices) ==
        _vdesc.getFaceIndexLeading(faceNext, indices);
}
inline bool UnorderedSubset::FaceIndicesMatchAcrossEdge(
    int facePrev, int faceNext, Index const indices[]) const {
    assert(_vdesc.isFinallized);
    return FaceIndicesMatchAtCorner(facePrev, faceNext, indices) &&
        FaceIndicesMatchAtEdgeEnd(facePrev, faceNext, indices);
}

inline int UnorderedSubset::getConnectedFaceNext(int face) const {
    return _faceEdgeNeighbors[2 * face + 1];
}
inline int UnorderedSubset::getConnectedFacePrev(int face) const {
    return _faceEdgeNeighbors[2 * face];
}

inline Index UnorderedSubset::VertexDescriptor::getFaceIndexAtCorner(int face, Index const indices[]) const {
    return indices[getFaceIndexOffset(face)];
}
inline Index UnorderedSubset::VertexDescriptor::getFaceIndexLeading(int face, Index const indices[]) const {
    return indices[getFaceIndexOffset(face) + 1];
}
inline Index UnorderedSubset::VertexDescriptor::getFaceIndexTrailing(int face, Index const indices[]) const {
    return indices[getFaceIndexOffset(face + 1) - 1];
}

inline int UnorderedSubset::VertexDescriptor::getFaceIndexOffset(int face) const {
    return commonFaceSize ? (face * commonFaceSize) : faceSizeOffsets[face];
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif // OPENSUBDIV3_TMR_UNOREDERED_SUBSET_H
