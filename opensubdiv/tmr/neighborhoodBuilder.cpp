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

#include "../tmr/neighborhoodBuilder.h"
#include "../tmr/unorderedSubset.h"
#include "../tmr/types.h"
#include "../far/topologyDescriptor.h"
#include "../vtr/fvarLevel.h"
#include "../vtr/level.h"

#include <limits>

// Re-orders creases & corners data in ascending order of vertex indices ; the arrays are
// assumed to be small, so a simple insertion sort should be optimal.
template <typename T> void sort(std::vector<T>& indices, std::vector<float>& sharpness) {
    
    assert(indices.size() == sharpness.size());
    
    int size = (int)indices.size();
    for (int i = 1; i < size; ++i) {
        int j = i - 1;        
        T key = indices[i];
        float s = sharpness[i];
        while (key < indices[j] && j > 0) {
            indices[j + 1] = indices[j];
            sharpness[j + 1] = sharpness[j];
            --j;
        }
        indices[j + 1] = key;
        sharpness[j + 1] = s;
    }
}

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

NeighborhoodBuilder::NeighborhoodBuilder(int maxValence) {

    reserve(maxValence);

    {
        using namespace Far;
        
        // Prevent accidental optimized builds w/ neighborhood debugging features still active
        // note : this guard-rail will be disabled if client code does not use CMake (or does
        // not define NDEBUG)
#ifndef NDEBUG
    #define NDEBUG 0
#endif
        static_assert((NDEBUG == 0) || (NDEBUG > 0 && kDebugNeighborhoods==false));


        // Make sure the Neighborhood vertex indexing types are binary-compatible with Far::TopologyDescriptor.
        // This allows to populate a local Far::TopologyRefiner directly from a Neighborhood without any
        // transcoding (at the cost of some memory though, since most of the neighborhood indices can be
        // represeneted with LocalIndices, with half the bits).
        static_assert(sizeof(Neighborhood) == (sizeof(uint32_t) * 4));

        static_assert(std::is_same_v<decltype(NeighborhoodBuilder::_faceVerts)::value_type, 
            std::remove_const_t<std::remove_pointer_t<decltype(TopologyDescriptor::vertIndicesPerFace)>>>);

        static_assert(std::is_same_v<decltype(TopologyDescriptor::vertIndicesPerFace), decltype(TopologyDescriptor::cornerVertexIndices)>);
        static_assert(std::is_same_v<decltype(TopologyDescriptor::vertIndicesPerFace), decltype(TopologyDescriptor::creaseVertexIndexPairs)>);
        static_assert(std::is_same_v<decltype(TopologyDescriptor::vertIndicesPerFace), decltype(TopologyDescriptor::holeIndices)>);
    }
}

void
NeighborhoodBuilder::clear() {

    _faces.clear();
    _faceVertCounts.clear();
    _faceVerts.clear();

    _creaseVerts.clear();
    _creaseSharpness.clear();
    
    _cornerVerts.clear();
    _cornerSharpness.clear();
    
    _vertRemaps.clear();
    _valueVerts.clear();
}

void 
NeighborhoodBuilder::reserve(int maxValence) {

    // First initialization guess ; does not have to be particularly
    // accurate as this scratch memory is persistent through re-uses
    // of the builder.
    // Neighborhoods with extremely high-valence vertices could cause
    // a fair amount of re-allocations from heavy push_back() use,
    // however those are discouraged and expected to be rare.

    int maxVerts = maxValence * 4;

    _faces.reserve(maxValence * 2);
    _faceVertCounts.reserve(maxValence * 2);
    _faceVerts.reserve(maxVerts);

    _cornerVerts.reserve(maxVerts);
    _cornerSharpness.reserve(maxVerts);

    _creaseVerts.reserve(maxVerts);
    _creaseVerts.reserve(maxVerts);

    _vertRemaps.reserve(maxVerts);
    _valueVerts.reserve(maxVerts);
}

// Remaps control point indices to an index-space local to the neighborhood.
// Assuming topology is always traversed in the same order, faces with
// identical topological configurations will have vertices with identical
// local space indices.

inline LocalIndex
NeighborhoodBuilder::findLocalIndex(Index vertexIndex) {
    for (LocalIndex i = 0; i < (LocalIndex)_vertRemaps.size(); ++i) {
        if (_vertRemaps[i] == vertexIndex)
            return i;
    }
    return INDEX_INVALID;
}

inline LocalIndex
NeighborhoodBuilder::mapToLocalIndex(Index vertexIndex) {
    for (LocalIndex i=0; i<(LocalIndex)_vertRemaps.size(); ++i) {
        if (_vertRemaps[i] == vertexIndex)
            return i;
    }
    _vertRemaps.push_back(vertexIndex);
    return (LocalIndex)_vertRemaps.size()-1;
}

inline LocalIndex
NeighborhoodBuilder::mapToLocalIndex(Index vertexIndex, Index valueIndex) {
    for (LocalIndex i=0; i<(LocalIndex)_vertRemaps.size(); ++i) {
        if (_vertRemaps[i] == valueIndex)
            return i;
    }
    _vertRemaps.push_back(valueIndex);
    _valueVerts.push_back(vertexIndex);
    return (LocalIndex)_vertRemaps.size()-1;
}

inline void
NeighborhoodBuilder::addFace(ConstIndexArray fverts, Index faceIndex, int startingEdge) {
    int nverts = fverts.size();
    _faceVertCounts.push_back(nverts);
    for (int vert = 0; vert < nverts; ++vert) {
        LocalIndex index = (vert + startingEdge) % nverts;
        _faceVerts.push_back(mapToLocalIndex(fverts[index]));
    }
    _faces.push_back(faceIndex);
}

inline void
NeighborhoodBuilder::addFace(ConstIndexArray fverts, ConstIndexArray fvalues, Index faceIndex, int startingEdge) {
    int nverts = fverts.size();
    _faceVertCounts.push_back(nverts);
    for (int vert = 0; vert < nverts; ++vert) {
        LocalIndex index = (vert + startingEdge) % nverts;
        _faceVerts.push_back(mapToLocalIndex(fverts[index], fvalues[index]));
    }
    _faces.push_back(faceIndex);
}

template <bool quash> inline void 
NeighborhoodBuilder::addCorner(Index v, float sharpness) {

    LocalIndex lv = findLocalIndex(v);
    assert(lv != LOCAL_INDEX_INVALID && sharpness > 0.f);

    int ncorners = (int)_cornerSharpness.size();
    for (int i = 0; i < ncorners; ++i)
        if (_cornerVerts[i] == lv) {
            if constexpr (quash) _cornerSharpness[i] = sharpness;
            return;
        }
    _cornerVerts.push_back(lv);
    _cornerSharpness.push_back(sharpness);
}

template <bool quash> inline void
NeighborhoodBuilder::addCrease(Index v0, Index v1, float sharpness) {

    LocalIndex lv0 = findLocalIndex(v0);
    LocalIndex lv1 = findLocalIndex(v1);
    assert(v0!=INDEX_INVALID && v1!= INDEX_INVALID && sharpness > 0.f);
    
    if (lv0 == LOCAL_INDEX_INVALID || lv1 == LOCAL_INDEX_INVALID)
        return;

    if (lv0 > lv1)
        std::swap(lv0, lv1);

    int ncreases = (int)_creaseSharpness.size();
    for (int i = 0; i < ncreases; ++i)
        if ((_creaseVerts[i].first == lv0) && (_creaseVerts[i].second == lv1)) {
            if constexpr (quash) _creaseSharpness[i] = sharpness;
            return;
        }
    _creaseVerts.push_back({ lv0, lv1 });
    _creaseSharpness.push_back(sharpness);
}


void NeighborhoodBuilder::addVertexFaces(Level const& level,
    Index baseFace, Index vertIndex, int vertInFace, int regFaceSize, bool unordered) {

    auto addVertexFace = [this, &level](Index faceIndex, Index vertIndex) {
        // skip faces that we have already indexed
        for (int i = 0; i < (int)_faces.size(); ++i)
            if (_faces[i] == faceIndex)
                return;

        ConstIndexArray fverts = level.getFaceVertices(faceIndex);
        int startingEdge = fverts.FindIndex(vertIndex);
        assert(startingEdge != INDEX_INVALID);
        this->addFace(fverts, faceIndex, startingEdge);
    };

    ConstIndexArray vfaces = level.getVertexFaces(vertIndex);

    if (!unordered) {
        
        // orient the enumeration starting from the 0-ring face
        Index firstFace = vfaces.FindIndex(baseFace);
        assert(firstFace != INDEX_INVALID);

        for (int face = 1; face < vfaces.size(); ++face) {
            Index findex = vfaces[(firstFace + face) % vfaces.size()];
            addVertexFace(findex, vertIndex);
        }
    
    } else {
        
        // need UnorderedSet helper to get an ordered traversal
        UnorderedSubset subset(level, baseFace, vertIndex, vertInFace, regFaceSize);
        assert(subset.IsValid());
        int fStart = subset.GetFaceStart();
        for (int f = subset.GetFaceNext(fStart); f >= 0; f = subset.GetFaceNext(f)) {
            if (f == fStart)
                goto done;
            addVertexFace(vfaces[f], vertIndex);
        }
        for (int f = subset.GetFacePrev(fStart); f >= 0; f = subset.GetFacePrev(f))
            addVertexFace(vfaces[f], vertIndex);
    done:        
        addCorner(vertIndex, Sdc::Crease::SHARPNESS_INFINITE);
    }
};

void 
NeighborhoodBuilder::gatherVertexTopology(BuildDescriptor const& desc) {

    Level const& level = desc.refiner.getLevel(0);

    Sdc::Options schemeOpts = desc.refiner.GetSchemeOptions();

    int regFaceSize = Sdc::SchemeTypeTraits::GetRegularFaceSize(desc.refiner.GetSchemeType());

    bool boundaryNone = schemeOpts.GetVtxBoundaryInterpolation() == Sdc::Options::VTX_BOUNDARY_NONE;
    bool sharpCorners = schemeOpts.GetVtxBoundaryInterpolation() == Sdc::Options::VTX_BOUNDARY_EDGE_AND_CORNER;

    auto addCrease = [this, &level, &boundaryNone](Index edgeIndex) {
        Level::ETag etag = level.getEdgeTag(edgeIndex);
        if (!boundaryNone && etag._boundary)
            return;
        float sharpness = etag._infSharp ?
            Sdc::Crease::SHARPNESS_INFINITE : level.getEdgeSharpness(edgeIndex);
        if (sharpness > 0.0f) {
            ConstIndexArray everts = level.getEdgeVertices(edgeIndex);
            this->addCrease(everts[0], everts[1], sharpness);
        }
    };

    clear();

    ConstIndexArray fverts = level.getFaceVertices(desc.faceIndex),
                    fedges = level.getFaceEdges(desc.faceIndex);
    
    assert(fverts.size() < VALENCE_LIMIT);

    // add 0-ring face to neighborhood
    addFace(fverts, desc.faceIndex, desc.startingEdge);

    for (int vert = 0; vert < fverts.size(); ++vert) {

        LocalIndex index = (vert + desc.startingEdge) % fverts.size();
        Index vertIndex = fverts[index];
        Index edgeIndex = fedges[index];

        Level::VTag vtag = level.getVertexTag(vertIndex);

        // add 0-ring edge sharpness (if any)
        if (vtag._rule > Sdc::Crease::RULE_SMOOTH)
            addCrease(edgeIndex);

        // add 1-ring faces & vertices around vert 'vertIndex'.
        bool unordered = vtag._nonManifold && (vtag._rule == Sdc::Crease::RULE_CORNER);
        addVertexFaces(level, desc.faceIndex, vertIndex, index, regFaceSize, unordered);

        // add 1-ring edge sharpness around vert 'vertIndex'
        // note: this excludes sharp edges along the 1-ring boundary as they
        // should not be relevant to the limit surface
        if (vtag._rule > Sdc::Crease::RULE_SMOOTH) {
            ConstIndexArray vedges = level.getVertexEdges(vertIndex);
            // order the enumeration starting from the 0-ring face
            Index firstEdge = vedges.FindIndex(edgeIndex);
            assert(firstEdge != INDEX_INVALID);
            for (int edge = 1; edge < vedges.size(); ++edge)
                addCrease(vedges[(firstEdge + edge) % vedges.size()]);
        }
    }

    // gather vertex sharpness values ('corners') if any
    // note: the 0-ring vertices are packed at the front of the remapping table.
    int numVertices = gatherOneRingCreaseAndCorners ? (int)_vertRemaps.size() : fverts.size();

    // drop corners if the 0-ring has degenerate vertex indices
    numVertices = std::min((int)_vertRemaps.size(), fverts.size());

    for (int vert = 0; vert < numVertices; ++vert) {
        Level::VTag vtag = level.getVertexTag(_vertRemaps[vert]);      
        if (vtag._nonManifold)
            continue;
        if (vtag._corner && sharpCorners)
            continue;
        float sharpness = vtag._infSharp ? 
            Sdc::Crease::SHARPNESS_INFINITE : level.getVertexSharpness(_vertRemaps[vert]);
        if (sharpness > 0.f) {
            _cornerVerts.push_back(vert);
            _cornerSharpness.push_back(sharpness);
        }
    }

    sort(_cornerVerts, _cornerSharpness);
    sort(_creaseVerts, _creaseSharpness);
}

inline bool
NeighborhoodBuilder::addVertexFace(
    FVarLevel const& fvlevel, Index faceIndex, Index vertIndex, Index valueIndex) {

    // skip faces that we have already indexed
    for (int i = 0; i < (int)_faces.size(); ++i)
        if (_faces[i] == faceIndex)
            return false;

    Vtr::internal::Level const& level = fvlevel.getLevel();

    ConstIndexArray fverts = level.getFaceVertices(faceIndex);
    ConstIndexArray fvalues = fvlevel.getFaceValues(faceIndex);

    // discard the face if its topology is not connected
    int vertInFace = fverts.FindIndex(vertIndex);
    assert(vertInFace != INDEX_INVALID);
    if (fvalues[vertInFace] != valueIndex)
        return false;

    addFace(fverts, fvalues, faceIndex, vertInFace);

    return true;
}

void 
NeighborhoodBuilder::gatherFVarTopology(BuildDescriptor const& desc) {

    assert(desc.fvarChannel >= 0 && desc.fvarChannel < desc.refiner.GetNumFVarChannels());

    clear();

    int regFaceSize = Sdc::SchemeTypeTraits::GetRegularFaceSize(desc.refiner.GetSchemeType());

    Sdc::SchemeType scheme = desc.refiner.GetSchemeType();

    Level const& level = desc.refiner.getLevel(0);
    FVarLevel const& fvlevel = level.getFVarLevel(desc.fvarChannel);

    Sdc::Options::FVarLinearInterpolation fvarOptions = fvlevel.getOptions().GetFVarLinearInterpolation();

    bool fvarCornersAreSharp = (fvarOptions > (scheme == Sdc::SCHEME_CATMARK ?
        Sdc::Options::FVAR_LINEAR_NONE : Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2));
    bool hasLinearBoundaries = (fvarOptions == Sdc::Options::FVAR_LINEAR_BOUNDARIES);

    ConstIndexArray fverts = level.getFaceVertices(desc.faceIndex);
    ConstIndexArray fedges = level.getFaceEdges(desc.faceIndex);
    ConstIndexArray fvalues = fvlevel.getFaceValues(desc.faceIndex);
    assert(fverts.size() < VALENCE_LIMIT);

    // add 0-ring face to neighborhood
    addFace(fverts, fvalues, desc.faceIndex, desc.startingEdge);

    Level::VTag::VTagSize vtagBits = 0;

    for (int vert = 0; vert < fverts.size(); ++vert) {

        LocalIndex index = (vert + desc.startingEdge) % fverts.size();

        Index vertIndex = fverts[index];
        //Index edgeIndex = fedges[index];
        Index valueIndex = fvalues[index];

        ConstIndexArray vfaces = level.getVertexFaces(vertIndex);
        Level::VTag vtag = level.getVertexTag(vertIndex);
        
        int fvindex = fvlevel.findVertexValueIndex(vertIndex, valueIndex);
        FVarLevel::ValueTag fvtag = fvlevel.getValueTag(fvindex);

        // orient the enumeration starting from the 0-ring face
        Index firstFace = vfaces.FindIndex(desc.faceIndex);
        assert(firstFace != INDEX_INVALID);

        // add 1-ring faces & vertices around vert 'vertIndex'.

        bool unordered = (vtag._nonManifold && (vtag._rule == Sdc::Crease::RULE_CORNER)) || fvtag._nonManifold;

        if (!unordered) {
            for (int face = 1; face < vfaces.size(); ++face) {
                Index findex = vfaces[(firstFace + face) % vfaces.size()];
                addVertexFace(fvlevel, findex, vertIndex, valueIndex);
            }
        } else {

            // need UnorderedSet helper to get an ordered traversal 
            // XXXX manuelk : ideally, would be nice to be able to re-use the 
            // subset from the vertex traversal instead of reconstructing here

            UnorderedSubset vtxSubset(level, desc.faceIndex, vertIndex, index, regFaceSize);
            assert(vtxSubset.IsValid());

            int numFaceVerts = vtxSubset.GetNumFaceVerts();

            Vtr::internal::StackBuffer<Index, 64, true> cFaceValueIndices;
            cFaceValueIndices.SetSize(numFaceVerts);

            int numFaceVertices = 
                getFaceVertexIncidentFaceVertexIndices(level, vertIndex, cFaceValueIndices, desc.fvarChannel);
            assert(numFaceVertices > 0);

            Index const* fvarIndices = &cFaceValueIndices[0];

            int cornerFace = vtxSubset.GetFaceStart();

            int numFacesAfter = vtxSubset.GetNumFacesAfter();

            bool done = false, boundary = true;
            if (numFacesAfter > 0) {
                int thisFace = cornerFace;
                int nextFace = vtxSubset.GetFaceNext(thisFace);
                for (int i = 0; i < numFacesAfter; ++i) {
                    if (!vtxSubset.FaceIndicesMatchAcrossEdge(thisFace, nextFace, fvarIndices))
                        break;
                    addVertexFace(fvlevel, vfaces[thisFace], vertIndex, valueIndex);
                    thisFace = nextFace;
                    nextFace = vtxSubset.GetFaceNext(thisFace);
                }
                if (nextFace == cornerFace) {
                    if (vtxSubset.FaceIndicesMatchAtEdgeEnd(thisFace, cornerFace, fvarIndices))
                        boundary = false;
                    done = true;
                }

            }

            if (!done) {
                int numFacesBefore = vtxSubset.GetNumFacesBefore();
                if (!vtxSubset.IsBoundary())
                    numFacesBefore += numFacesAfter - numFacesBefore;
                if (numFacesBefore > 0) {
                    int thisFace = cornerFace;
                    int prevFace = vtxSubset.GetFacePrev(thisFace);
                    for (int i = 0; i < numFacesBefore; ++i) {
                        if (!vtxSubset.FaceIndicesMatchAcrossEdge(prevFace, thisFace, fvarIndices))
                            break;
                        else
                            addVertexFace(fvlevel, vfaces[prevFace], vertIndex, valueIndex);
                        thisFace = prevFace;
                        prevFace = vtxSubset.GetFacePrev(thisFace);
                    }
                }
            }
            addCorner(valueIndex, Sdc::Crease::SHARPNESS_INFINITE);
        }

        if (fvarOptions == Sdc::Options::FVAR_LINEAR_CORNERS_PLUS1 ||
            fvarOptions == Sdc::Options::FVAR_LINEAR_CORNERS_PLUS2) {
            if (fvtag._mismatch) {
                if (fvtag._depSharp && fvtag._semiSharp) {
                    // dependent sharpness : see FVarLevel::completeTopologyFromFaceValues
                    vtag._infSharp = true;
                } else {
                    if (!fvtag._semiSharp)
                        vtag._infSharp |= !fvtag._crease;
                }
            }
        }

        if (vtag._semiSharp || vtag._infSharp) {
            if (fvarCornersAreSharp) {
                // if fvarCornersAreSharp, the assumption is that boundaries are implicitly
                // infinitely sharp and should not be tagged.
                if (!vtag._boundary) {
                    addCorner(valueIndex, vtag._infSharp ?
                        Sdc::Crease::SHARPNESS_INFINITE : level.getVertexSharpness(vertIndex));
                }
            } else {
                addCorner(valueIndex, vtag._infSharp ?
                    Sdc::Crease::SHARPNESS_INFINITE : level.getVertexSharpness(vertIndex));
            }
        }

        vtagBits |= vtag.getBits();
    }

    // gather sharpness tags for 1-ring edges
    bool traverseOneRing = true;  
    if (!hasLinearBoundaries && !gatherOneRingCreaseAndCorners) {
        Level::VTag compVTag(vtagBits);
        traverseOneRing = compVTag._semiSharpEdges || compVTag._infSharpEdges;
    }
    if (traverseOneRing) {

        for (int face = 0; face < (int)_faces.size(); ++face) {
        
            ConstIndexArray fedges = level.getFaceEdges(_faces[face]);
            ConstIndexArray fvalues = fvlevel.getFaceValues(_faces[face]);
        
            for (int edge = 0; edge < fedges.size(); ++edge) {

                Index edgeIndex = fedges[edge];

                // skip edges that are not connected to at least one vertex of the 0-ring face
                if (!gatherOneRingCreaseAndCorners) {
                    ConstIndexArray everts = level.getEdgeVertices(edgeIndex);
                    if (fverts.FindIndex(everts[0]) == INDEX_INVALID &&
                        fverts.FindIndex(everts[1]) == INDEX_INVALID) {
                        continue;
                    }
                }

                Level::ETag etag = level.getEdgeTag(edgeIndex);

                // if 'fvarCornersAreSharp', we assume that the topology map has vertex boundary
                // rules of type 'edge & corner', where boundary edges are implicitly sharp: in
                // case, adding a crease tag is redundant and bloats the topology map.
                if (fvarCornersAreSharp && etag._boundary)
                    continue;

                if (etag._semiSharp || etag._infSharp) {
                    float sharpness = etag._infSharp ? 
                        Sdc::Crease::SHARPNESS_INFINITE : level.getEdgeSharpness(fedges[edge]);
                    assert(sharpness > 0.f);
                    Index fv0 = fvalues[edge];
                    Index fv1 = fvalues[edge < (fedges.size() - 1) ? edge + 1 : 0];
                    addCrease(fv0, fv1, sharpness);
                }
            }
        }
    }

    // sharpen vertices on linear boundaries
    if (hasLinearBoundaries) {
        constexpr bool quash = true; // override semi-sharp tags if they already exist
        int nverts = (int)_vertRemaps.size();
        for (int vert = 0; vert < nverts; ++vert) {
            Index valueIndex = _vertRemaps[vert];
            FVarLevel::ValueTag vtag = fvlevel.getValueTag(fvlevel.findVertexValueIndex(_valueVerts[vert], valueIndex));
            if (vtag._infIrregular)
                addCorner<quash>(valueIndex, Sdc::Crease::SHARPNESS_INFINITE);
        }
    }

    sort(_cornerVerts, _cornerSharpness);
    sort(_creaseVerts, _creaseSharpness);
}

void
NeighborhoodBuilder::populateData(uint8_t* buf, size_t buf_size, int startingEdge, bool remapVerts) const {

    assert(!_faceVertCounts.empty() && !_faceVerts.empty());
    assert(_cornerVerts.size() == _cornerSharpness.size());
    assert(_creaseVerts.size() == _creaseSharpness.size());

    assert(_faceVertCounts.size() < std::numeric_limits<decltype(Neighborhood::_faceCount)>::max());
    assert(_faceVerts.size() < std::numeric_limits<decltype(Neighborhood::_faceVertsCount)>::max());
    assert(_cornerSharpness.size() < std::numeric_limits<decltype(Neighborhood::_cornersCount)>::max());
    assert(_creaseSharpness.size() < std::numeric_limits<decltype(Neighborhood::_creasesCount)>::max());
    assert(_vertRemaps.size() < std::numeric_limits<decltype(Neighborhood::_controlPointsCount)>::max());

    // rotate the starting edge to the opposite side
    int faceSize = _faceVertCounts[0];
    startingEdge = (faceSize - startingEdge) % faceSize;
    assert(startingEdge < std::numeric_limits<decltype(Neighborhood::_startingEdge)>::max());

    Neighborhood::populateData(buf, buf_size, {
        .faceVerts = _faceVerts.data(),
        .nfaceVerts = (uint16_t)_faceVerts.size(),
        .faceVertCounts = _faceVertCounts.data(),
        .nfaceVertCounts = (uint16_t)_faceVertCounts.size(),
        .cornerVerts = _cornerVerts.data(),
        .cornerSharpness = _cornerSharpness.data(),
        .ncorners = (uint16_t)_cornerVerts.size(),
        .creaseVerts = (int*)_creaseVerts.data(),
        .creaseSharpness = _creaseSharpness.data(),
        .ncreases = (uint16_t)_creaseSharpness.size(),
        .controlPoints = remapVerts ? _vertRemaps.data() : nullptr,
        .ncontrolPoints = remapVerts ? (uint16_t)_vertRemaps.size() : uint16_t(0),
        .startingEdge = (uint16_t)startingEdge,
    });
}

size_t
NeighborhoodBuilder::computeByteSize(bool remapVerts) const {
    return Neighborhood::byteSize(
        (int)_faceVertCounts.size(),
        (int)_faceVerts.size(),
        (int)_cornerSharpness.size(),
        (int)_creaseSharpness.size(),
        remapVerts ? (int)_vertRemaps.size() : 0);
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv

