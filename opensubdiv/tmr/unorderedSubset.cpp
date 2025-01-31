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

#include "../tmr/unorderedSubset.h"
#include "../vtr/level.h"
#include "../vtr/stackBuffer.h"

#include <algorithm>
#include <map>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

using namespace Vtr::internal;

namespace Tmr {

//  Specifying vertex and face-varying indices around a face-vertex
unsigned int getFaceVertexIncidentFaceVertexIndices(
    Level const& level, Index vIndex, Index indices[], int fvarChannel) {

    //Index vIndex = level.getFaceVertices(baseFace)[cornerVertex];

    ConstIndexArray vFaces = level.getVertexFaces(vIndex);
    ConstLocalIndexArray vInFace = level.getVertexFaceLocalIndices(vIndex);

    int nIndices = 0;
    for (int i = 0; i < vFaces.size(); ++i) {

        ConstIndexArray srcIndices = (fvarChannel < 0) ?
            level.getFaceVertices(vFaces[i]) : level.getFaceFVarValues(vFaces[i], fvarChannel);

        int srcStart = vInFace[i];
        int srcCount = srcIndices.size();
        for (int j = srcStart; j < srcCount; ++j)
            indices[nIndices++] = srcIndices[j];
        for (int j = 0; j < srcStart; ++j)
            indices[nIndices++] = srcIndices[j];
    }
    return nIndices;
}

//  Main and supporting internal datatypes and methods to connect unordered
//  faces and allow for topological traversals of the incident faces:
//
//  The fundamental element of this process is the following definition of
//  an Edge. It is lightweight and only stores a state (boundary, interior,
//  or non-manifold) along with the one or two faces for a manifold edge.
//  It is initialized as a boundary when first created and is then modified
//  by adding additional incident faces.

struct UnorderedSubset::Edge {
    
    Index endVertex = INDEX_INVALID;   
    Index prevFace = INDEX_INVALID;
    Index nextFace = INDEX_INVALID;

    union {
        struct {
            uint8_t boundary    : 1;
            uint8_t interior    : 1;
            uint8_t nonManifold : 1;
            uint8_t trailing    : 1;
            uint8_t degenerate  : 1;
            uint8_t duplicate   : 1;
            uint8_t infSharp    : 1;
            uint8_t semiSharp   : 1;
        };
        uint8_t flags = 0;
    };

    void setInterior() { boundary = false; interior = true; }
    void setNonManifold() { boundary = interior = false; nonManifold = true; }
    void setDegenerate() { setNonManifold(); degenerate = true; }
    void setDuplicate() { setNonManifold(); duplicate = true; }

    void setFace(Index face, bool t) { *((trailing = t) ? &prevFace : &nextFace) = face; }

    void addFace(Index face, bool t) {
        //  Update the state of the Edge based on the added incident face:
        if (boundary) {
            if (t == trailing) 
                setNonManifold(); //  Edge is reversed
            else if (face == (trailing ? prevFace : nextFace))
                setNonManifold(); //  Edge is repeated in the face
            else {
                setInterior();    //  Edge is manifold thus far -- promote to interior
                setFace(face, t);
            }
        } else if (interior) {
            setNonManifold();     //  More than two incident faces -- make non-manifold
        }
    }

    void setSharpness(float sharpness) {
        if (sharpness > 0.f) {
            if (Sdc::Crease::IsInfinite(sharpness)) 
                infSharp = true;
            else 
                semiSharp = true;
        }
    }
};

bool UnorderedSubset::VertexDescriptor::initialize(int numIncidentFaces) {
    //  Mark invalid if too many or too few incident faces specified:
    isValid = (numIncidentFaces > 0) && (numIncidentFaces <= Vtr::VALENCE_LIMIT);
    numFaces = isValid ? (uint16_t)numIncidentFaces : 0;
    isInitialized = isValid;
    return isInitialized;
}

bool UnorderedSubset::VertexDescriptor::finalize() {

    if (!isValid)
        return false;

    //  Test for valid face size assignments while converting the sizes
    //  to offsets. Also detect if the faces are all the same size -- in
    //  which case, ignore the explicit assignments:
    if (hasFaceSizes) {
        int  size0 = faceSizeOffsets[0];
        bool sameSizes = true;

        int sum = 0;
        for (int i = 0; i < numFaces; ++i) {
            int faceSize = faceSizeOffsets[i];
            if ((faceSize < 3) || (faceSize > Vtr::VALENCE_LIMIT)) {
                isValid = false;
                return false;
            }
            sameSizes &= (faceSize == size0);

            faceSizeOffsets[i] = sum;
            sum += faceSize;
        }
        faceSizeOffsets[numFaces] = sum;

        //  No need to make use of explicit face sizes and offsets:
        if (sameSizes)
            hasFaceSizes = false;
    }
    return (isFinallized = true);
}

void UnorderedSubset::VertexDescriptor::setIncidentFaceSize(int incFaceIndex, int faceSize) {
    if (!hasFaceSizes) {
        faceSizeOffsets.SetSize(numFaces + 1);
        std::fill(&faceSizeOffsets[0], &faceSizeOffsets[numFaces + 1], 0);
        hasFaceSizes = true;
    }
    faceSizeOffsets[incFaceIndex] = faceSize;
}

inline void UnorderedSubset::VertexDescriptor::setIncidentFaceEdgeSharpness(
    int faceIndex, float leadingEdgeSharpness, float trailingEdgeSharpness) {

    if (!hasEdgeSharpness) {
        faceEdgeSharpness.SetSize(numFaces * 2);
        std::fill(&faceEdgeSharpness[0], &faceEdgeSharpness[numFaces * 2], 0.0f);
        hasEdgeSharpness = true;
    }
    faceEdgeSharpness[2 * faceIndex] = leadingEdgeSharpness;
    faceEdgeSharpness[2 * faceIndex + 1] = trailingEdgeSharpness;
}

inline int UnorderedSubset::VertexDescriptor::getFaceSize(int face) const {
    return commonFaceSize ? commonFaceSize : (faceSizeOffsets[face + 1] - faceSizeOffsets[face]);
}

inline int UnorderedSubset::GetFaceAfter(int step) const {
    assert(step >= 0);
    if (step == 1) {
        return getConnectedFaceNext(_faceInRing);
    } else if (step == 2) {
        return getConnectedFaceNext(getConnectedFaceNext(_faceInRing));
    } else {
        int face = _faceInRing;
        for (; step > 0; --step)
            face = getConnectedFaceNext(face);
        return face;
    }
}

inline int UnorderedSubset::GetFaceBefore(int step) const {
    assert(step >= 0);
    if (step == 1) {
        return getConnectedFacePrev(_faceInRing);
    } else if (step == 2) {
        return getConnectedFacePrev(getConnectedFacePrev(_faceInRing));
    } else {
        int face = _faceInRing;
        for (; step > 0; --step)
            face = getConnectedFacePrev(face);
        return face;
    }
}

int UnorderedSubset::createUnOrderedEdges(Edge* edges, short* feEdges, Index const* fvIndices) const {

    // Optional map to help construction for high valence:
    std::map<Index, int> edgeMap;
    bool useMap = (_vdesc.numFaces > 16);

    // Iterate through the face-edge pairs to find connecting edges:
    Index vCorner = _vdesc.getFaceIndexAtCorner(0, fvIndices);

    int numFaceEdges = 2 * _vdesc.numFaces;
    int numEdges = 0;

    //  Don't rely on the tag yet to determine presence of sharpness:
    bool hasSharpness = _vdesc.hasEdgeSharpness;

    for (int feIndex = 0; feIndex < numFaceEdges; ++feIndex) {
        Index vIndex = (feIndex & 1) ?
            _vdesc.getFaceIndexTrailing((feIndex >> 1), fvIndices) :
            _vdesc.getFaceIndexLeading((feIndex >> 1), fvIndices);

        int eIndex = -1;
        if (vIndex != vCorner) {
            if (useMap) {
                if (auto eFound = edgeMap.find(vIndex); eFound != edgeMap.end())
                    eIndex = eFound->second;
                else
                    edgeMap[vIndex] = numEdges;
            } else {
                for (int j = 0; j < numEdges; ++j)
                    if (edges[j].endVertex == vIndex) {
                        eIndex = j;
                        break;
                    }
            }

            // Update an existing edge or create a new one
            if (eIndex >= 0) {
                edges[eIndex].addFace(feIndex >> 1, feIndex & 1);
            } else {
                // Index of the new (pre-allocated, but not initialized) edge:
                eIndex = numEdges++;
                Edge& e = edges[eIndex] = {};
                e.endVertex = vIndex;
                e.boundary = true;
                e.setFace(feIndex >> 1, feIndex & 1);
                if (hasSharpness)
                    e.setSharpness(_vdesc.faceEdgeSharpness[feIndex]);
            }
        } else {
            //  If degenerate, create unique edge (non-manifold)
            eIndex = numEdges++;
            edges[eIndex].endVertex = vIndex;
            edges[eIndex].setDegenerate();
        }
        assert(eIndex >= 0);
        feEdges[feIndex] = (short)eIndex;
    }
    return numEdges;
}


void UnorderedSubset::markDuplicateEdges(Edge* edges, short* feEdges, Index const* fvIndices) const {

    //  The edge assignment thus far does not correctly detect the presence
    //  of all edges repeated or duplicated in the same face, e.g. for quad
    //  with vertices {A, B, A, C} the edge AB occurs both as AB and BA.
    //  When the face is oriented relative to corner B, we have {B, A, C, A}
    //  and edge BA will be detected as non-manifold -- but not from corner
    //  A or C.
    //
    //  So look for repeated instances of the corner vertex in the face and
    //  inspect its neighbors to see if they match the leading or trailing
    //  edges.
    //
    //  This is a trivial test for a quad:  if the opposite vertex matches
    //  the corner vertex, both the leading and trailing edges will be
    //  duplicated and so can immediately be marked non-manifold.  So deal
    //  with the common case of all neighboring quads separately.
    if (_vdesc.commonFaceSize == 3)
        return;

    Index vCorner = fvIndices[0];
    int numFaces = _vdesc.numFaces;

    if (_vdesc.commonFaceSize == 4) {
        Index const* fvOpposite = fvIndices + 2;
        for (int face = 0; face < numFaces; ++face, fvOpposite += 4) {
            if (*fvOpposite == vCorner) {
                edges[feEdges[2 * face]].setDuplicate();
                edges[feEdges[2 * face + 1]].setDuplicate();
            }
        }
    } else {
        Index const* fv = fvIndices;

        for (int face = 0; face < numFaces; ++face) {
            int faceSize = _vdesc.getFaceSize(face);

            if (faceSize == 4) {
                if (fv[2] == vCorner) {
                    edges[feEdges[2 * face]].setDuplicate();
                    edges[feEdges[2 * face + 1]].setDuplicate();
                }
            } else {
                for (int j = 2; j < (faceSize - 2); ++j) {
                    if (fv[j] == vCorner) {
                        if (fv[j - 1] == fv[1])
                            edges[feEdges[2 * face]].setDuplicate();
                        if (fv[j + 1] == fv[faceSize - 1])
                            edges[feEdges[2 * face + 1]].setDuplicate();
                    }
                }
            }
            fv += faceSize;
        }
    }
}

void UnorderedSubset::assignUnOrderedFaceNeighbors(Edge* edges, short* feEdges) {

    int numFaceEdges = 2 * _vdesc.numFaces;

    for (int i = 0; i < numFaceEdges; ++i) {
        assert(feEdges[i] >= 0);

        Edge const& E = edges[feEdges[i]];
        bool edgeIsSingular = E.nonManifold || E.boundary;
        if (edgeIsSingular)
            _faceEdgeNeighbors[i] = -1;
        else
            _faceEdgeNeighbors[i] = (i & 1) ? E.nextFace : E.prevFace;
    }
}

void UnorderedSubset::finalizeUnOrderedTags(Edge const* edges, int numEdges) {

    //  Summarize properties of the corner given the number and nature of
    //  the edges around its vertex and initialize remaining members or
    //  tags that depend on them.
    //
    //  First, take inventory of relevant properties from the edges:
    int numNonManifoldEdges = 0;
    int numInfSharpEdges = 0;
    int numSemiSharpEdges = 0;
    int numSingularEdges = 0;

    bool hasBoundaryEdges = false;
    bool hasBoundaryEdgesNotSharp = false;
    bool hasDegenerateEdges = false;
    bool hasDuplicateEdges = false;

    for (int i = 0; i < numEdges; ++i) {
        Edge const& E = edges[i];

        if (E.interior) {
            numInfSharpEdges += E.infSharp;
            numSemiSharpEdges += E.semiSharp;
        } else if (E.boundary) {
            hasBoundaryEdges = true;
            hasBoundaryEdgesNotSharp |= !E.infSharp;
        } else {
            ++numNonManifoldEdges;
            hasDegenerateEdges |= E.degenerate;
            hasDuplicateEdges |= E.duplicate;
        }
        // Singular edges include all that are effectively inf-sharp:
        numSingularEdges += E.nonManifold || E.boundary || E.infSharp;
    }

    bool isNonManifold = false;
    bool isNonManifoldCrease = false;

    // Next determine whether manifold or not.  Some obvious tests quickly
    // indicate if the corner is non-manifold, but ultimately it will be
    // necessary to traverse the faces to confirm that they form a single
    // connected set (e.g. two cones sharing their apex vertex may appear
    // manifold to this point but as two connected sets are non-manifold).

    if (numNonManifoldEdges) {
        isNonManifold = true;
        if (!hasDegenerateEdges && !hasDuplicateEdges && !hasBoundaryEdges)
            // Special crease case that avoids sharpening: two interior
            // non-manifold edges radiating more than two sets of faces:
            isNonManifoldCrease = (numNonManifoldEdges == 2) && (_vdesc.numFaces > numEdges);
    } else {
        //  Mismatch between number of incident faces and edges:
        isNonManifold = ((numEdges - _vdesc.numFaces) != (int)hasBoundaryEdges);
        if (!isNonManifold) {
            // If all faces are not connected, the set is non-manifold:
            int numFacesInSubset = findConnectedSubsetExtent();
            if (numFacesInSubset < _vdesc.numFaces) {
                isNonManifold = true;
            }
        }
    }

    // XXXX manuelk : much of this work is redundant, as we know the vertex to be non-manifold
    // leaving this logic in for debugging purposes (for now)
    assert(isNonManifold);

    // Assign tags and other members related to the inventory of edges
    // (boundary status is relevant if non-manifold as it can affect
    // the presence of the limit surface):

    _tag.nonManifoldVerts = isNonManifold;

    _tag.boundaryVerts = hasBoundaryEdges;
    _tag.boundaryNonSharp = hasBoundaryEdgesNotSharp;

    _tag.infSharpEdges = (numInfSharpEdges > 0);
    _tag.semiSharpEdges = (numSemiSharpEdges > 0);
    _tag.infSharpDarts = (numInfSharpEdges == 1) && !hasBoundaryEdges;

    //  Conditions effectively making the vertex sharp, include the usual
    //  excess of inf-sharp edges plus some non-manifold cases:
    if ((numSingularEdges > 2) || (isNonManifold && !isNonManifoldCrease))
        _isImpInfSharp = true;
    else if ((numSingularEdges + numSemiSharpEdges) > 2)
        _isImpSemiSharp = true;

    //  Mark the vertex inf-sharp if implicitly inf-sharp:
    if (!_isExpInfSharp && _isImpInfSharp) {
        _tag.infSharpVerts = true;
        _tag.semiSharpVerts = false;
    }
}


void UnorderedSubset::connectUnOrderedFaces(Index const* fvIndices) {

    //  There are two transient sets of data needed here:  a set of Edges
    //  that connect adjoining faces, and a set of indices (one for each
    //  of the 2*N face-edges) to identify the Edge for each face-edge.
    //
    //  IMPORTANT -- since these later edge indices are of the same type
    //  and size as the internal face-edge neighbors, we'll use that array
    //  to avoid a separate declaration (and possible allocation) and will
    //  update it in place later.

    int numFaceEdges = _vdesc.numFaces * 2;

    _faceEdgeNeighbors.SetSize(numFaceEdges);

    //  Allocate and populate the edges and indices referring to them.
    //  Initialization fails to detect some "duplicate" edges in a face,
    //  so post-process to catch these before continuing:
    StackBuffer<Edge, 32, true> edges(numFaceEdges);

    short* feEdges = &_faceEdgeNeighbors[0];

    int numEdges = createUnOrderedEdges(edges, feEdges, fvIndices);

    markDuplicateEdges(edges, feEdges, fvIndices);
   
    //  Use the connecting edges to assign neighboring faces (overwriting
    //  our edge indices) and finish initializing the tags retaining the
    //  properties of the corner:
    assignUnOrderedFaceNeighbors(edges, feEdges);

    finalizeUnOrderedTags(edges, numEdges);
}

int UnorderedSubset::finalizeVertexTag(int faceInVertex) {

    assert(_vdesc.isFinallized);

    _faceInRing = (short)faceInVertex;

    int numFaceVerts;
    if (!_vdesc.hasFaceSizes) 
        numFaceVerts = _vdesc.numFaces * _vdesc.commonFaceSize;
    else {
        // face sizes are available as differences between offsets:
        _vdesc.commonFaceSize = 0;
        numFaceVerts = _vdesc.faceSizeOffsets[_vdesc.numFaces];
    }
    
    _isExpInfSharp = Sdc::Crease::IsInfinite(_vdesc.vertSharpness);
    _isExpSemiSharp = Sdc::Crease::IsSemiSharp(_vdesc.vertSharpness);

    //  Initialize tags from VertexDescriptor and other members
    //
    //  Note that not all tags can be assigned at this point if the vertex
    //  is defined by a set of unordered faces. In such cases, the tags
    //  will be assigned later when the connectivity between incident faces
    //  is determined. Those that can be assigned regardless of ordering
    //  are set here -- splitting the assignment of those remaining between
    //  ordered and unordered cases.
    _tag.bits = 0;
    _tag.unCommonFaceSizes = _vdesc.hasFaceSizes;
    _tag.irregularFaceSizes = (_vdesc.commonFaceSize != _regFaceSize);
    _tag.infSharpVerts = _isExpInfSharp;
    _tag.semiSharpVerts = _isExpSemiSharp;

    return numFaceVerts;
}

float UnorderedSubset::getImplicitVertexSharpness() const {
    if (_isImpInfSharp)
        return Sdc::Crease::SHARPNESS_INFINITE;
    assert(_isImpSemiSharp);

    //  Since this will be applied at an inf-sharp crease, there will be
    //  two inf-sharp edges in addition to the semi-sharp, so we only
    //  need find the max of the semi-sharp edges and whatever explicit
    //  vertex sharpness may have been assigned.  Iterate through all
    //  faces and inspect the sharpness of each leading interior edge:
    float sharpness = _vdesc.vertSharpness;
    for (int i = 0; i < _vdesc.numFaces; ++i)
        if (GetFacePrev(i) >= 0)
            sharpness = std::max(sharpness, _vdesc.faceEdgeSharpness[2 * i]);
    return sharpness;
}

int UnorderedSubset::findConnectedSubsetExtent() const {
    int extent = 0, fStart = _faceInRing;
    for (int f = GetFaceNext(fStart); f >= 0; f = GetFaceNext(f)) {
        //  periodic -- return
        if (f == fStart)
            return extent;
        ++extent;
    }
    for (int f = GetFacePrev(fStart); f >= 0; f = GetFacePrev(f))
        ++extent;
    return extent;
}

int UnorderedSubset::populateFaceVertexDescriptor(
    Level const& level, Index faceIndex, Index vertIndex, int vertInface) {

    //  Identify the vertex index for the specified corner of the face
    //  and topology information related to it:

    _vdesc.commonFaceSize = level.getNumFaceVertices(faceIndex);

    Index vIndex = vertIndex;

    ConstIndexArray vFaces = level.getVertexFaces(vIndex);
    int nFaces = vFaces.size();

    ConstLocalIndexArray vInFace = level.getVertexFaceLocalIndices(vIndex);

    Level::VTag vTag = level.getVertexTag(vIndex);

    // Initialize, assign and finalize the vertex topology:
    // Note there is no need to check valence or face sizes with any
    // max here as TopologyRefiner construction excludes extreme cases.
    _vdesc.initialize(nFaces);

    if (vTag._incidIrregFace) {
        for (int i = 0; i < nFaces; ++i) {
            _vdesc.setIncidentFaceSize(i, level.getFaceVertices(vFaces[i]).size());
        }
    }

    if (vTag._semiSharp || vTag._infSharp)
        _vdesc.vertSharpness = level.getVertexSharpness(vIndex);

    if (vTag._semiSharpEdges || vTag._infSharpEdges) {
        // Must use face-edges and identify next/prev edges in face:
        for (int i = 0; i < nFaces; ++i) {

            ConstIndexArray fEdges = level.getFaceEdges(vFaces[i]);

            int eLeading = vInFace[i];
            int eTrailing = (eLeading ? eLeading : fEdges.size()) - 1;

            _vdesc.setIncidentFaceEdgeSharpness(i,
                level.getEdgeSharpness(fEdges[eLeading]),
                level.getEdgeSharpness(fEdges[eTrailing]));
        }
    }

    _vdesc.finalize();

    //  Return the index of the base face around the vertex:
    //  Remember that for some non-manifold cases the face may occur
    //  multiple times around this vertex, so make sure to identify the
    //  instance that matches the specified corner of the face.
    for (int i = 0; i < vFaces.size(); ++i) {
        if ((vFaces[i] == faceIndex) && (vInFace[i] == vertInface)) {
            return i;
        }
    }
    assert("Cannot identify face-vertex around non-manifold vertex." == 0);
    return -1;
}

UnorderedSubset::UnorderedSubset(Level const& level, 
    Index faceIndex, Index vertIndex, int vertInFace, int regFaceSize) :
        _faceIndex(faceIndex), _vertIndex(vertIndex), _vertInFace(vertInFace), _regFaceSize(regFaceSize) {

    int faceInRing = populateFaceVertexDescriptor(level, faceIndex, vertIndex, vertInFace);

    if (faceInRing >= 0) {

        int numFaceVerts = finalizeVertexTag(faceInRing);
        assert(numFaceVerts > 0);

        StackBuffer<Index, 64, true> cFaceVertIndices;
        cFaceVertIndices.SetSize(numFaceVerts);

        int numFaceVertices = getFaceVertexIncidentFaceVertexIndices(level, vertIndex, cFaceVertIndices);
        assert(numFaceVertices > 0);

        connectUnOrderedFaces(cFaceVertIndices);

        if (_tag.boundaryNonSharp) {
            _vdesc.isValid = false;
            return;
        }
    }
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
