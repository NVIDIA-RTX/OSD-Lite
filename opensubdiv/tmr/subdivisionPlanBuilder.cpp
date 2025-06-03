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

#include "../tmr/subdivisionPlanBuilder.h"
#include "../tmr/neighborhood.h"
#include "../tmr/topologyMap.h"

#include "../far/patchDescriptor.h"
#include "../far/patchBuilder.h"
#include "../far/primvarRefiner.h"
#include "../far/sparseMatrix.h"
#include "../far/topologyDescriptor.h"
#include "../far/topologyRefinerFactory.h"
#include "../far/types.h"

#include "../vtr/array.h"
#include "../vtr/level.h"
#include "../vtr/fvarLevel.h"
#include "../vtr/refinement.h"
#include "../vtr/stackBuffer.h"

#include <cassert>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

using Vtr::internal::Refinement;
using Vtr::internal::Level;

//
//  Helper functions:
//

// convert CCW winding to match bitwise ^= traversal
//
//  Sequential      ^ Bitwise     traversal pseudo-code
//  +---+---+       +---+---+     corner = 0
//  | 3 | 2 |       | 2 | 3 |     while (recursive)
//  +---+---+  ==>  +---+---+         if s>0.5 then corner^=1; s=1-s;
//  | 0 | 1 |       | 0 | 1 |         if t>0.5 then corner^=2; t=1-t;
//  +---+---+       +---+---+         s *= 2
//                                    t *= 2
static uint8_t
permuteWinding(LocalIndex i) {
   static int const permuteWinding[4] = { 0, 1, 3, 2 };
   return permuteWinding[i];
}

static void
offsetIndices(ConstIndexArray indices, Index offset, Index result[]) {
    if (offset) {
        for (int i = 0; i < indices.size(); ++i)
            result[i] = offset + indices[i];
    } else {
        std::memcpy(result, indices.begin(), indices.size() * sizeof(Index));
    }
}

// Terminal node helpers : copy indices from 16-wide bicubic basis into
// 25-wide terminal node. ('X' = extraordinary vertex)

// evIndex
//         0           1            2            3
//    X . . . .    . . . . X    + + + + .    . + + + +
//    . + + + +    + + + + .    + + + + .    . + + + +
//    . + + + +    + + + + .    + + + + .    . + + + +
//    . + + + +    + + + + .    + + + + .    . + + + +
//    . + + + +    + + + + .    . . . . X    X . . . .
inline void
copyDiagonalIndices(int evIndex, Index const* src, Index* dst) {
    static int offsets[4] = { 6, 5, 0, 1 };
    Index* cornerPtr = dst + offsets[evIndex];
    memcpy(cornerPtr + 0,  src + 0,  4 * sizeof(Index));
    memcpy(cornerPtr + 5,  src + 4,  4 * sizeof(Index));
    memcpy(cornerPtr + 10, src + 8,  4 * sizeof(Index));
    memcpy(cornerPtr + 15, src + 12, 4 * sizeof(Index));
}

// evIndex
//        0            1            2            3
//    X + + + +    + + + + X    . . . . .    . . . . .
//    . . . . .    . . . . .    . . . . .    . . . . .
//    . . . . .    . . . . .    . . . . .    . . . . .
//    . . . . .    . . . . .    . . . . .    . . . . .
//    . . . . .    . . . . .    + + + + X    X + + + +
inline void
copyRowIndices(int evIndex, Index const* src, Index* dst) {
    static int srcOffsets[4] = { 0, 0, 12, 12 },
               dstOffsets[4] = { 1, 0, 20, 21 };
    src += srcOffsets[evIndex];
    dst += dstOffsets[evIndex];
    memcpy(dst, src, 4 * sizeof(Index));
}

// evIndex
//        0            1            2            3
//    X . . . .    . . . . X    . . . . +    + . . . .
//    + . . . .    . . . . +    . . . . +    + . . . .
//    + . . . .    . . . . +    . . . . +    + . . . .
//    + . . . .    . . . . +    . . . . +    + . . . .
//    + . . . .    . . . . +    . . . . X    X . . . .
inline void
copyColIndices(int evIndex, Index const* src, Index* dst) {
    static int srcOffsets[4] = { 0, 3, 3, 0 },
               dstOffsets[4] = { 5, 9, 4, 0 };
    for (int i=0; i<4; ++i, src+=4, dst+=5) {
        *(dst+dstOffsets[evIndex]) = *(src+srcOffsets[evIndex]);
    }
}


inline Far::PatchBuilder::BasisType
convertEndCapType(EndCapType type) {
    switch (type) {
        using enum EndCapType;
        case ENDCAP_BILINEAR_BASIS : return Far::PatchBuilder::BASIS_LINEAR;
        case ENDCAP_BSPLINE_BASIS : return Far::PatchBuilder::BASIS_REGULAR;
        case ENDCAP_GREGORY_BASIS : return Far::PatchBuilder::BASIS_GREGORY;
        default :
            return Far::PatchBuilder::BASIS_UNSPECIFIED;
    }
}

static Far::PatchBuilder*
createPatchBuilder(Far::TopologyRefiner const& faceRefiner, EndCapType endcaps) {
    
    Far::PatchBuilder::Options opts;
    opts.regBasisType = Far::PatchBuilder::BASIS_REGULAR;
    opts.irregBasisType = convertEndCapType(endcaps);
    opts.fillMissingBoundaryPoints = true;
    opts.approxInfSharpWithSmooth = false;
    opts.approxSmoothCornerWithSharp = false;

    return Far::PatchBuilder::Create(faceRefiner, opts);
}

//
//  When accessing a "row" for a control point, the only non-zero
//  entry is that at the index, with a value of 1, so just store
//  that index so the StencilRows can combine it:
//
struct ControlRow {
    ControlRow(int index) : _index(index) { }
    ControlRow() { }

    ControlRow operator[] (int index) const {
        return ControlRow(index);
    }

    //  Members:
    int _index;
};

//  A "row" for each stencil is just our typical vector of variable
//  size that needs to support [].
//
//  For the first level, there are no source rows for the control
//  points so combine with the proxy ControlRow defined above.  All
//  other levels will accumulate StencilRows as weighted combinations
//  of other StencilRows.
//
//  WIP - consider combining StencilRows to exploit SSE/AVX vectorization
//      - we can (in future) easily guarantee both are 4-word aligned
//      - we can also pad the rows to a multiple of 4
//      - prefer writing the combination in a portable way that makes
//        use of auto-vectorization
//
template <typename REAL=float>
struct StencilRow {

    StencilRow() : _data(0), _size(0) { }
    StencilRow(REAL* data, int size) :
        _data(data), _size(size) { }
    StencilRow(REAL const* data, int size) :
        _data(const_cast<REAL*>(data)), _size(size) { }

    void Clear() {
        for (int i = 0; i < _size; ++i) {
            _data[i] = 0.0f;
        }
    }

    void AddWithWeight(ControlRow const& src, REAL weight) {
        assert(src._index >= 0);
        _data[src._index] += weight;
    }

    void AddWithWeight(StencilRow const& src, REAL weight) {
        assert(src._size == _size);
        //  Weights passed here by PrimvarRefiner should be non-zero
        //  WIP - see note on potential/future auto-vectorization above
        for (int i = 0; i < _size; ++i) {
            _data[i] += weight * src._data[i];
        }
    }

    StencilRow operator[](int index) const {
        return StencilRow(_data + index * _size, _size);
    }

    //  Members:
    REAL* _data;
    int   _size;
};


//
// SubdivisionPlan builder implementation
//

SubdivisionPlanBuilder::SubdivisionPlanBuilder() {

    // reserve memory for the proto-nodes store ; this progression curve should be fine
    // for most 'non-pathological' topologies

    static int const levelSizes[kNumLevels] = { 1, 4, 16, 64, 128, 128, 128, 128, 128, 128, 128 };

    for (int level = 0, levelVertOffset = 0; level < kNumLevels; ++level)
        _protoNodeStore[level].reserve(levelSizes[level]);
}

SubdivisionPlanBuilder::~SubdivisionPlanBuilder() { }

inline SubdivisionPlanBuilder::ProtoNode const&
SubdivisionPlanBuilder::getProtoNodeChild(ProtoNode const& pn, LocalIndex childIndex) const {
    return const_cast<SubdivisionPlanBuilder*>(this)->getProtoNodeChild(pn, childIndex);
}

inline SubdivisionPlanBuilder::ProtoNode&
SubdivisionPlanBuilder::getProtoNodeChild(ProtoNode const& pn, LocalIndex childIndex) {
    return _protoNodeStore[pn.levelIndex+1][pn.children[childIndex]];
}

void
SubdivisionPlanBuilder::setFaceTags(FaceTags& tags, int levelIndex, LocalIndex faceIndex) const {

    tags.clear();

    tags.hasPatch = _patchBuilder->IsFaceAPatch(levelIndex, faceIndex)
         && _patchBuilder->IsFaceALeaf(levelIndex, faceIndex);

    if (!tags.hasPatch)
        return;

    tags.isRegular = _patchBuilder->IsPatchRegular(levelIndex, faceIndex);
    if (tags.isRegular) {

        tags.boundaryMask =
            _patchBuilder->GetRegularPatchBoundaryMask(levelIndex, faceIndex);

        if (tags.boundaryMask == 0 && _faceRefiner->GetAdaptiveOptions().useSingleCreasePatch) {

            Far::PatchBuilder::SingleCreaseInfo info;
            if (_patchBuilder->IsRegularSingleCreasePatch(levelIndex, faceIndex, info)) {

                unsigned int isolationLevel = _faceRefiner->GetAdaptiveOptions().isolationLevel;

                float sharpness = std::min(info.creaseSharpness, float(isolationLevel - levelIndex));

                if (sharpness>0.f) {
                    tags.isSingleCrease = true;
                    tags.boundaryIndex = info.creaseEdgeInFace;
                    tags.sharpness = sharpness;
                }
            }
        }
    }
}

void
SubdivisionPlanBuilder::initializeProtoNodes(LocalIndex subfaceIndex) {

    assert(_faceRefiner && _patchBuilder);

    // recursive lambda traverses the topology of the neighborhood in the _faceRefiner
    // and gathers proto-nodes.

    auto traverse = [this](int levelIndex, LocalIndex faceIndex, auto&& traverse) -> int {

        using enum NodeType;

        int index = (int)_protoNodeStore[levelIndex].size();

        ProtoNode& pn = _protoNodeStore[levelIndex].emplace_back();

        pn.active = true;
        pn.levelIndex = levelIndex;
        pn.faceIndex = faceIndex;
        pn.hasIrregularPatch = false;

        setFaceTags(pn.faceTags, levelIndex, faceIndex);

        if (pn.faceTags.hasPatch) {
            if (pn.faceTags.isRegular)
                pn.nodeType = uint8_t(NODE_REGULAR);
            else {
                pn.hasIrregularPatch = true;
                pn.nodeType = uint8_t(NODE_END);
            }
            pn.numChildren = 0;
        } else {
            ConstIndexArray children = _faceRefiner->GetLevel(levelIndex).GetFaceChildFaces(faceIndex);
            for (int i=0; i<children.size(); ++i) 
                pn.children[i] = traverse(levelIndex+1, children[i], traverse);
            pn.numChildren = (int)children.size();
            pn.hasIrregularPatch = _useDynamicIsolation && levelIndex > 0;
            pn.nodeType = uint8_t(NODE_RECURSIVE);
        }
        return index;
    };

    if (_faceSize == _regularFaceSize) {
        assert(_faceRefiner->GetLevel(0).GetFaceVertices(0).size() == _regularFaceSize);
        traverse(0, 0, traverse);
    } else {
        assert(_faceRefiner->GetLevel(0).GetFaceVertices(0).size() != _regularFaceSize);
        assert(subfaceIndex < _faceRefiner->GetLevel(0).GetFaceVertices(0).size());
        traverse(1, subfaceIndex, traverse);
    }
}

unsigned int
SubdivisionPlanBuilder::finalizeProtoNodes() {

    using Node = SubdivisionPlan::Node;

    int numLevels = _faceRefiner->GetNumLevels();

    // total tree size = tree descriptor size + nodes sizes
    unsigned int treeSize = NodeBase::rootNodeOffset();

    _numPatchPointsTotal = 0;
    _numIrregularPatchesTotal = 0;

    for (short level=0; level < numLevels; ++level) {

        LevelOffsets& offsets = _levelOffsets[level];

        offsets = {};

        for (ProtoNode& pn : _protoNodeStore[level]) {

            if (!pn.active)
                continue;

            pn.treeOffset = treeSize;
            pn.firstPatchPoint = _numPatchPointsTotal;

            int numIrregularPatches = 0;

            switch (pn.getType()) {
                using enum NodeType;
                case NODE_REGULAR: _numPatchPointsTotal += _regularPatchSize; break;
                case NODE_END: ++numIrregularPatches; break;
                case NODE_TERMINAL: _numPatchPointsTotal += _terminalPatchSize;
                case NODE_RECURSIVE:
                    if (pn.hasIrregularPatch)
                        ++numIrregularPatches;
                    break;
                default:
                    assert(0);
            }

            _numPatchPointsTotal += (numIrregularPatches * _irregularPatchSize);

            treeSize += Node::getNodeSize(pn.getType(), pn.faceTags.isSingleCrease);

            offsets.numIrregularPatches += numIrregularPatches;  
        }
        _numIrregularPatchesTotal += offsets.numIrregularPatches;
    }

    // generate per-level indexing information for the stencil matrix

    Index indexBase = 0;
   
    for (int level = 0, pointOffset = 0; level < numLevels; ++level) {

        LevelOffsets& offsets = _levelOffsets[level];

        int numVertices = _faceRefiner->GetLevel(level).GetNumVertices();
        
        offsets.levelOffset = pointOffset;       
        pointOffset += numVertices;

        offsets.regularIndexBase = indexBase;       
        indexBase += numVertices;

        if (_reorderStencilMatrix) {
            offsets.irregularIndexBase = indexBase;
            indexBase += offsets.numIrregularPatches * _irregularPatchSize;
        } else
            offsets.irregularIndexBase = _numControlPoints + _numRefinedPoints;
    }
    return treeSize;
}

void
SubdivisionPlanBuilder::identifyTerminalNodes() {

    // returns true if the children of the proto-node meet the requirements of a
    // terminal node:
    //     - 3 of the children proto-nodes must be regular
    //     - 1 of the children proto-nodes is irregular with no boundaries & no creases
    auto isTerminal = [this](ProtoNode const& pn, LocalIndex& evIndex) -> bool {

        if (pn.numChildren != 4)
            return false;
  
        int regularCount = 0, 
            irregularIndex = 0;

        for (int i = 0; i < (int)pn.numChildren; ++i) {

            ProtoNode const& child = getProtoNodeChild(pn, i);

            if (child.faceTags.isRegular) {

                assert(child.faceTags.hasPatch);
                ++regularCount;

            } else {

                // trivial rejection for boundaries or creases
                if ((child.faceTags.boundaryCount > 0) || child.faceTags.isSingleCrease)
                    return false;

                // complete check
                Level const& level = _faceRefiner->getLevel(pn.levelIndex);
                ConstIndexArray fverts = level.getFaceVertices(pn.faceIndex);
                assert(fverts.size() == _regularFaceSize);

                Level::VTag vt = level.getFaceCompositeVTag(fverts);
                if (vt._semiSharp || vt._semiSharpEdges || vt._rule != Sdc::Crease::RULE_SMOOTH)
                    return false;

                evIndex = i;
            }
        }
        return regularCount == 3 ? true : false;
    };

    unsigned int isolationLevel = _faceRefiner->GetAdaptiveOptions().isolationLevel;

    for (uint8_t level = 0; level < isolationLevel; ++level) {

        for (int pnIndex=0; pnIndex<(int)_protoNodeStore[level].size(); ++pnIndex) {

            ProtoNode& pn = _protoNodeStore[level][pnIndex];

            LocalIndex evIndex=LOCAL_INDEX_INVALID;

            if (isTerminal(pn, evIndex)) {

                assert(evIndex != LOCAL_INDEX_INVALID);

                // de-activate the regular child nodes that the terminal node replaces
                for (int i=0; i<(int)pn.numChildren; ++i) {
                    if (i == evIndex)
                        continue;
                    getProtoNodeChild(pn, i).active = false;
                }

                pn.evIndex = evIndex;
                pn.nodeType = uint8_t(NodeType::NODE_TERMINAL);
            }
        }
    }
}

bool
SubdivisionPlanBuilder::computeSubPatchDomain(
    int levelIndex, Index faceIndex, short& u, short& v) const {

    switch (_faceRefiner->GetSchemeType()) {
        case Sdc::SCHEME_CATMARK:
            return computeCatmarkSubPatchDomain(levelIndex, faceIndex, u, v);
        case Sdc::SCHEME_LOOP:
            return computeLoopSubPatchDomain(levelIndex, faceIndex, u, v);
        default:
            assert(0);
    }
    return false;
}

bool
SubdivisionPlanBuilder::computeCatmarkSubPatchDomain(
    int levelIndex, Index faceIndex, short& u, short& v) const {

    assert(_faceRefiner->GetSchemeType() == Sdc::SchemeType::SCHEME_CATMARK);

    int const regularFaceSize = 
        Sdc::SchemeTypeTraits::GetRegularFaceSize(Sdc::SchemeType::SCHEME_CATMARK);

    u = v = 0;

    Level const& level = _faceRefiner->getLevel(levelIndex);

    bool irregular = (level.getFaceVertices(faceIndex).size() != regularFaceSize);

    // Move up the hierarchy accumulating u,v indices to the coarse level:

    for (int i = levelIndex, ofs = 1; i > 0; --i) {

        Refinement const& refinement  = _faceRefiner->getRefinement(i-1);
        Level const&      parentLevel = _faceRefiner->getLevel(i-1);

        Index parentFaceIndex  = refinement.getChildFaceParentFace(faceIndex);

        int childIndexInParent = refinement.getChildFaceInParentFace(faceIndex);

        ConstIndexArray fverts = parentLevel.getFaceVertices(parentFaceIndex);

        if (fverts.size() == regularFaceSize) {
            switch (childIndexInParent) {
                case 0 :                     break; // CCW winding
                case 1 : { u+=ofs;         } break;
                case 2 : { u+=ofs; v+=ofs; } break;
                case 3 : {         v+=ofs; } break;
            }
            ofs = (unsigned short)(ofs << 1);
        } else {
            irregular = true;
        }
        faceIndex = parentFaceIndex;
    }
    assert(u < 1024 && v < 1024);
    return irregular;
}

bool
SubdivisionPlanBuilder::computeLoopSubPatchDomain(
    int levelIndex, Index faceIndex, short& u, short& v) const {

    assert(_faceRefiner->GetSchemeType() == Sdc::SchemeType::SCHEME_LOOP);

    int const regularFaceSize =
        Sdc::SchemeTypeTraits::GetRegularFaceSize(Sdc::SchemeType::SCHEME_LOOP);

    u = v = 0;

    Level const& level = _faceRefiner->getLevel(levelIndex);

    int ofs = 1;
    bool rotatedTriangle = false;

    // Move up the hierarchy accumulating u,v indices to the coarse level:

    // For triangle refinement, the parameterization is rotated at
    // the fourth triangle subface at each level. The u and v values
    // computed for rotated triangles will be negative while we are
    // walking through the refinement levels.

    // For now, we don't consider irregular faces for triangle refinement.

    for (int i = levelIndex; i > 0; --i) {

        Refinement const& refinement = _faceRefiner->getRefinement(i - 1);
        Level const& parentLevel = _faceRefiner->getLevel(i - 1);

        Index parentFaceIndex = refinement.getChildFaceParentFace(faceIndex);

        int childIndexInParent = refinement.getChildFaceInParentFace(faceIndex);

        ConstIndexArray fverts = parentLevel.getFaceVertices(parentFaceIndex);

        if (rotatedTriangle) {
            switch (childIndexInParent) {
                case 0:                         break;
                case 1: { u -= ofs;           } break;
                case 2: {           v -= ofs; } break;
                case 3: { u += ofs; v += ofs; rotatedTriangle = false; } break;
            }
        } else {
            switch (childIndexInParent) {
                case 0:                         break;
                case 1: { u += ofs;           } break;
                case 2: {           v += ofs; } break;
                case 3: { u -= ofs; v -= ofs; rotatedTriangle = true; } break;
            }
        }
        ofs = (unsigned short)(ofs << 1);
        faceIndex = parentFaceIndex;
    }

    if (rotatedTriangle) {
        // If the triangle is tagged as rotated at this point then the
        // computed u and v parameters will both be negative and we map
        // them onto positive values in the opposite diagonal of the
        // parameter space.
        u += ofs;
        v += ofs;
    }

    assert(u < 1024 && v < 1024);
    return false;
}

void SubdivisionPlanBuilder::gatherIrregularPatchPoints(int levelIndex, Index* patchPoints) {

    // Irregular patches have a unique set of 'local' patch points that are factorized from
    // the regular patch points - we need to inject new rows in the stencil matrix for these
    // 'local points' ; the number of rows is given by _irregularPatchSize (REGULAR_BASIS adds
    // 16 points, GREGORY_BASIS 20, ...). 
    // 
    // By default these rows are located at the bottom of the stencil matrix, unless the rows
    // are re-ordered to be placed right after the stencils of the regular points of a given
    // level.

    Index indexBase = _levelOffsets[levelIndex].irregularIndexBase;
    
    unsigned int& indexOffset = _levelOffsets[_reorderStencilMatrix ? levelIndex : 0].irregularIndexOffset;

    for (int i = 0; i < _irregularPatchSize; ++i)
        patchPoints[i] = (indexBase + indexOffset++);
}

void
SubdivisionPlanBuilder::encodeRegularNode(
    ProtoNode const& pn, uint32_t* tree, Index* patchPoints) {

    unsigned int level = pn.levelIndex;
    Index face = pn.faceIndex;
    int boundaryMask = pn.faceTags.boundaryMask;

    short u, v;
    bool irregRoot = computeSubPatchDomain(level, face, u, v);

    Index patchVerts[16];
    _patchBuilder->GetRegularPatchPoints(level, face, boundaryMask, patchVerts);

    int offset = _levelOffsets[level].regularIndexBase;
    offsetIndices({ patchVerts, _regularPatchSize }, offset, patchPoints + pn.firstPatchPoint);

    bool singleCrease = pn.faceTags.isSingleCrease;
    if (singleCrease)
        boundaryMask = 1 << pn.faceTags.boundaryIndex;

    tree += pn.treeOffset;

    NodeDescriptor* descr = reinterpret_cast<NodeDescriptor*>(tree);

    descr->SetRegular(singleCrease, level, boundaryMask, u, v);

    tree[1] = pn.firstPatchPoint;

    if (singleCrease)
        *((float *)&tree[2]) = pn.faceTags.sharpness;
}

void
SubdivisionPlanBuilder::encodeEndCapNode(
    ProtoNode const& pn, uint32_t* tree, Index* patchPoints) {

    assert(pn.hasIrregularPatch);

    unsigned int level = pn.levelIndex;
    Index face = pn.faceIndex;

    short u, v;
    bool irregRoot = computeSubPatchDomain(level, face, u, v);

    gatherIrregularPatchPoints(level, patchPoints + pn.firstPatchPoint);

    tree += pn.treeOffset;

    NodeDescriptor* descr = reinterpret_cast<NodeDescriptor*>(tree);

    descr->SetEnd(level, 0, u, v);

    tree[1] = pn.firstPatchPoint;
}

void
SubdivisionPlanBuilder::encodeRecursiveNode(
    ProtoNode const& pn, uint32_t* tree, Index* patchPoints) {

    tree += pn.treeOffset;

    NodeDescriptor* descr = reinterpret_cast<NodeDescriptor*>(tree);

    if (pn.hasIrregularPatch) {

        unsigned int level = pn.levelIndex;

        gatherIrregularPatchPoints(level, patchPoints + pn.firstPatchPoint);

        short u, v;
        bool irregRoot = computeSubPatchDomain(level, pn.faceIndex, u, v);

        descr->SetRecursive(level, u, v, pn.hasIrregularPatch);

        tree[1] = pn.firstPatchPoint;

    } else {

        descr->SetRecursive(pn.levelIndex, 0, 0, pn.hasIrregularPatch);

        tree[1] = INDEX_INVALID;
    }

    tree += 2;
    assert((int)pn.numChildren==4);
    for (int i=0; i<(int)pn.numChildren; ++i)
        tree[permuteWinding(i)] = getProtoNodeChild(pn, i).treeOffset;
}

void
SubdivisionPlanBuilder::encodeTerminalNode(
    ProtoNode const& pn, uint32_t* tree, Index* patchPoints) {

    // xxxx manuelk: right now terminal nodes are recursive : paper suggests
    // a single node packing 25 * level patch points. Considering memory gains
    // not worth it against impact on code readability. 

    assert(pn.evIndex!=INDEX_INVALID);

    int childLevelIndex = pn.levelIndex + 1;

    patchPoints += pn.firstPatchPoint;

    for (int i=0; i<(int)pn.numChildren; ++i) {

        // gather patch point indices for the 3 regular sub-patches

        ProtoNode const& childNode = getProtoNodeChild(pn, i);

        if (pn.evIndex!=i) {

            // regular patch children nodes should have been de-activated when
            // this node was identified as Terminal
            assert(childNode.active==false);

            // XXXX can we have a terminal node w/ boundaries ? non-0 boundary mask work ?
            Index localVerts[16], patchVerts[16];
            _patchBuilder->GetRegularPatchPoints(childNode.levelIndex, childNode.faceIndex, 0, localVerts);

            // copy data to node buffers
            int offset = _levelOffsets[childNode.levelIndex].regularIndexBase;
            offsetIndices({ localVerts, 16 }, offset, patchVerts);

            // merge non-overlapping indices into a 5x5 array
                   if (i == ((pn.evIndex+2)%4)) {
                copyDiagonalIndices(pn.evIndex, patchVerts, patchPoints);
            } else if (i == (5-pn.evIndex)%4) {
                copyRowIndices(pn.evIndex, patchVerts, patchPoints);
            } else if (i == (3-pn.evIndex)%4) {
                copyColIndices(pn.evIndex, patchVerts, patchPoints);
            } else {
                assert(0);
            }
        }
        // not a recursive traversal : the xordinary child proto-node will be
        // handled when its proto-node is processed
    }

    // set index of the patch point in the EV corner to INVALID
    static int emptyIndices[4] = {0, 4, 24, 20 };
    patchPoints[emptyIndices[pn.evIndex]] = INDEX_INVALID;


    if (pn.hasIrregularPatch)
        gatherIrregularPatchPoints(pn.levelIndex, patchPoints + NodeBase::catmarkTerminalPatchSize());

    short u, v;
    bool irregRoot = computeSubPatchDomain(childLevelIndex, getProtoNodeChild(pn, 0).faceIndex, u, v);

    tree += pn.treeOffset;

    NodeDescriptor* descr = reinterpret_cast<NodeDescriptor*>(tree);

    descr->SetTerminal(pn.levelIndex, permuteWinding(pn.evIndex), u, v, pn.hasIrregularPatch);

    tree[1] = pn.firstPatchPoint;
    tree[2] = getProtoNodeChild(pn, pn.evIndex).treeOffset;
}

void
SubdivisionPlanBuilder::initializePatches(SubdivisionPlan& plan) {

    plan._patchPoints.resize(_numPatchPointsTotal);

    int numLevels = _faceRefiner->GetNumLevels();

    // encode nodes in the tree & collect patch point indices
    for (int level = 0; level < numLevels; ++level) {

        for (int pnIndex=0; pnIndex<(int)_protoNodeStore[level].size(); ++pnIndex) {

            ProtoNode const& pn = _protoNodeStore[level][pnIndex];

            if (!pn.active)
                continue;

            switch (pn.getType()) {
                using enum NodeType;
                case NODE_REGULAR:
                    encodeRegularNode(pn, plan._tree.data(), plan._patchPoints.data()); break;
                case NODE_END:
                    encodeEndCapNode(pn, plan._tree.data(), plan._patchPoints.data()); break;
                case NODE_TERMINAL:
                    encodeTerminalNode(pn, plan._tree.data(), plan._patchPoints.data()); break;
                case NODE_RECURSIVE:
                    encodeRecursiveNode(pn, plan._tree.data(), plan._patchPoints.data()); break;
            }
        }
    }

    // tabulate some offsets & sums that will help patch encoding
    TreeDescriptor* treeDesc = const_cast<TreeDescriptor *>(&plan.GetTreeDescriptor());
    for (int level = 0; level < numLevels; ++level) {

        LevelOffsets const& offsets = _levelOffsets[level];

        size_t npoints = 0;
        if (_reorderStencilMatrix)
            npoints = offsets.irregularIndexBase +
                (offsets.numIrregularPatches * _irregularPatchSize) - _numControlPoints;
        else
            npoints = (_numRefinedPoints + (_numIrregularPatchesTotal * _irregularPatchSize));

        // XXXX mk: can we use uint16_t ?
        assert(npoints < std::numeric_limits<std::remove_reference<decltype(*TreeDescriptor::numPatchPoints)>::type>::max());

        treeDesc->numPatchPoints[level] = uint32_t(npoints);
    }

    for (int level = numLevels; level < (int)std::size(treeDesc->numPatchPoints); ++level)
        treeDesc->numPatchPoints[level] = treeDesc->numPatchPoints[numLevels - 1];
}

void
SubdivisionPlanBuilder::getIrregularPatchConversion(ProtoNode const& pn) {

    //  The topology of an irregular patch is determined by its four corners:
    Level::VSpan cornerSpans[4];
    _patchBuilder->GetIrregularPatchCornerSpans(pn.levelIndex, pn.faceIndex, cornerSpans);

    //  Compute the conversion matrix from refined/source points to the
    //  set of points local to this patch:
    _patchBuilder->GetIrregularPatchConversionMatrix(pn.levelIndex, pn.faceIndex, cornerSpans, _conversion.matrix);

    //  Identify the refined/source points for the patch and append stencils
    //  for the local patch points in terms of the source points:
    int numSourcePoints = _conversion.matrix.GetNumColumns();

    _conversion.sourcePoints.resize(numSourcePoints);

    _patchBuilder->GetIrregularPatchSourcePoints(
        pn.levelIndex, pn.faceIndex, cornerSpans, _conversion.sourcePoints.data());

    int sourceIndexOffset = _levelOffsets[pn.levelIndex].levelOffset;
    for (int i = 0; i < numSourcePoints; ++i)
        _conversion.sourcePoints[i] += sourceIndexOffset;
}


void 
SubdivisionPlanBuilder::appendConversionStencilsToMatrix(
    std::vector<float>& stencilMatrix, int stencilBaseIndex) {

    //  Each row of the sparse conversion matrix corresponds to a row
    //  of the stencil matrix -- which will be computed from the weights
    //  and indices of stencils indicated by the SparseMatrix row:

    int numControlPoints = _numControlPoints;
    int numPatchPoints = _conversion.matrix.GetNumRows();

    StencilRow<float> srcStencils(&stencilMatrix[0], numControlPoints);
    StencilRow<float> dstStencils = srcStencils[stencilBaseIndex];

    for (int i = 0; i < numPatchPoints; ++i) {
        StencilRow<float> dstStencil = dstStencils[i];
        dstStencil.Clear();

        int const* rowIndices = &_conversion.matrix.GetRowColumns(i)[0];
        float const* rowWeights = &_conversion.matrix.GetRowElements(i)[0];
        int rowSize = _conversion.matrix.GetRowSize(i);

        for (int j = 0; j < rowSize; ++j) {
            float srcWeight = rowWeights[j];
            int srcIndex = _conversion.sourcePoints[rowIndices[j]];

            //  Simply increment single weight if this is a control point
            if (srcIndex < numControlPoints) {
                dstStencil._data[srcIndex] += srcWeight;
            } else {
                int srcStencilIndex = srcIndex - numControlPoints;

                StencilRow<float> srcStencil = srcStencils[srcStencilIndex];

                dstStencil.AddWithWeight(srcStencil, srcWeight);
            }
        }
    }
}

void
SubdivisionPlanBuilder::initializeStencilMatrix(SubdivisionPlan& plan) {

    int numPointStencils =
        _numRefinedPoints + (_numIrregularPatchesTotal * _irregularPatchSize);

    std::vector<float>& stencilMatrix = plan._stencilMatrix;

    stencilMatrix.resize(numPointStencils * _numControlPoints, 0.f);

    //  For refined points, initialize successive rows of the stencil matrix
    //  a level at a time using the PrimvarRefiner to accumulate contributing
    //  rows:
    int numLevels = _faceRefiner->GetNumLevels();

    if (numLevels > 1) {

        Far::PrimvarRefinerReal<float> primvarRefiner(*_faceRefiner);

        StencilRow<float> dstRow(&stencilMatrix[0], _numControlPoints);
        primvarRefiner.Interpolate(1, ControlRow(-1), dstRow);

        for (int level = 2; level < numLevels; ++level) {
            StencilRow<float> srcRow = dstRow;
            dstRow = srcRow[_faceRefiner->getLevel(level-1).getNumVertices()];
            primvarRefiner.Interpolate(level, srcRow, dstRow);
        }
    }

    //  For irregular patch points, append rows for each irregular patch:
    if (_numIrregularPatchesTotal > 0) {

        Index stencilIndexBase = _numRefinedPoints;

        for (int level = 0; level < numLevels; ++level) {

            if (_levelOffsets[level].numIrregularPatches == 0)
                continue;

            for (ProtoNode const& pn : _protoNodeStore[level]) {
               
                if (!pn.active || !pn.hasIrregularPatch)
                    continue;

                getIrregularPatchConversion(pn);

                appendConversionStencilsToMatrix(stencilMatrix, stencilIndexBase);

                stencilIndexBase += _irregularPatchSize;
            }
        }
    }
}

// Re-orders the stencil matrix rows: moves the irregular patch point stencil 
// rows from the end of the matrix, and re-inserts them at the correct row for 
// their isolation level ; the patchPoints indices are already set to match 
// this organization. This ordering allows for bulk serial evaluation of only 
// the stencils required for a given level of dynamic isolation.

void SubdivisionPlanBuilder::reorderStencilMatrix(SubdivisionPlan& plan) {

    assert(_levelOffsets[0].numIrregularPatches == 0);

    std::vector<float> const& inMatrix = plan._stencilMatrix;

    std::vector<float> outMatrix(inMatrix.size());

    int numLevels = _faceRefiner->GetNumLevels();
   
    int rowSize = _numControlPoints;
    int rowByteSize = rowSize * sizeof(decltype(outMatrix)::value_type);

    float const* srcRegularPtr = inMatrix.data();
    float const* srcIrregularPtr = inMatrix.data() + (_numRefinedPoints * rowSize);

    float* destPtr = outMatrix.data();

    for (int level = 1, srcOffset = 0; level < numLevels; ++level) {
    
        LevelOffsets const& offsets = _levelOffsets[level];

        int numRows = _faceRefiner->GetLevel(level).GetNumVertices();

        std::memcpy(destPtr, srcRegularPtr, numRows * rowByteSize);

        srcRegularPtr += numRows * rowSize;
        destPtr += numRows * rowSize;

        if (offsets.numIrregularPatches > 0) {

            numRows = offsets.numIrregularPatches * _irregularPatchSize;

            std::memcpy(destPtr, srcIrregularPtr, numRows * rowByteSize);

            srcIrregularPtr += numRows * rowSize;
            destPtr += numRows * rowSize;
        }
    }
    
    plan._stencilMatrix = std::move(outMatrix);
}

std::unique_ptr<SubdivisionPlan>
SubdivisionPlanBuilder::Create(Sdc::SchemeType scheme, Sdc::Options schemeOptions, 
    Options const& options, Neighborhood const& neighborhood, LocalIndex subfaceIndex) {

    for (auto& store : _protoNodeStore)
        store.clear();

    if (subfaceIndex == 0) {

        _faceSize = neighborhood.GetBaseFaceSize();

        _regularFaceSize = Sdc::SchemeTypeTraits::GetRegularFaceSize(scheme);

        uint8_t isolationLevel = std::min(options.isolationLevel, kMaxIsolationLevel);
         
        // force isolation level smooth feature <= sharp feature
        uint8_t isolationLevelSecondary = std::min(options.isolationLevelSecondary, isolationLevel);

        // force 1 level of isolation if the face is irregular
        if (_faceSize != _regularFaceSize) {
            isolationLevel = std::max(uint8_t(1), isolationLevel);
            isolationLevelSecondary = std::max(uint8_t(1), isolationLevelSecondary);
        }          

        bool useSingleCreasePatch = options.useSingleCreasePatch && (scheme != Sdc::SCHEME_LOOP);

        _faceRefiner = neighborhood.CreateRefiner(scheme, schemeOptions);

        Far::TopologyRefiner::AdaptiveOptions adaptiveOptions(isolationLevel);
        adaptiveOptions.secondaryLevel = isolationLevelSecondary;
        adaptiveOptions.useSingleCreasePatch = useSingleCreasePatch;
        adaptiveOptions.useInfSharpPatch = options.useInfSharpPatch;

        Index faceAtRoot = 0;
        ConstIndexArray baseFaceArray(&faceAtRoot, 1);

        _faceRefiner->RefineAdaptive(adaptiveOptions, baseFaceArray);

        _numControlPoints = _faceRefiner->GetLevel(0).GetNumVertices();
        _numRefinedPoints = _faceRefiner->GetNumVerticesTotal() - _numControlPoints;

        _patchBuilder.reset(
            createPatchBuilder(*_faceRefiner, options.endCapType));

        _regularPatchSize = Far::PatchDescriptor::GetNumControlVertices(_patchBuilder->GetRegularPatchType());
        _irregularPatchSize = Far::PatchDescriptor::GetNumControlVertices(_patchBuilder->GetIrregularPatchType());

        _useDynamicIsolation = options.useDynamicIsolation;
        _reorderStencilMatrix = options.orderStencilMatrixByLevel && options.useDynamicIsolation;

    } else
        // sub-face of an irregular face can re-use the builder state, but make sure the basic
        // details match at least, in case of horrible user-error
        assert(_faceSize == neighborhood.GetBaseFaceSize() && _faceSize != _regularFaceSize);

    assert(_patchBuilder && _faceRefiner);
    assert(neighborhood.GetBaseFaceSize() == _faceRefiner->GetLevel(0).GetFaceVertices(0).size());
    assert(neighborhood.GetNumControlPoints() == _numControlPoints);

    // Traverse refined TopologyRefiner & collect proto-nodes

    initializeProtoNodes(subfaceIndex);

    // Terminal (TER) nodes optimization is not supported w/ Loop scheme yet
    if (options.useTerminalNode && (scheme == Sdc::SCHEME_CATMARK))
        identifyTerminalNodes();

    int treeSize = finalizeProtoNodes();

    // Build the plan patches & stencils

    assert(treeSize > 0 && _numPatchPointsTotal > 0);

    std::unique_ptr<SubdivisionPlan> plan(new SubdivisionPlan(
        treeSize, _numControlPoints, _numPatchPointsTotal, _regularFaceSize, _faceSize, subfaceIndex));

    initializePatches(*plan);

    initializeStencilMatrix(*plan);

    if (_reorderStencilMatrix)
        reorderStencilMatrix(*plan);    

    return plan;
}

bool SubdivisionPlanBuilder::Options::operator == (Options const& other) const {
    return std::memcmp(this, &other, sizeof(SubdivisionPlanBuilder::Options)) == 0;
}


} // end namespace Tmr
} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv

