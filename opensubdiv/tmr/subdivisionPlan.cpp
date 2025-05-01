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

#include "../tmr/subdivisionPlan.h"
#include "../tmr/neighborhoodBuilder.h"
#include "../tmr/subdivisionPlanBuilder.h"
#include "../tmr/topologyMap.h"

#include "../far/patchBasis.h"
#include "../far/patchParam.h"
#include "../far/patchDescriptor.h"

#include "../sdc/scheme.h"

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

//
 // The combination of dynamic isolation and single-crease patch optimizations
 // can introduce discontinuities in the limit surface along boundaries with
 // end-cap patches. The `SingleCreaseDynamicISolation` enum controls the
 // desired behavior of the limit surface as follows:
 // 
 //  - Smooth mode : local edge sharpness of single-crease patches is forcibly
 //    lowered in order to maintain continuity with the limit surface of
 //    neighboring end-cap patches. Effectively, the limit surface starts
 //    ignoring semi-sharp tags in areas where the level of feature isolation
 //    is dynamically lowered.
 // 
 //  - Sharp mode : edge sharpness of single-crease patches is always enforced,
 //    regardless of dynamic isolation. The consequence is that single-crease 
 //    patches neighboring an end-cap patch from a lower level of isolation will
 //    cause a discontinuity in the limit surface. Semi-sharp edge tags are
 //    maintained, but neighboring dynamically isolated patches will show
 //    discontinuities.
 //
 // note: dynamic isolation should always be performed at the micro-vertex
 // or edge-level. Inconsistent dynamic isolation levels across neighboring
 // surfaces will cause limit surface discontinuities.
 //
 // note: this option is currently set at compile-time and can be set with
 // the definition of TMR_SINGLE_CREASE_DYNAMIC_ISOLATION. This decision may be
 // revisited if valid use-cases can be made for a run-time selection.
 
 #ifndef TMR_SINGLE_CREASE_DYNAMIC_ISOLATION
     #define TMR_SINGLE_CREASE_DYNAMIC_ISOLATION SMOOTH
 #endif
 
 enum class SingleCreaseDynamicIsolation : uint8_t { SHARP = 0, SMOOTH = 1, };    
 
 static constexpr SingleCreaseDynamicIsolation const
     single_crease_dynamic_isolation = SingleCreaseDynamicIsolation::TMR_SINGLE_CREASE_DYNAMIC_ISOLATION;
 
 template <typename REAL> REAL
 computeSingleCreaseSharpness(
     SubdivisionPlan::Node const& n, NodeDescriptor const& desc, short depth, int level) {
 
     using enum SingleCreaseDynamicIsolation;
 
     if constexpr (single_crease_dynamic_isolation == SHARP) {
         return desc.HasSharpness() ? n.GetSharpness() : (REAL)0.0;
     } else if constexpr (single_crease_dynamic_isolation == SMOOTH) {
         if (desc.HasSharpness()) {
             REAL sharpness = n.GetSharpness();
             // single-crease patches require a non-null boundary mask and sharpness > 0.f
             // std::numeric_limits::min() ensures EvalBasisBSpline() evaluates the crease
             // matrix
             return std::max(std::numeric_limits<REAL>::min(),
                 + std::min(sharpness, (REAL)level - (REAL)depth));
         }
         return REAL(0);
     } else
         static_assert(false, "Invalid single crease dynamic isolation mode");
 }

//
// Subdivision Plan Node
//

int
SubdivisionPlan::Node::GetPatchSize(int quadrant, unsigned short maxLevel) const {

    Sdc::SchemeType scheme = GetSubdivisionPlan()->getSchemeType();

    int regularPatchSize = 0;
    int irregularPatchSize = 0;

    if (scheme == Sdc::SCHEME_CATMARK) {
        regularPatchSize = catmarkRegularPatchSize();
        irregularPatchSize = catmarkIrregularPatchSize(plan->getEndCapType());    
    } else if (scheme == Sdc::SCHEME_LOOP) {
        regularPatchSize = loopRegularPatchSize();
        irregularPatchSize = loopIrregularPatchSize(plan->getEndCapType());
    }

    NodeDescriptor desc = GetDescriptor();

    int numPatchPoints = 0;
    switch (desc.GetType()) {

        using enum NodeType;

        case NODE_REGULAR: 
            numPatchPoints = regularPatchSize; break;

        case NODE_END:
            numPatchPoints = irregularPatchSize; break;

        case NODE_RECURSIVE: 
            numPatchPoints = desc.HasEndcap() ? irregularPatchSize : 0; break;

        // note : terminal nodes hold 25 patch points, but the public-facing API
        // client code selects which regular patch to use and automatically
        // reduces the 25-set to 16 indices
        case NODE_TERMINAL:
            if (desc.GetDepth() >= maxLevel) {
                if (desc.HasEndcap())
                    numPatchPoints = irregularPatchSize;
            } else if ((int)desc.GetEvIndex() != quadrant)
                numPatchPoints = regularPatchSize;
            break;
    }
    return numPatchPoints;
}

Index
SubdivisionPlan::Node::GetPatchPoint(int pointIndex, int quadrant, unsigned short maxLevel) const {

    uint32_t const* tree = plan->GetPatchTreeData().data();

    Index offset = tree[patchPointsOffset()];

    if (offset == INDEX_INVALID)
        return INDEX_INVALID;

    NodeDescriptor desc = GetDescriptor();
    switch (desc.GetType()) {
        using enum NodeType;
        case NODE_REGULAR:
        case NODE_END:
            offset += pointIndex; break;
        case NODE_TERMINAL: {
            if (desc.GetDepth() >= maxLevel) {
                // if we hit the dynamic isolation cap, return the end-cap patch point
                // indices that are stored after the 5x5 grid of regular patch
                // point indices.
                if (desc.HasEndcap())
                    offset += pointIndex + catmarkTerminalPatchSize();
            } else {
                // note : Z winding order is 0, 1, 3, 2 !
                static int permuteTerminal[4][16] = {
                    {0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18},
                    {1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19},
                    {5, 6, 7, 8, 10, 11, 12, 13, 15, 16, 17, 18, 20, 21, 22, 23},
                    {6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 24},
                };
                assert(quadrant != (int)desc.GetEvIndex() && quadrant >= 0 && quadrant < 4);
                offset += permuteTerminal[quadrant][pointIndex];
            }
        } break;
        case NODE_RECURSIVE:
            offset = (desc.GetDepth() >= maxLevel) && desc.HasEndcap() ?
                offset + pointIndex : INDEX_INVALID;
    }
    assert(offset != INDEX_INVALID);
    return plan->_patchPoints[offset];
}

Vtr::ConstArray<float>
SubdivisionPlan::Node::GetPatchPointWeights(int pointIndex) const {
    
    std::vector<float> const& weights = plan->GetStencilMatrix();
    
    if (weights.empty())
        return Vtr::ConstArray<float>(nullptr, 0);

    int numControlPoints = plan->GetNumControlPoints();

    return Vtr::ConstArray<float>(weights.data() + (pointIndex * numControlPoints), numControlPoints);
}

//
// Subdivision Plan
//

SubdivisionPlan::SubdivisionPlan(int treeSize,
    int numControlPoints, int numPatchPoints, int regFaceSize, int faceSize, int subfaceIndex) { 

    static_assert(sizeof(TreeDescriptor) ==
        NodeBase::rootNodeOffset() * sizeof(decltype(SubdivisionPlan::_tree)::value_type));

    static_assert(NodeBase::maxIsolationLevel() == kMaxIsolationLevel);

    _tree.resize(treeSize);

    TreeDescriptor* tdesc = reinterpret_cast<TreeDescriptor*>(_tree.data());

    tdesc->Set(faceSize == regFaceSize, faceSize, subfaceIndex, numControlPoints);

    _patchPoints.resize(numPatchPoints);
}

SubdivisionPlan::~SubdivisionPlan() {
    // Plans own the lifespan of their neighborhoods : we need to
    // delete them
    clearNeighborhoods();
}

SubdivisionPlan::Node
SubdivisionPlan::GetNode(float s, float t, unsigned char& quadrant, unsigned short level) const {
    switch (getSchemeType()) {
        case Sdc::SchemeType::SCHEME_CATMARK: 
            return Node(this, NodeBase::traverse<NodeBase::CatmarkWinding>(_tree.data(), s, t, quadrant, level).treeOffset);
        case Sdc::SchemeType::SCHEME_LOOP:
            return Node(this, NodeBase::traverse<NodeBase::LoopWinding>(_tree.data(), s, t, quadrant, level).treeOffset);
        default:
            break;
    }
    assert(0);
    return {};
}

// translators for Tmr => Far types
inline Far::PatchDescriptor::Type 
regularBasisType(Sdc::SchemeType scheme) {

    using enum Far::PatchDescriptor::Type;       
    switch (scheme) {
        case Sdc::SCHEME_CATMARK: return REGULAR;
        case Sdc::SCHEME_LOOP: return LOOP;
        default:
            break;
    }
    return NON_PATCH;
}

inline Far::PatchDescriptor::Type
irregularBasisType(Sdc::SchemeType scheme, EndCapType endcap) {

    using enum Far::PatchDescriptor::Type;
    if (scheme == Sdc::SCHEME_CATMARK) {
        switch (endcap) {
            case EndCapType::ENDCAP_BILINEAR_BASIS: return QUADS;
            case EndCapType::ENDCAP_BSPLINE_BASIS: return REGULAR;
            case EndCapType::ENDCAP_GREGORY_BASIS: return GREGORY_BASIS;
            default:
                break;
        }
    } else if (scheme == Sdc::SCHEME_LOOP) {
        switch (endcap) {
            case EndCapType::ENDCAP_BILINEAR_BASIS: return TRIANGLES;
            case EndCapType::ENDCAP_BSPLINE_BASIS: return LOOP;
            case EndCapType::ENDCAP_GREGORY_BASIS: return GREGORY_TRIANGLE;
            default:
                break;
        }
    }
    return Far::PatchDescriptor::NON_PATCH;
}

template <typename REAL> SubdivisionPlan::Node
SubdivisionPlan::evaluateBasis(REAL s, REAL t, REAL wP[],
    REAL wDs[], REAL wDt[], REAL wDss[], REAL wDst[], REAL wDtt[], unsigned char* subpatch, short level) const {

    Sdc::SchemeType const scheme = getSchemeType();
    Far::PatchDescriptor::Type const regularBasis = regularBasisType(scheme);
    Far::PatchDescriptor::Type const irregularBasis = irregularBasisType(scheme, getEndCapType());

    bool nonQuad = !GetTreeDescriptor().IsRegularFace();

    unsigned char quadrant = 0;
    Node n = GetNode(float(s), float(t), quadrant, level);

    NodeDescriptor desc = n.GetDescriptor();

    NodeType nodeType = desc.GetType();
    int depth = desc.GetDepth();

    using enum NodeType;
    bool dynamicIsolation = (nodeType == NODE_RECURSIVE || nodeType == NODE_TERMINAL) && (depth >= level) && desc.HasEndcap();

    unsigned short u = desc.GetU();
    unsigned short v = desc.GetV();

    Far::PatchParam param;

    if (dynamicIsolation) {

        if (nodeType == NODE_TERMINAL) {
            u = u >> 1;
            v = v >> 1;
        }
        param.Set(INDEX_INVALID, u, v, depth, nonQuad, 0, 0, true);
        Far::internal::EvaluatePatchBasis<REAL>(irregularBasis, param, s, t, wP, wDs, wDt, wDss, wDst, wDtt);

    } else {

        switch (nodeType) {

            case NODE_REGULAR : {
                param.Set(INDEX_INVALID, u, v, depth, nonQuad, desc.GetBoundaryMask(), 0, true);
                REAL sharpness = computeSingleCreaseSharpness<REAL>(n, desc, depth, level);   
                Far::internal::EvaluatePatchBasis<REAL>(regularBasis, param, s, t, wP, wDs, wDt, wDss, wDst, wDtt, sharpness);
            } break;

            case NODE_END : {
                param.Set(INDEX_INVALID, u, v, depth, nonQuad, desc.GetBoundaryMask(), 0, true);
                Far::internal::EvaluatePatchBasis<REAL>(irregularBasis, param, s, t, wP, wDs, wDt, wDss, wDst, wDtt);
            } break;

            case NODE_TERMINAL : {
                assert(scheme == Sdc::SCHEME_CATMARK);
                switch (quadrant) {
                    case 0 :                 break;
                    case 1 : { u+=1;       } break;
                    case 2 : {       v+=1; } break; // Z order '^' bitwise winding
                    case 3 : { u+=1; v+=1; } break;
                }
                param.Set(INDEX_INVALID, u, v, depth+1, nonQuad, 0, 0, true);
                Far::internal::EvaluatePatchBasis<REAL>(regularBasis, param, s, t, wP, wDs, wDt, wDss, wDst, wDtt);
            } break;

            default:
                assert(0);
        }
    }
    if (subpatch)
        *subpatch = quadrant;

    return n;
}

template SubdivisionPlan::Node
SubdivisionPlan::evaluateBasis<float>(float s, float t, float wP[],
    float wDs[], float wDt[], float wDss[], float wDst[], float wDtt[], unsigned char* subpatch, short level) const;

template SubdivisionPlan::Node
SubdivisionPlan::evaluateBasis<double>(double s, double t, double wP[],
    double wDs[], double wDt[], double wDss[], double wDst[], double wDtt[], unsigned char* subpatch, short level) const;

inline Sdc::SchemeType 
SubdivisionPlan::getSchemeType() const {
    return GetTopologyMap()->GetTraits().getSchemeType(); }

inline EndCapType
SubdivisionPlan::getEndCapType() const {
    return GetTopologyMap()->GetTraits().getEndCapType();
}

//
// Neighborhoods
//

void
SubdivisionPlan::reserveNeighborhoods(int count) {
    assert(_neighborhoods.GetSize() == 0);
    _neighborhoods.Reserve(count);
#if defined(DEBUG) || defined(_DEBUG)
    std::memset(&_neighborhoods[0], 0, count * sizeof(decltype(_neighborhoods[0])));
#endif
}

void
SubdivisionPlan::addNeighborhood(std::unique_ptr<uint8_t const[]> neighborhoodData, int /*startEdge*/) {

    Neighborhood const& neighborhood = *reinterpret_cast<Neighborhood const*>(neighborhoodData.get());

    for (int i = 0; i < GetNumNeighborhoods(); ++i)
        if (neighborhood.IsEquivalent(GetNeighborhood(i)))
            return;

    Index index = _neighborhoods.GetSize();

    // SetSize is a no-op since since 'reserve' above used the face valence to make sure
    // we have enough space (and revert to dynamic allocation for 'big valence' plans.
    // note : we are abusing the StackBuffer API since it does not expose 'push_back'
    // or 'resize' functions - we have to make sure that we always 'Reserve' enough space,
    // resort to storing a sparse vector with nullptrs or revert back to a more conventional
    // std::vector
    _neighborhoods.SetSize(index + 1);
    _neighborhoods[index] = neighborhoodData.release();
}

void
SubdivisionPlan::clearNeighborhoods() {
    for (int i = 0; i < (int)_neighborhoods.GetSize(); ++i) {
        delete[] _neighborhoods[i];
#if defined(DEBUG) || defined(_DEBUG)
        _neighborhoods[i] = nullptr;
#endif
    }
}

bool 
SubdivisionPlan::IsTopologicalMatch(Neighborhood const& neighborhood, int& startingEdge) const {
    for (int i = 0; i < GetNumNeighborhoods(); ++i) {
        if (neighborhood.IsEquivalent(GetNeighborhood(i))) {
            startingEdge = GetNeighborhood(i).GetStartingEdge();
            return true;
        }
    }
    return false;
}

//
//
//

size_t
SubdivisionPlan::GetByteSize(bool includeNeighborhoods) const {

    auto getByteSize = []<typename T>(std::vector<T> const& vec) -> size_t {
        return vec.size() * sizeof(T);
    };

    size_t size = sizeof(SubdivisionPlan);
    size += getByteSize(GetPatchTreeData());
    size += getByteSize(GetPatchPoints());
    size += getByteSize(GetStencilMatrix());

    if (includeNeighborhoods) {
        for (int nIndex = 0; nIndex < GetNumNeighborhoods(); ++nIndex) {
            Neighborhood const& n = GetNeighborhood(nIndex);
            if constexpr (!Tmr::NeighborhoodBuilder::kDebugNeighborhoods)
                assert(!n.HasControlPoints());
            size += n.GetByteSize();
        }
    }
    return size;
}

//
// Debug functions
//

void SubdivisionPlan::Node::PrintStencil() const {

    std::vector<uint32_t> const& tree = plan->GetPatchTreeData();
    std::vector<Index> const& patchPoints = plan->GetPatchPoints();
    std::vector<float> const& stencilMatrix = plan->GetStencilMatrix();

    int numControlPoints = plan->GetNumControlPoints();

    NodeDescriptor descr = GetDescriptor();

    auto print = [&](int offset, int numPatchPoints, char const* type) {

        printf("node (type='%s' fsup=%d depth=%d u=%d v=%d) = {\n", 
            type, offset, descr.GetDepth(), descr.GetU(), descr.GetV());
        for (int i = 0; i < numPatchPoints; ++i) {
            if (Index index = patchPoints[offset + i]; index < numControlPoints)
                printf("    control vertex : %3d\n", index);
            else {
                index -= numControlPoints;
                printf("    matrix row     : %3d (%4d) = {", index, index * numControlPoints);
                for (int i = 0; i < numControlPoints; ++i)
                    printf("% .06f ", stencilMatrix[index * numControlPoints + i]);
                printf("}\n");
            }
        }
        printf("}\n");
    };

    constexpr int regularSize = catmarkRegularPatchSize();
    constexpr int terminalSize = catmarkTerminalPatchSize();
    int irregularSize = catmarkIrregularPatchSize(plan->getEndCapType());

    int patchPointOffset = tree[patchPointsOffset()];

    switch (descr.GetType()) {
        using enum NodeType;      
        case NODE_REGULAR: 
            print(patchPointOffset, regularSize, "REG"); 
            break;
        case NODE_END:
            print(patchPointOffset, irregularSize, "END");
            break;
        case NODE_RECURSIVE:
            if (descr.HasEndcap())
                print(patchPointOffset, irregularSize, "REC");
            break;
        case NODE_TERMINAL:
            print(patchPointOffset, terminalSize, "TERM - base");
            if (descr.HasEndcap())
                print(patchPointOffset + terminalSize, irregularSize, "TERM - end");
            break;
    }
}
void
printNodeIndices(FILE* fout, ConstIndexArray indices) {

    if (indices.empty())
        return;

    int stride = (indices.size() == 16) ? 4 : 5;
    fprintf(fout, "\\n");
    for (int i = 1; i <= indices.size(); ++i) {
        if (((i-1) % stride) != 0)
            fprintf(fout, " ");
        fprintf(fout, "%d ", indices[i-1]);
        if ((i % stride) == 0)
            fprintf(fout, "\\n");
    }
}

void
SubdivisionPlan::WriteTreeDigraph(FILE* fout, int planIndex, bool showIndices, bool isSubgraph) const {

    constexpr int const regularPatchSize = NodeBase::catmarkRegularPatchSize();
    constexpr int const terminalPatchSize = NodeBase::catmarkTerminalPatchSize();
    int const irregularPatchSize = NodeBase::catmarkIrregularPatchSize(getEndCapType());
    
    auto hashNodeId = [&planIndex](Node const& n) -> uint64_t {
        return (uint64_t(planIndex) << 32) + n.treeOffset;
    };

    if (isSubgraph) {
        fprintf(fout, "subgraph cluster_%d {\n", planIndex);
        fprintf(fout, "  label = \"Plan %d\"; style=filled; color=lightgrey;\n", planIndex);
    } else {
        fprintf(fout, "digraph {\n");
    }

    int count = 0;
    for (Node node = GetRootNode(); node.IsValid(); ++node, ++count) {

        NodeDescriptor const& desc = node.GetDescriptor();

        uint64_t nodeID = hashNodeId(node);

        int treeOffset = node.treeOffset;
        int patchPointsOffset = _tree[node.patchPointsOffset()];

        switch (desc.GetType()) {

            using enum NodeType;

            case NODE_REGULAR: {
                fprintf(fout, "  %zu [label=\"%d REG (%d)\\ntofs=%d fsup=%d",
                    nodeID, count, desc.GetDepth(), treeOffset, patchPointsOffset);
                if (showIndices)
                    printNodeIndices(fout, { &_patchPoints[patchPointsOffset], regularPatchSize});
                if (desc.HasSharpness()) {
                    fprintf(fout, "\\n\\nsharp=%f", node.GetSharpness());
                }
                fprintf(fout, "\", shape=box, style=filled, color=%s]\n", desc.HasSharpness() ? "darksalmon" : "white");
            } break;

            case NODE_END: {
                fprintf(fout, "  %zu [label=\"%d END (%d)\\ntofs=%d fsup=%d",
                    nodeID, count, desc.GetDepth(), treeOffset, patchPointsOffset);
                if (showIndices)
                    printNodeIndices(fout, { &_patchPoints[patchPointsOffset], irregularPatchSize });
                fprintf(fout, "\", shape=box, style=filled, color=darkorange]\n");
            } break;

            case NODE_RECURSIVE: {
                fprintf(fout, "  %zu [label=\"%d REC (%d)\\ntofs=%d fsup=%d",
                    nodeID, count, desc.GetDepth(), treeOffset, patchPointsOffset);
                if (showIndices && desc.HasEndcap())
                    printNodeIndices(fout, { &_patchPoints[patchPointsOffset], irregularPatchSize });
                fprintf(fout, "\", shape=square, style=filled, color=dodgerblue]\n");
                for (int i = 0; i < 4; ++i)
                    fprintf(fout, "  %zu -> %zu [label=\"%d\"]\n", nodeID, hashNodeId(node.GetChild(i)), i);
            } break;

            case NODE_TERMINAL: {
                fprintf(fout, "  %zu [label=\"%d TRM (%d)\\ntofs=%d fsup=%d",
                    nodeID, count, desc.GetDepth(), treeOffset, patchPointsOffset);
                if (showIndices) {
                    printNodeIndices(fout, { &_patchPoints[patchPointsOffset], terminalPatchSize} );
                    printNodeIndices(fout, { &_patchPoints[patchPointsOffset + terminalPatchSize], irregularPatchSize });
                }
                fprintf(fout, "\", shape=box, style=filled, color=grey]\n");
                fprintf(fout, "  %zu -> %zu\n", nodeID, hashNodeId(node.GetChild()));
            } break;
            default:
                assert(0);
        }
    }
    fprintf(fout, "}\n");
}


} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
