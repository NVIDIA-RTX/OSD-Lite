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

#ifndef OPENSUBDIV3_TMR_SUBDIVISION_PLAN_BUILDER_H
#define OPENSUBDIV3_TMR_SUBDIVISION_PLAN_BUILDER_H

#include "../version.h"

#include "../tmr/subdivisionPlan.h"
#include "../tmr/types.h"

#include "../far/sparseMatrix.h"

#include <array>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Far {
    class PatchBuilder;
    class TopologyRefiner;
}

namespace Tmr {

class Neighborhood;

//
// A specialized builder for subdivision plan hierarchies
//

class SubdivisionPlanBuilder {

public:

    SubdivisionPlanBuilder();
    ~SubdivisionPlanBuilder();

    struct Options {

        // Type of patch used to approximate limit surface areas covered by Extraordinay nodes
        EndCapType endCapType = EndCapType::ENDCAP_BSPLINE_BASIS;

        // Limits adaptive isolation to a given level
        // note : given that infinite sharpness is set to 10.f isolation
        // levels are clamped to a maximum value of kMaxLevel = 10
        // See Far::TopologyRefiner::AdaptiveOptions for more details
        uint8_t isolationLevel = 6;
        uint8_t isolationLevelSecondary = 2;

        // see Far::AdaptiveOptions::useSingleCreasePatch
        bool useSingleCreasePatch = true;

        // see Far::AdaptiveOptions::useInfSharpPatch
        bool useInfSharpPatch = true;

        // Enable Terminal nodes optimization: terminal nodes replace 3 adjancent
        // Regular nodes surrounding an Extraordinary node with a Terminal node.
        // This optimization reduces the number of patch points from 16 x 3 = 48 to 25
        bool useTerminalNode = true; 

        // Enable Dynamic isolation of extraordinary vertices: will insert an Extraordinary
        // node at every level of the tree instead of just the deepest one. This allows
        // traversal to stop at an arbitrary level and have a fall-back approximation to 
        // cover the limit surface areas covered by Recursive nodes.
        bool useDynamicIsolation = true; 

        // By default, the stencil matrix positions the stencil rows for irregular patch
        // points at the end of the matrix. This setting forces a re-ordering of the rows
        // so that all the stencils for a given level of isolation are contiguous in the
        // matrix (requires useDynamicIsolation = true)
        bool orderStencilMatrixByLevel = false;

        // legacy behaviors

        // Generate sharp regular patches at smooth corners (legacy)
        bool generateLegacySharpCornerPatches = false; 

        bool operator == (Options const& other) const;
    };

    // Returns a subdivision plan for the given topological neighborhood.
    // note: the builder maintains state between build calls ; sub-faces of an irregular
    // face should be built in sequence to take advantage of the builder's cached state
    // (subfaceIndex = 0 invalidates the cache, subfaceIndex > 0 re-uses the cache)
    std::unique_ptr<SubdivisionPlan> Create(Sdc::SchemeType scheme, Sdc::Options schemeOptions, 
        Options const& options, Neighborhood const& neighborhood, LocalIndex subfaceIndex = 0);

private:

    // Face Tags

    struct FaceTags {

        uint32_t hasPatch        : 1,
                 isRegular       : 1,
                 boundaryMask    : 5,
                 boundaryIndex   : 2,
                 boundaryCount   : 3,
                 hasBoundaryEdge : 3,
                 isSingleCrease  : 1;

        float sharpness; // single crease edge sharpness

        void clear() { std::memset(this, 0, sizeof(FaceTags)); }
    };

    void setFaceTags(FaceTags& tags, int levelIndex, LocalIndex faceIndex) const;

    // Proto Nodes

    struct ProtoNode {

        LocalIndex faceIndex; // index of face in Vtr::level

        FaceTags faceTags;

        uint32_t active            : 1,  // proto-node will not generate a node if false
                 nodeType          : 2,  // NodeType of generated node
                 numChildren       : 3,  // number of children linked to the node
                 hasIrregularPatch : 1,  // node will contain an irregular patch
                 levelIndex        : 4,  // index of Vtr::level
                 evIndex           : 2;  // index of ev for terminal nodes

        uint32_t treeOffset; // offset to corresponding Node in plan's Nodes array
        uint32_t firstPatchPoint; // index of first patch point for the node

        uint32_t children[4]; // indices of children nodes in the node-store

        NodeType getType() const { return NodeType(nodeType); }
        void setType(NodeType type) { nodeType = uint8_t(type); }
    };

    // returns a child of the node (undetermined if child does not exist !)
    ProtoNode& getProtoNodeChild(ProtoNode const& node, LocalIndex childIndex);
    ProtoNode const& getProtoNodeChild(ProtoNode const& node, LocalIndex childIndex) const;

    // encodes proto-nodes into SubdivisionPlan::Nodes
    void encodeRegularNode(ProtoNode const& pn, uint32_t* treePtr, Index* patchPoints);
    void encodeEndCapNode(ProtoNode const& pn, uint32_t* treePtr, Index* patchPoints);
    void encodeTerminalNode(ProtoNode const& pn, uint32_t* treePtr, Index* patchPoints);
    void encodeRecursiveNode(ProtoNode const& pn, uint32_t* treePtr, Index* patchPoints);

    void initializeProtoNodes(LocalIndex subfaceIndex);
    void identifyTerminalNodes();
    unsigned int finalizeProtoNodes();

private:    

    // misc. helpers

    bool computeLoopSubPatchDomain(int levelIndex, Index faceIndex, short& s, short& t) const;
    bool computeCatmarkSubPatchDomain(int levelIndex, Index faceIndex, short& s, short& t) const;
    bool computeSubPatchDomain(int levelIndex, Index faceIndex, short& s, short& t) const;

    void gatherIrregularPatchPoints(int levelIndex, Index* patchPoints);

    void getIrregularPatchConversion(ProtoNode const& pn);
    void appendConversionStencilsToMatrix(std::vector<float>& stencilMatrix, int stencilBaseIndex);

    void initializePatches(SubdivisionPlan& plan);
    void initializeStencilMatrix(SubdivisionPlan& plan);
    void reorderStencilMatrix(SubdivisionPlan& plan);

private:

    static constexpr int const kNumLevels = kMaxIsolationLevel + 1;

    std::vector<ProtoNode> _protoNodeStore[kNumLevels];

    std::unique_ptr<Far::TopologyRefiner> _faceRefiner;
    std::unique_ptr<Far::PatchBuilder> _patchBuilder;

    struct LevelOffsets {
        Index regularIndexBase = INDEX_INVALID;
        Index irregularIndexBase = INDEX_INVALID;
        unsigned int irregularIndexOffset = 0;
        unsigned int levelOffset = 0;
        unsigned int numIrregularPatches = 0;
    } _levelOffsets[kNumLevels];

    struct IrregularPatchConversion {
        Far::SparseMatrix<float> matrix;
        std::vector<Index> sourcePoints;
    } _conversion;

    bool _useDynamicIsolation = false;
    bool _reorderStencilMatrix = false;

    int _regularFaceSize = 0;
    int _regularPatchSize = 0;
    int _irregularPatchSize = 0;
    int _terminalPatchSize = NodeBase::catmarkTerminalPatchSize();

    int _faceSize = 0;
    int _numControlPoints = 0;
    int _numRefinedPoints = 0;
    int _numPatchPointsTotal = 0;
    int _numIrregularPatchesTotal = 0;
};

} // end namespace Tmr
} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_SUBDIVISION_PLAN_BUILDER_H */

