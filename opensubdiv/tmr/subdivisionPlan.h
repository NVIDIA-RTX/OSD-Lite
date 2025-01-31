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


#ifndef OPENSUBDIV3_TMR_SUBDIVISION_PLAN_H
#define OPENSUBDIV3_TMR_SUBDIVISION_PLAN_H

#include "../version.h"

#include "../tmr/nodeDescriptor.h"
#include "../tmr/nodeBase.h"
#include "../tmr/treeDescriptor.h"
#include "../tmr/types.h"
#include "../vtr/stackBuffer.h"
#include "../sdc/scheme.h"

#include <array>
#include <memory>
#include <vector>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

class SubdivisionPlanBuilder;
class TopologyMap;
class Neighborhood;

///
///  \brief Stores a subdivision plan
///
/// The subdivision plan is a data structure that represents the feature-adaptive
/// isolation hierarchy of sub-patches for a face, down to some fixed, maximum depth.
///
/// Specifically, it comprises:
///     - an optimized quad-tree representing the hierarchy of sub-patches
///       that describe the limit surface of this face.
///     - a stencil matrix for evaluation of the patch-points supporting the
///       sub-patches.
///
/// * The plan for a face depends only on the configuration of elements that
///   can exert an influence on the limit surface within its local domain. This
///   includes the topology of the face and its 1-ring, sharpness tags for
///   incident edges and vertices, and boundary rules. (see Tmr::NeighborhoodBuilder)
///
/// * The quad-tree represents the hierarchy of feature adaptive sub-patches. Its
///   nodes are encoded in a flat array of integers, with offsets pointing to children
///   nodes. The encoding is designed to be device-friendly. (see Tmr::NodeBase)
///
///  For more details, see:
///  "Efficient GPU Rendering of SUbdivision Surfaces using Adaptive Quadtrees"
///    W. Brainer, T. Foley, M. Mkraemer, H. Moreton, M. Niessner - Siggraph 2016
///
///  http://www.graphics.stanford.edu/~niessner/papers/2016/4subdiv/brainerd2016efficient.pdf
///
class SubdivisionPlan {

public:

    ~SubdivisionPlan();

    /// \brief Returns the map this plan belongs to
    TopologyMap const* GetTopologyMap() const { return _topologyMap; }

    TreeDescriptor const& GetTreeDescriptor() const { return *reinterpret_cast<TreeDescriptor const*>(_tree.data()); }

    /// \brief Returns true if this topology plan belongs to a regular face
    bool IsRegularFace() const { return GetTreeDescriptor().IsRegularFace(); }

    /// \briaf Returns the valence of the face
    int GetFaceSize() const { return GetTreeDescriptor().GetFaceSize(); }

    /// \briaf Returns the index of the sub-patch if the face is not a quad
    /// and INDEX_INVALID if the face is a quad
    int GetSubfaceIndex() const { return GetTreeDescriptor().GetSubfaceIndex(); }

    /// \brief Returns the number of control vertices in the plan's neighborhood
    unsigned short GetNumControlPoints() const { return GetTreeDescriptor().GetNumControlPoints(); }

    /// \brief Returns the number of patch points that need to be evaluated in
    /// order to support all the sub-patches up to the given isolation level
    /// note: 'level' equires enabling both the dynamic isolation and stencil
    /// matrix reordering features in SurfaceTableFactory::Options
    unsigned int GetNumPatchPoints(int level = kMaxIsolationLevel) const;

    /// \brief Returns true if the topology that the plan was built from is equivalent
    /// to that of the input neighborhood (and if so, what the starting edge is)
    bool IsTopologicalMatch(Neighborhood const& neighborhood, int& startingEdge) const;

    /// \brief Returns the number of (rotated) neighborhoods
    int GetNumNeighborhoods() const { return (int)_neighborhoods.GetSize(); }

    /// \brief Returns the neighborhood at 'index'
    Neighborhood const& GetNeighborhood(int index) const;

public:

    struct Node;

    /// \brief Evaluate basis functions for position and first derivatives at a
    /// given (s,t) parametric location of a patch.
    ///
    /// @param s             Patch coordinate (in coarse face normalized space)
    /// @param t             Patch coordinate (in coarse face normalized space)
    ///
    /// @param wP            Weights (evaluated basis functions) for the position
    /// @param wDs           Weights (evaluated basis functions) for derivative wrt s
    /// @param wDt           Weights (evaluated basis functions) for derivative wrt t
    ///
    /// @param[out] quadrant Domain quadrant containing (s,t) location
    ///                      (required in order to obtain the correct patch points for
    ///                      a terminal node)
    ///
    /// @param[int] level    Dynamic level of isolation (if possible)
    ///
    /// @return              The leaf node pointing to the sub-patch evaluated

    template <typename REAL>
    Node evaluateBasis(REAL s, REAL t, REAL wP[],
        REAL wDs[], REAL wDt[], REAL wDss[], REAL wDst[], REAL wDtt[],
        unsigned char* quadrant, short level = kMaxIsolationLevel) const;

    Node EvaluateBasis(float s, float t, float wP[],
        float wDs[], float wDt[], float wDss[], float wDst[], float wDtt[],
            unsigned char* quadrant, short level = kMaxIsolationLevel) const;

    Node EvaluateBasis(double s, double t, double wP[],
        double wDs[], double wDt[], double wDss[], double wDst[], double wDtt[],
        unsigned char* quadrant, short level = kMaxIsolationLevel) const;

    ///
    /// Patch points data interpolation
    /// 
    /// \note EvaluatePatchPoints templates both the source and destination
    ///      data buffer classes. Client-code is expected to provide interfaces
    ///      that implement the functions specific to its primitive variable
    ///      data layout. Template APIs must implement the following:
    ///      <br><br> \code{.cpp}
    ///
    ///      class MySource {
    ///          MySource & operator[](int index);
    ///      };
    ///
    ///      class MyDestination {
    ///          void Clear();
    ///          void AddWithWeight(MySource const & value, float weight);
    ///          void AddWithWeight(MyDestination const & value, float weight);
    ///      };
    ///
    ///      \endcode
    ///    
    ///      It is possible to implement a single interface only and use it as
    ///      both source and destination.
    /// 
    template <typename T, typename U> int EvaluatePatchPoints(T const* controlPoints, 
        ConstIndexArray controlPointIndices, U* patchPoints, int level = kMaxIsolationLevel) const;

public:


    /// \brief Node in the tree of sub-patches
    ///
    /// note : The burden is on the client to check whether a particular accessor
    ///        method can be applied on a given node. If the node is of the wrong
    ///        type, behavior will be "undefined".

    struct Node : public NodeBase {

        Node(SubdivisionPlan const* plan = nullptr, int treeOffset = -1) : plan(plan) { this->treeOffset = treeOffset; }

        /// \brief Returns false if un-initialized
        bool IsValid() const { return plan && treeOffset >= 0; }

        /// \brief Returns a pointer to the plan that owns this node
        SubdivisionPlan const* GetSubdivisionPlan() const { return plan; }

        /// \brief Returns the node's descriptor
        inline NodeDescriptor GetDescriptor() const;

        /// \brief Returns the crease sharpness of the node
        /// The node is expected to be 'Regular' and flagged as a 'single-crease' patch
        inline float GetSharpness() const;

        /// \brief Returns the number of children of the node
        inline int GetNumChildren() const;

        /// \brief Returns the child of the node
        inline Node GetChild(int childIndex=0) const;

        /// \brief Returns the number of patch points points supporting the sub-patch
        int GetPatchSize(int quadrant, unsigned short level = kMaxIsolationLevel + 1) const;

        /// \brief Returns the index of the requested patch point in the sub-patch
        /// described by the node (pointIndex is in range [0, GetPatchSize()])
        /// 
        /// @oaram pointIndex  Index of the patch point 
        /// @param quadrant    One of the 3 quadrants of a Terminal patch that are
        ///                    not the quadrant containing the extraordinary vertex
        /// @param level       Desired dynamic level of isolation
        Index GetPatchPoint(int pointIndex, int quadrant = 0, unsigned short level = kMaxIsolationLevel + 1) const;

        /// \brief Returns the row corresponding to the patch point in the
        /// dense stencil matrix
        ConstFloatArray GetPatchPointWeights(int pointIndex) const;

        /// \brief Returns the next node in the tree (serial traversal). Check 
        //  IsValid() for iteration past end-of-tree
        Node operator ++ ();

        /// \brief Returns true if the nodes are identical
        bool operator == (Node const & other) const;

        void PrintStencil() const;

        SubdivisionPlan const* plan = nullptr;
    };

    /// \brief Returns a pointer to the root node of the sub-patches tree
    Node GetRootNode() const { return Node(this, 2); }

    /// \brief Returns the node corresponding to the sub-patch at the given (s,t) location
    Node GetNode(float s, float t, 
        unsigned char& quadrant, unsigned short level = kMaxIsolationLevel + 1) const;

public:

    // Raw data accessor
    std::vector<uint32_t> const& GetPatchTreeData() const { return _tree; }
    std::vector<Index> const& GetPatchPoints() const { return _patchPoints; }
    std::vector<float> const& GetStencilMatrix() const { return _stencilMatrix; }

    /// \brief Returns the number of stencils (ie. the number of stencil matrix rows)
    size_t GetNumStencils() const;

    /// \brief Returns the total amount of host memory used by the plan
    size_t GetByteSize(bool includeNeighborhoods = false) const;

    /// \brief Writes a GraphViz 'dot' digraph of the tree
    void WriteTreeDigraph(FILE* fout,
        int planIndex=0, bool showIndices=true, bool isSubgraph=false) const;

private:

    friend class SubdivisionPlanBuilder;
    friend class SurfaceTableFactory;
    friend class TopologyMap;

    SubdivisionPlan(int treeSize, int numControlPoints, int numPatchPoints,
        int regFaceSize, int faceSize, int subfaceIndex = INDEX_INVALID);


    Sdc::SchemeType getSchemeType() const;
    EndCapType getEndCapType() const;

    void reserveNeighborhoods(int count);
    void addNeighborhood(std::unique_ptr<uint8_t const[]> neighborhoodData, int startEdge);
    void clearNeighborhoods();

private:

    TopologyMap const* _topologyMap = nullptr;

    // Topology description of the 1-ring control vertices : we assume that the
    // hash function can potentially generate collisions, so we keep the 
    // compacted 1-ring neighborhood descriptions for exact comparison.
    Vtr::internal::StackBuffer<uint8_t const*, 4, true> _neighborhoods;

    // The sub-patch "tree" is stored as a linear buffer of integers for
    // efficient look-up & traversal on a GPU. Use the Node class to traverse
    // the tree and access each node's data.
    std::vector<uint32_t> _tree;
    static_assert(sizeof(NodeDescriptor) == sizeof(decltype(_tree)::value_type));

    // Stencil indices into the _weights
    std::vector<Index> _patchPoints;

    // Stencil matrix for computing patch points from control points:
    // - columns contain 1 scalar weight per control point of the 1-ring 
    // - rows contain a stencil of weights for each patch point
    std::vector<float> _stencilMatrix;
};

//
// Inline implementation
//

inline SubdivisionPlan::Node
SubdivisionPlan::EvaluateBasis(float s, float t, float wP[], float wDs[], float wDt[],
    float wDss[], float wDst[], float wDtt[], unsigned char* quadrant, short level) const {
    return evaluateBasis<float>(s, t, wP, wDs, wDt, wDss, wDst, wDtt, quadrant, level);
}

inline SubdivisionPlan::Node
SubdivisionPlan::EvaluateBasis(double s, double t, double wP[], double wDs[], double wDt[],
    double wDss[], double wDst[], double wDtt[], unsigned char* quadrant, short level) const {
    return evaluateBasis<double>(s, t, wP, wDs, wDt, wDss, wDst, wDtt, quadrant, level);
}

template <typename T, typename U> int
SubdivisionPlan::EvaluatePatchPoints(
    T const* controlPoints, ConstIndexArray controlPointIndices, U* patchPoints, int level) const {

    assert(level <= kMaxIsolationLevel);

    std::vector<float> const& stencilMatrix = GetStencilMatrix();

    unsigned int numPatchPoints = GetNumPatchPoints(level);
    int numControlPoints = controlPointIndices.size();

    assert(numControlPoints == GetNumControlPoints());

    for (unsigned int point = 0; point < numPatchPoints; ++point) {
        
        U& patchPoint = patchPoints[point];

        patchPoint.Clear();

        for (short i = 0; i < numControlPoints; ++i) {

            float w = stencilMatrix[point * numControlPoints + i];
        
            T const& controlPoint = controlPoints[controlPointIndices[i]];

            patchPoint.AddWithWeight(controlPoint, w);
        }
    }
    return numPatchPoints;
}

inline Neighborhood const& 
SubdivisionPlan::GetNeighborhood(int index) const {
    assert(index < (int)_neighborhoods.GetSize());
    return *reinterpret_cast<Neighborhood const*>(_neighborhoods[index]);
}

inline unsigned int 
SubdivisionPlan::GetNumPatchPoints(int level) const {
    return GetTreeDescriptor().GetNumPatchPoints(level);
}

inline size_t SubdivisionPlan::GetNumStencils() const {
    return _stencilMatrix.size() / GetNumControlPoints();
}

inline NodeDescriptor 
SubdivisionPlan::Node::GetDescriptor() const {
    return NodeDescriptor{static_cast<uint32_t>(plan->GetPatchTreeData()[descriptorOffset()])};
}

inline int 
SubdivisionPlan::Node::GetNumChildren() const {
    return getNumChilren(GetDescriptor().GetType());
}

inline SubdivisionPlan::Node
SubdivisionPlan::Node::GetChild(int childIndex) const {
    
    assert(childIndex < 4);

    uint32_t const* tree = plan->GetPatchTreeData().data();
    
    return Node(plan, tree[childOffset(childIndex)]);
}

inline float
SubdivisionPlan::Node::GetSharpness() const {

    assert(GetDescriptor().GetType()==NodeType::NODE_REGULAR
        && GetDescriptor().HasSharpness());

    uint32_t const* tree = plan->GetPatchTreeData().data();

    return *reinterpret_cast<float const*>(tree + sharpnessOffset());
}

inline SubdivisionPlan::Node
SubdivisionPlan::Node::operator ++ () {

    treeOffset += getNodeSize(GetDescriptor());

    if (treeOffset >= (int)GetSubdivisionPlan()->GetPatchTreeData().size())
        *this = {};

    return *this;
}

inline bool
SubdivisionPlan::Node::operator == (SubdivisionPlan::Node const & other) const {
    return (plan == other.plan) && (treeOffset == other.treeOffset);
}


} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_SUBDIVISION_PLAN_H */
