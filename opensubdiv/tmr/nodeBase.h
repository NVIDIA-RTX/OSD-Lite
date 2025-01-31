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

#ifndef OPENSUBDIV3_TMR_NODE_BASE_H
#define OPENSUBDIV3_TMR_NODE_BASE_H

#ifdef __cplusplus
namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {
namespace Tmr {
#endif

// Base class for SubdivisionPlan patch tree nodes. The goal is to push the decoding
// logic of tree nodes into a single cross-API header that can be included in both
// host & device source code.
//
// The tree contains 4 different types of nodes:
//
//   * Regular nodes (REG): describe a sub-domain of the limit surface that consists
//     of a 'regular' sub-patch.
//
//   * End / Extraordinary nodes (END): describe the surface sub-domain around
//     isolated irregular features. All "End-cap" patches in a subdivision plan 
//     must be of the same type (bilinear, 'regular' or 'gregory').
//
//   * Recursive nodes (REC): split the domain into 4 child sub-domains.
//     (see the domain winding notes)
//
//   * Terminal nodes (TER): an optimization for the sub-patches tree that allows
//     the collapsing of 4 children nodes into a single "terminal" node. Faces 
//     that contain a single extraordinary vertex, but are otherwise "regular"
//     (no boundaries, no creases...) generate 3 regular sub-patches for each level
//     of adaptive isolation.
//     Each of these regular sub-patches require a set patch points, however many
//     rows & columns overlap. Ex: out of 48 catmark patch points, only 24 are not
//     redundant. Terminal nodes store those patch point indices more efficiently.
//
// Node layouts by type :
//
//   Field            | Size | Content
//  ------------------|:----:|----------------------------------------------------
//  REGULAR           |      |
//  0. descriptor     | 1    | see NodeDescriptor
//  1. patch points   | 1    | offset to the first patch point index
//  2. sharpness      | 1    | crease sharpness (single crease nodes only)
//  ------------------|:----:|----------------------------------------------------
//  EMD               |      |
//  0. descriptor     | 1    | see NodeDescriptor
//  1. patch points   | 1    | offset to the first patch point index
//  ------------------|:----:|----------------------------------------------------
//  RECURSIVE         |      |
//  0. descriptor     | 1    | see NodeDescriptor
//  1. patch points   | 1    | offset to the first patch point index (dyn. isolation only)
//  2. child nodes    | 4    | offsets to the 4 children nodes
//  ------------------|:----:|----------------------------------------------------
//  TERMINAL          |      |
//  0. descriptor     | 1    | see NodeDescriptor
//  1. patch points   | 1    | offset to the first patch point index
//  2. child node     | 1    | offset to the child node (END or TERMINAL)
//  ------------------|:----:|----------------------------------------------------
//  (sizes are in 'integers' format, ie. 4 bytes)
//
// Notes:
//
//   * Domain winding
//     Sub-domain winding in topology trees does not follow the general
//     sequential winding order used elsewhere in OpenSubdiv. Instead, the
//     domains of sub-patches are stored with a "Z" winding order.
// 
//     Winding patterns:
//       Sequential     Z-order          Sequential        Z-order
//       +---+---+     +---+---+              +               +
//       | 3 | 2 |     | 2 | 3 |             /2\             /3\
//       +---+---+     +---+---+            /___\           /___\
//       | 0 | 1 |     | 0 | 1 |           /0\3/1\         /0\2/1\
//       +---+---+     +---+---+          +---+---+       +---+---+
// 
//     This winding pattern allows for faster traversal by using simple
//     '^' bitwise operators in the Quad case.
//
//  For more details, see:
//  "Efficient GPU Rendering of SUbdivision Surfaces using Adaptive Quadtrees"
//    W. Brainer, T. Foley, M. Mkraemer, H. Moreton, M. Niessner - Siggraph 2016
//
//  http://www.graphics.stanford.edu/~niessner/papers/2016/4subdiv/brainerd2016efficient.pdf
//
//

// XXXX manuelk: look into possible memory alignment gains (ie. node sizes are either 4 or 8 ui32)

struct NodeBase {

    static constexpr int maxIsolationLevel() { return 10; }

    // patch points

    static constexpr int catmarkRegularPatchSize() { return 16; };
    static constexpr int catmarkIrregularPatchSize(EndCapType type);
    static constexpr int catmarkTerminalPatchSize() { return 25; };

    static constexpr int loopRegularPatchSize() { return 12; };
    static constexpr int loopIrregularPatchSize(EndCapType type);

    // node sizes (in 'ints', not bytes)

    static constexpr int regularNodeSize(bool singleCrease) { return singleCrease ? 3 : 2; }
    static constexpr int endCapNodeSize() { return 2; }
    static constexpr int terminalNodeSize() { return 3; }
    static constexpr int recursiveNodeSize() { return 6; }

    static int getNodeSize(NodeType type, bool singleCrease);
    static int getNodeSize(NodeDescriptor desc);

    static constexpr int getNumChilren(NodeType type) {
        switch (type) {
            case NodeType::NODE_TERMINAL: return 1;
            case NodeType::NODE_RECURSIVE: return 4;
            default: return 0;
        }
    }

    static constexpr int rootNodeOffset() { return 14; }

    // internal node offsets in tree array

    constexpr int descriptorOffset() const { return treeOffset; }
    constexpr int sharpnessOffset() const { return treeOffset + 2; }
    constexpr int patchPointsOffset() const { return treeOffset + 1; }
    constexpr int childOffset(int childIndex) const { return treeOffset + 2 + childIndex; }

#ifdef __cplusplus
    // winding operators : both Loop & Catmark quadrant traversals
    // expect Z-curve winding (see subdivisionPlanBuilder for details)
    struct CatmarkWinding {
        inline unsigned char operator()(float& u, float& v);
        unsigned char quadrant = 0;

    };
    struct LoopWinding {
        inline unsigned char operator()(float& u, float& v);
        bool rotated = false;
        float median = 1.f;
    };

    // traverse the tree & returns the node & rotation quadrant for the given (s,t)
    // note: maxLevel will stop the traversal at the given tree depth
    template <typename Winding> static inline NodeBase traverse(uint32_t const* tree,
        float s, float t, unsigned char& quadrant, unsigned short maxLevel = maxIsolationLevel());
#endif
    int treeOffset = -1;
};


inline int
NodeBase::getNodeSize(NodeType type, bool singleCrease) {
    switch (type) {
        case NodeType::NODE_REGULAR: return regularNodeSize(singleCrease);
        case NodeType::NODE_END: return endCapNodeSize();
        case NodeType::NODE_TERMINAL: return terminalNodeSize();
        case NodeType::NODE_RECURSIVE: return recursiveNodeSize();
    }
    return 0;
}

inline int
NodeBase::getNodeSize(NodeDescriptor desc) {
    return getNodeSize(desc.GetType(), desc.HasSharpness());
}

constexpr int
NodeBase::catmarkIrregularPatchSize(EndCapType type) {
    switch (type) {
        case EndCapType::ENDCAP_BILINEAR_BASIS: return 4;
        case EndCapType::ENDCAP_BSPLINE_BASIS: return 16;
        case EndCapType::ENDCAP_GREGORY_BASIS: return 20;
        case EndCapType::ENDCAP_NONE: return 0;
    }
    return 0;
}

constexpr int
NodeBase::loopIrregularPatchSize(EndCapType type) {
    switch (type) {
        case EndCapType::ENDCAP_BILINEAR_BASIS: return 3;
        case EndCapType::ENDCAP_BSPLINE_BASIS: return 12;
        case EndCapType::ENDCAP_GREGORY_BASIS: return 18;
        case EndCapType::ENDCAP_NONE: return 0;
    }
    return 0;
}

#ifdef __cplusplus
template<typename Winding>
inline NodeBase NodeBase::traverse(
    uint32_t const* tree, float s, float t, unsigned char& quadrant, unsigned short maxLevel)
{
    using enum NodeType;

    NodeBase node = { .treeOffset = NodeBase::rootNodeOffset() };
    NodeDescriptor desc = { .field0 = tree[node.descriptorOffset()] };
    NodeType type = desc.GetType();

    Winding applyWinding;

    while (type == NODE_RECURSIVE || type == NODE_TERMINAL) {

        if ((short)desc.GetDepth() == maxLevel)
            break;

        quadrant = applyWinding(s, t);

        if (type == NODE_RECURSIVE) {
            // traverse to child node
            node = { .treeOffset = (int)tree[node.childOffset(quadrant)] };
        } else if (type == NODE_TERMINAL) {
            if (quadrant == desc.GetEvIndex()) {
                // traverse to end-cap patch (stored in the child node)
                node = { .treeOffset = (int)tree[node.childOffset(0)] };
            } else
                // regular sub-patch : exit
                break;
        }
        desc = { .field0 = tree[node.descriptorOffset()] };
        type = desc.GetType();
    }
    return node;
}

// Z-curve catmark winding : use bit operators to wind around the quadrants
//
//   1,0       1,1
//     +---+---+
//     | 2 | 3 |
//     +---+---+
//     | 0 | 1 |
//     +---+---+
//   0,0       1,0
//
inline unsigned char NodeBase::CatmarkWinding::operator()(float& u, float& v) {
    if (u >= 0.5f) { quadrant ^= 1; u = 1.0f - u; }
    if (v >= 0.5f) { quadrant ^= 2; v = 1.0f - v; }
    u *= 2.0f;
    v *= 2.0f;
    return quadrant;
}

// note: Z-winding of triangle faces rotates sub-domains every subdivision level,
// but the center face is always at index (2)
// 
//                0,1                                    0,1
//                 *                                      *
//               /   \                                  /   \
//              /     \                                /  3  \
//             /       \                              /       \
//            /         \           ==>        0,0.5 . ------- . 0.5,0.5
//           /           \                          /   \ 2 /   \
//          /             \                        /  0  \ /  1  \
//         * ------------- *                      * ----- . ----- *
//      0,0                 1,0                0,0      0.5,0      1,0
//
// XXXX manuelk: optimize this with bit operators too
inline unsigned char NodeBase::LoopWinding::operator()(float& u, float& v) {
    median *= 0.5f;
    if (!rotated) {
        if (u >= median) { u -= median; return 1; }
        if (v >= median) { v -= median; return 3; }
        if ((u + v) >= median) { rotated = true; return 2; }
    } else {
        if (u < median) { v -= median; return 1; }
        if (v < median) { u -= median; return 3; }
        u -= median;
        v -= median;
        if ((u + v) < median) { rotated = true; return 2; }
    }
    return 0;
}

} // end namespace Tmr
} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
#endif

#endif // OPENSUBDIV3_TMR_NODE_BASE_H
