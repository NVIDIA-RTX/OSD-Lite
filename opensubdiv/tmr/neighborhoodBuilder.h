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

#ifndef OPENSUBDIV3_TMR_NEIGHBORHOOD_BUILDER_H
#define OPENSUBDIV3_TMR_NEIGHBORHOOD_BUILDER_H

#include "../version.h"

#include "../tmr/neighborhood.h"
#include "../far/topologyRefiner.h"
#include "../vtr/stackBuffer.h"

#include <cassert>
#include <memory>
#include <vector>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Vtr::internal {
    class Level;
    class FVarLevel;
}

namespace Tmr {

///
/// \brief Neighborhood build helper
///
/// Traverses a Vtr::Level mesh to gather the topology components of the 1-ring
/// around a given face. Topology parsing requires generating at least one
/// Neighborhood for each base face, so the builders maintains a persistent
/// cache of scratch memory in order to avoid high-frequency allocations / de-
/// allocations. 
/// 
/// The public interface also exposes both a creation mode that allocates heap memory,
/// and a 'populate' mode that can output the Neighborhood data into a Vtr::StackBuffer.
/// The latter is useful if the neighborhood data is transient and only used for a fast
/// topology comparison against an existing dictionary.
/// 
/// note: it may eventually be necessary to expand the builder to support a memory-pool
/// allocator. Neighborhoods are generally relatively small data packets and a large
/// dictionary may allocate large quantities, which could cause fragmentation.
/// 
class NeighborhoodBuilder {

public:

    using Level = Vtr::internal::Level;
    using FVarLevel = Vtr::internal::FVarLevel;
    template <unsigned int SIZE> using StackBuffer = Vtr::internal::StackBuffer<uint8_t, SIZE>;

    static constexpr bool const kDebugNeighborhoods = false;

public:


    /// The expected maximum valence in a neighborhood is used as a guess
    /// to initialize persistent scratch memory. The value does not have to be
    /// particularly accurate, as potential reallocations should be rare and
    /// amortized over re-uses of the builder.
    NeighborhoodBuilder(int maxValence=128);

    struct BuildDescriptor {
        Far::TopologyRefiner const& refiner;
        Index faceIndex = INDEX_INVALID;
        int startingEdge = 0;
        int fvarChannel = -1;
    };

    // Allocates from the heap & return the unique binary blob representing the
    // neighborhood around faceIndex (use reinterpret_cast<Neighborhood const*>
    // to decode the information)
    std::unique_ptr<uint8_t const []> Create(BuildDescriptor const& desc);

    /// Encodes the topology data of the neighborhood into a persistent stackbuffer.
    /// (returns the size of the neighborhood, 0 on failure)
    /// 
    /// note: the vast majority of neighborhoods fit in less than 512 bytes ; the 
    /// largest neighborhood in the OpenSubdiv regression suite (as of this writing) 
    /// uses about 9 Kb, and is caused by a vertex of valence 360. "512 bytes should
    /// be enough for everyone" (famous last words)
    template <unsigned int SIZE> size_t Populate(
        StackBuffer<SIZE>& data, BuildDescriptor const& desc, bool remapVerts = true);

private:

    // Vertex and edge sharpness of the outer 1-ring should have no influence on 
    // the limit surface and its derivatives, so we should be able to ignore them
    static constexpr bool const gatherOneRingCreaseAndCorners = false;

    void reserve(int maxValence);
    void clear();

    LocalIndex findLocalIndex(Index index);
    LocalIndex mapToLocalIndex(Index vertexIndex);
    LocalIndex mapToLocalIndex(Index vertexIndex, Index valueIndex);
    
    void addFace(ConstIndexArray fverts, Index faceIndex, int startingEdge);
    void addFace(ConstIndexArray fverts, ConstIndexArray fvalues, Index faceIndex, int startingEdge);
    template <bool quash=false> void addCrease(Index v0, Index v1, float sharpness);
    template <bool quash=false> void addCorner(Index v, float sharpness);

private:

    // vertex topology
    void addVertexFaces(Level const& level, 
        Index baseFace, Index vertIndex, int vertInFace, int regFaceSize, bool unordered);
    void gatherVertexTopology(BuildDescriptor const& desc);

    // face-varying topology
    bool addVertexFace(FVarLevel const& fvlevel,
        Index faceIndex, Index vertIndex, Index valueIndex);
    void gatherFVarTopology(BuildDescriptor const& desc);

private:

    // pack topology data into binary form
    size_t computeByteSize(bool remapVerts) const;
    void populateData(uint8_t* buf, size_t buf_size, int startingEdge, bool remapVerts) const;

    typedef std::pair<Index, Index> Index2;

    std::vector<int> _faceVertCounts;    

    std::vector<Index> _faces;
    std::vector<Index> _faceVerts;
    std::vector<Index> _vertRemaps;
    
    std::vector<Index> _valueVerts;

    std::vector<Index> _cornerVerts;
    std::vector<float> _cornerSharpness;

    std::vector<Index2> _creaseVerts;
    std::vector<float> _creaseSharpness;
};

std::unique_ptr<uint8_t const []> inline
NeighborhoodBuilder::Create(BuildDescriptor const& desc) {

    assert(desc.faceIndex != INDEX_INVALID && desc.startingEdge >= 0);

    if (desc.fvarChannel < 0)
        gatherVertexTopology(desc);
    else
        gatherFVarTopology(desc);

    constexpr bool remapVerts = kDebugNeighborhoods;

    std::unique_ptr<uint8_t[]> neighborhoodData;
    if (size_t byteSize = computeByteSize(remapVerts)) {
        neighborhoodData = std::make_unique<uint8_t[]>(byteSize);
        populateData(neighborhoodData.get(), byteSize, desc.startingEdge, remapVerts);
    }
    return neighborhoodData;
}

template <unsigned int SIZE> inline size_t NeighborhoodBuilder::Populate(
    StackBuffer<SIZE>& data, BuildDescriptor const& desc, bool remapVerts) {

    assert(desc.faceIndex != INDEX_INVALID && desc.startingEdge >= 0);

    if (desc.fvarChannel < 0)
        gatherVertexTopology(desc);
    else
        gatherFVarTopology(desc);

    size_t byteSize = computeByteSize(remapVerts);
    data.SetSize((unsigned int)byteSize);
    if (byteSize > 0)
        populateData(data, data.GetSize(), desc.startingEdge, remapVerts);
    return byteSize;
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_NEIGHBORHOOD_BUILDER_H */

