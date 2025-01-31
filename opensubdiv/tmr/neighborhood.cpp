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

#include "../tmr/neighborhood.h"

#include "../far/topologyDescriptor.h"
#include "../far/topologyRefiner.h"
#include "../far/topologyRefinerFactory.h"
#include "../sdc/options.h"
#include "../sdc/types.h"

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

using hashkey_type = Neighborhood::hashkey_type;

std::array<uint8_t, 196> const Neighborhood::_regularCatmarkData = []() {
    
    // Pur regular Catmark neighborhood: no irregular vertex, no boundaries, no tags.
    //
    //      15 ------ 14 ------ 13 ------ 12
    //       |         |         |         |
    //       |    8    |    7    |    6    |
    //       |         |         |         |
    //       4 ------  3 ------  2 ------ 11
    //       |         |         |         |
    //       |    1    |    0    |    5    |
    //       |         |         |         |
    //       5 ------  0 ------  1 ------ 10
    //       |         |         |         |
    //       |    2    |    3    |    4    |
    //       |         |         |         |
    //       6 ------  7 ------  8 ------  9
    //    
    // This particular topology is canonical of the topology traversal implemented
    // in the NeighborhoodBuilder and must be revised if the algorithm changes!

    constexpr int const faceVertCounts[] = { 4, 4, 4, 4, 4, 4, 4, 4, 4 };

    constexpr int const faceVerts[] = {
        0,  1,  2,  3, /**/ 0,  3,  4,  5, /**/ 0,  5,  6, 7,
        0,  7,  8,  1, /**/ 1,  8,  9, 10, /**/ 1, 10, 11, 2,
        2, 11, 12, 13, /**/ 2, 13, 14,  3, /**/ 3, 14, 15, 4, };

    std::array<uint8_t, 196> neighborhoodData;

    populateData(neighborhoodData.data(), neighborhoodData.size(), {
        .faceVerts = faceVerts, 
        .nfaceVerts = (int)std::size(faceVerts),
        .faceVertCounts = faceVertCounts, 
        .nfaceVertCounts = (int)std::size(faceVertCounts), 
    });

    assert(reinterpret_cast<Neighborhood const*>(neighborhoodData.data())->GetHash() == 0x6D2C4B37);

    //reinterpret_cast<Neighborhood const*>(neighborhoodData.data())->Print();

    return neighborhoodData;
}();

std::array<uint8_t, 224> const Neighborhood::_regularLoopData = []() {

    // Pure regular loop neighborhood: no irregular vertex, no boundaries, no tags.
    //
    //               11 --------- 10
    //             /    \       /    \ 
    //            /  12  \  11 /  10  \       
    //           /        \   /        \        
    //          3 --------- 2 --------- 9 
    //        /   \       /   \       /   \
    //       /  2  \  1  /  0  \  9  /  8  \      
    //      /       \   /       \   /       \      
    //     4 -------- 0 --------- 1 -------- 8
    //      \       /   \       /   \       /      
    //       \  3  /  4  \  5  /  6  \  7  /      
    //        \   /       \   /       \   /
    //          5 --------- 6 --------- 7      
    //     
    // This particular topology is canonical of the topology traversal implemented
    // in the NeighborhoodBuilder and must be revised if the algorithm changes!

    constexpr int faceVertCounts[] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

    constexpr int faceVerts[] = {
        0,  1,  2, /**/  0,  2,  3, /**/ 0,  3,  4,
        0,  4,  5, /**/  0,  5,  6, /**/ 0,  6,  1,
        1,  6,  7, /**/  1,  7,  8, /**/ 1,  8,  9,
        1,  9,  2, /**/  2,  9, 10, /**/ 2, 10, 11,
        2, 11,  3
    };

    std::array<uint8_t, 224> neighborhoodData;

    populateData(neighborhoodData.data(), neighborhoodData.size(), {
        .faceVerts = faceVerts,
        .nfaceVerts = (int)std::size(faceVerts),
        .faceVertCounts = faceVertCounts,
        .nfaceVertCounts = (int)std::size(faceVertCounts),
        });

    assert(reinterpret_cast<Neighborhood const*>(neighborhoodData.data())->GetHash() == 0x30D16FAF);

    //reinterpret_cast<Neighborhood const*>(neighborhoodData.data())->Print();

    return neighborhoodData;
}();

void 
Neighborhood::populateData(uint8_t* dest, size_t dest_size, DataDescriptor const& desc)
{
    Neighborhood* n = reinterpret_cast<Neighborhood*>(dest);

    assert(dest_size >= byteSize(desc.nfaceVertCounts, desc.nfaceVerts, desc.ncorners, desc.ncreases, desc.ncontrolPoints));

    n->_faceCount = desc.nfaceVertCounts;
    n->_faceVertsCount = desc.nfaceVerts;
    n->_cornersCount = desc.ncorners;
    n->_creasesCount = desc.ncreases;
    n->_controlPointsCount = desc.ncontrolPoints;

    auto populate = []<typename T>(T* dest, T const* src, int src_count) {
        std::memcpy(dest, src, src_count * sizeof(T));
    };

    populate(n->faceVertCountsPtr(), desc.faceVertCounts, desc.nfaceVertCounts);
    populate(n->faceVertsPtr(), desc.faceVerts, desc.nfaceVerts);
    if (desc.ncorners > 0) {
        populate(n->cornerVertsPtr(), desc.cornerVerts, desc.ncorners);
        populate(n->cornerSharpnessPtr(), desc.cornerSharpness, desc.ncorners);
    }
    if (desc.ncreases) {
        populate(n->creaseVertsPtr(), desc.creaseVerts, desc.ncreases * 2);
        populate(n->creaseSharpnessPtr(), desc.creaseSharpness, desc.ncreases);
    }
    if (desc.ncontrolPoints)
        populate(n->controlPointsPtr(), desc.controlPoints, desc.ncontrolPoints);

    n->_hash = n->computeHashKey();
    n->_startingEdge = desc.startingEdge;
}

// actual size of a Neighborhood in bytes
uint32_t
Neighborhood::byteSize(
    int faceCount, int faceVertsCount, int cornersCount, int creasesCount, int controlPointsCount) {

    int size = sizeof(Neighborhood);
    size += faceCount * sizeof(int);
    size += faceVertsCount * sizeof(int);
    size += cornersCount * sizeof(int) + cornersCount * sizeof(float);
    size += creasesCount * 2 * sizeof(int) + creasesCount * sizeof(float);
    size += controlPointsCount * sizeof(int);
    return size;    
}

uint32_t
Neighborhood::GetByteSize() const {
    return byteSize(_faceCount, _faceVertsCount, _cornersCount, _creasesCount, _controlPointsCount);
}

inline hashkey_type
hashBytes(void const* bytes, size_t size) {

    // FNV-1 hash (would FNV1-a be better ?)
    uint32_t hvalue = 0;
    uint8_t const* start = (unsigned char const*)bytes;
    uint8_t const* end = start + size;
    while (start<end) {
        hvalue *= (uint32_t)0x01000193;
        hvalue ^= (uint32_t)*start;
        ++start;
    }
    return hvalue;
}

hashkey_type
Neighborhood::computeHashKey() const {

    // topology information starts with the first component count (_faceCount)
    // and ends at the last component value (crease sharpness values)

    uint8_t const* start = byteStart() + offsetof(Neighborhood, _faceCount);
    uint8_t const* end = byteEnd();
    
    size_t size = end - start;

    return hashBytes(start, size);
}

std::unique_ptr<Far::TopologyRefiner> 
Neighborhood::CreateRefiner(int schemeType, Sdc::Options const& schemeOptions) const {

    if (!HasControlPoints())
        return nullptr;

    Far::TopologyDescriptor descr;
    descr.numVertices = GetNumControlPoints();
    descr.numFaces = GetNumFaces();
    descr.numVertsPerFace = GetFaceVertCounts().begin();
    descr.vertIndicesPerFace = GetFaceVerts().begin();

    if (HasSharpVertices()) {

        Far::ConstIndexArray cornerIndices = GetCornerVertIndices();
        Vtr::ConstArray<float> cornerSharpness = GetCornerSharpness();

        int numCorners = cornerSharpness.size();
        assert(cornerIndices.size() == numCorners);

        descr.numCorners = numCorners;
        descr.cornerVertexIndices = cornerIndices.begin();
        descr.cornerWeights = cornerSharpness.begin();
    }

    if (HasSharpEdges()) {

        Far::ConstIndexArray creaseIndices = GetCreaseEdgeIndices();
        Vtr::ConstArray<float> creaseSharpness = GetCreaseSharpness();

        int numCreases = creaseSharpness.size();
        assert(creaseIndices.size() == (numCreases * 2));

        descr.numCreases = numCreases;
        descr.creaseVertexIndexPairs = creaseIndices.begin();
        descr.creaseWeights = creaseSharpness.begin();
    }

    using RefinerFactory = Far::TopologyRefinerFactory<Far::TopologyDescriptor>;

    RefinerFactory::Options refinerOptions(Sdc::SchemeType(schemeType), schemeOptions);

    return std::unique_ptr<Far::TopologyRefiner>(RefinerFactory::Create(descr, refinerOptions));
}

void
Neighborhood::Print() const {

    auto dumpArray = []<typename T>(char const* name, Vtr::ConstArray<T> const& a) {
        printf("%s = (%d) { ", name, a.size());
        for (int i = 0; i < a.size(); ++i) {
            if (i > 0) {
                printf(", ");
            }
            if constexpr (std::is_same<T, Vtr::Index>::value)
                printf("%d", a[i]);
            else
                printf("%.3f", a[i]);
        }
        printf(" }\n");
    };

    printf("Neighborhood hash = 0x%X edge=%d {\n", GetHash(), GetStartingEdge());
    dumpArray("\tface vert counts", GetFaceVertCounts());
    dumpArray("\tface verts", GetFaceVerts());
    dumpArray("\tcorner verts", GetCornerVertIndices());
    dumpArray("\tcorner sharpness", GetCornerSharpness());
    dumpArray("\tcrease verts", GetCreaseEdgeIndices());
    dumpArray("\tcrease sharpness", GetCreaseSharpness());
    if (HasControlPoints()) {
        
        ConstIndexArray remapTable(controlPointsPtr(), _controlPointsCount);
        auto faceVerts = GetFaceVerts();
        printf("\tcontrol points = (%d) { ", faceVerts.size());
        for (int i = 0; i < faceVerts.size(); ++i) {
            if (i > 0)
                printf(", ");
            printf("%d", remapTable[faceVerts[i]]);
        }
        printf(" }\n");
        dumpArray("\tremap table", remapTable);
    }
    printf("}\n");

    fflush(stdout);
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
