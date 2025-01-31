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

#ifndef OPENSUBDIV3_TMR_NEIGHBORHOOD_H
#define OPENSUBDIV3_TMR_NEIGHBORHOOD_H

#include "../version.h"

#include "../sdc/types.h"
#include "../tmr/types.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Sdc {
    class Options;
}

namespace Far {
    class TopologyRefiner;
}

namespace Tmr {

/// \brief Topoloogy neighborhood descriptor
/// 
/// Uniquely identifies the topological configuration of a 1-ring, including crease
/// tags. All component indices (face, verts, ...) are strictly local to the topology
/// of the neighborhood, which makes comparison for topological equivalence possible.
/// 
/// Neighborhoods can also optionally store an array of control point indices from
/// their mesh of origin. The local indices of an equivalent neighborhood can thus be
/// remapped to that mesh (specifically, this feature is used to map the columns 
/// of  the stencil matrix of a subdivision plan back to the control points of an 
/// arbirtrary mesh).
/// 
/// note: the amount of information stored in a neighborhood can be very variable, so 
/// we reduce memory fragmentation by allocating the total footprint of the data in one
/// block of bytes. The Neighborhood interface allows safe & accurate re-indexing of
/// the components data, but Neighborhood instances are not PODs !
///
/// Memory layout:
/// 
///   - header  'Neighborhood' struct
///   - int     faceValences[faceValencesCount]
///   - int     faceVerts[faceVertsCount]
///   - int     cornerIndices[cornersCount]
///   - float   cornerSharpness[cornersCount]
///   - int2    creaseIndices[creasesCount]
///   - float   creaseSharpness[creasesCount]
///   - int     controlPoints[controlPointsCount] (optional)
/// 
class Neighborhood {

public:

    /// \brief Returns the number of vertices of the center face
    int GetBaseFaceSize() const { return GetFaceVertCounts()[0]; }

    /// \brief Index of sub-face, if base-face is irregular
    int GetStartingEdge() const { return _startingEdge; }

    int GetNumFaces() const { return _faceCount; }
    int GetNumControlPoints() const { return _controlPointsCount; }
    ConstIndexArray GetFaceVertCounts() const;
    ConstIndexArray GetFaceVerts() const;

    bool HasSharpVertices() const { return _cornersCount > 0; }
    ConstFloatArray GetCornerSharpness() const;
    ConstIndexArray GetCornerVertIndices() const;

    bool HasSharpEdges() const { return _creasesCount > 0; }
    ConstFloatArray GetCreaseSharpness() const;
    ConstIndexArray GetCreaseEdgeIndices() const;

    bool HasControlPoints() const { return _controlPointsCount > 0; }
    ConstIndexArray GetControlPoints() const;

    /// \brief Remaps a vertex index to the local neighborhood
    /// note: control point indices are optional, so this will return
    /// INDEX_INVALID if they were not stored in this instance.
    LocalIndex FindLocalIndex(Index vertIndex) const;

    /// \brief Returns total memory use in bytes
    uint32_t GetByteSize() const;

    typedef uint32_t hashkey_type;

    /// \brief Returns a hash of the neighborhood data (using FNV-1). This
    /// hash is extremely fast & cheap, so collisions should be expected.
    hashkey_type GetHash() const { return _hash; }

    /// \brief Returns true if both neighborhoods are topologically equivalent
    bool IsEquivalent(Neighborhood const& other) const;

    /// \brief Returns true if the neighborhood is topologically equivalent
    /// to that of a pure regular surface (no irregular vertex, no boundaries,
    /// no tags).
    bool IsRegularNeighborhood(Sdc::SchemeType scheme) const;

    // \brief Returns the unique canonical neighborhood of a pure regular 
    // surface (no irregular vertex, no boundaries, no tags).
    static Neighborhood const* GetRegularNeighborhood(Sdc::SchemeType scheme);

    /// \brief Rebuild the complete topolgoy of the neighborhood ; requires that
    /// the neighborhood stores its control point indices (ie. HasControlPoints() 
    /// must return true)
    std::unique_ptr<Far::TopologyRefiner> CreateRefiner(
        int schemeType, Sdc::Options const& schemeOptions) const;

    // Debug printout
    void Print() const;

private:

    friend class NeighborhoodBuilder;

    // Neighborhood is a 'descriptor' and does not own its memory: prevent C++ "accidents"
    Neighborhood(Neighborhood&&) = delete;
    Neighborhood(Neighborhood const&) = delete;
    Neighborhood& operator = (Neighborhood const&) = delete;
    ~Neighborhood() { assert(false); }

    hashkey_type computeHashKey() const;

    int* faceVertCountsPtr() const { return (int *)(byteStart()+sizeof(Neighborhood)); }
    int* faceVertsPtr() const { return faceVertCountsPtr() + _faceCount; }

    int* cornerVertsPtr() const { return faceVertsPtr() + _faceVertsCount; }
    float* cornerSharpnessPtr() const { return (float*)(cornerVertsPtr() + _cornersCount); }

    int* creaseVertsPtr() const { return (int*)(cornerSharpnessPtr() + _cornersCount); }
    float* creaseSharpnessPtr() const { return (float*)(creaseVertsPtr() + (2 * _creasesCount)); }

    int* controlPointsPtr() const { return (int*)(creaseSharpnessPtr() + _creasesCount); }

    // total size in bytes of a neighborhood with the given component counts
    static uint32_t byteSize(int faceCount, int faceVertsCount, 
        int cornersCount, int creasesCount, int controlPointsCount);

    uint8_t const* byteStart() const; // first byte of the neighborhood data
    uint8_t const* byteEnd() const;   // last byte of the topology data (excludes control points)

    struct DataDescriptor {
        int const* faceVerts = nullptr;
        uint16_t nfaceVerts = 0;
        int const* faceVertCounts = nullptr;
        uint16_t nfaceVertCounts = 0;
        int const* cornerVerts = nullptr;
        float const* cornerSharpness = nullptr;
        uint16_t ncorners = 0;
        int const* creaseVerts = nullptr;
        float const* creaseSharpness = nullptr;
        uint16_t ncreases = 0;
        int const* controlPoints = nullptr;
        uint16_t ncontrolPoints = 0;
        uint16_t startingEdge = 0;
    };

    static void populateData(uint8_t* dest, size_t dest_size, DataDescriptor const& desc);

private:

    static std::array<uint8_t, 196> const _regularCatmarkData;
    static std::array<uint8_t, 224> const _regularLoopData;

    // Some information (like vertex indices) is not part of a neighborhood hash
    // and has to be excluded from hashing or comparisons. For simplicity, the 
    // mandatory neighborhood information is arranged to be contiguous in memory.
    // note : C++ guarantees the order of struct/class members in memory, but the 
    // compiler is free to adjust data sizes for memory alignment

    // optional fields (excluded from hashing)
    uint16_t _startingEdge;
    uint16_t _controlPointsCount;
    uint32_t _hash;

    // topological information (hashing starts here)
    uint16_t _faceCount;
    uint16_t _faceVertsCount;
    uint16_t _cornersCount;
    uint16_t _creasesCount;
};


inline ConstIndexArray
Neighborhood::GetFaceVertCounts() const {
    return ConstIndexArray(faceVertCountsPtr(), _faceCount);
}

inline ConstIndexArray
Neighborhood::GetFaceVerts() const {
    return ConstIndexArray(faceVertsPtr(), _faceVertsCount);
}

inline ConstIndexArray
Neighborhood::GetCornerVertIndices() const {
    return ConstIndexArray(cornerVertsPtr(), _cornersCount);
}

inline ConstFloatArray
Neighborhood::GetCornerSharpness() const {
    return Vtr::ConstArray<float>(cornerSharpnessPtr(), _cornersCount);
}

inline ConstIndexArray
Neighborhood::GetCreaseEdgeIndices() const {
    return ConstIndexArray(creaseVertsPtr(), _creasesCount * 2);
}

inline ConstFloatArray
Neighborhood::GetCreaseSharpness() const {
    return Vtr::ConstArray<float>(creaseSharpnessPtr(), _creasesCount);
}

inline ConstIndexArray
Neighborhood::GetControlPoints() const {
    return ConstIndexArray(controlPointsPtr(), _controlPointsCount);
}

inline bool
Neighborhood::IsEquivalent(Neighborhood const& other) const {

    // because we do no trust our hash function, run a bit-wise comparison
    // of the topology data defined in contiguous memory (excluding other
    // non-topological information)

    // note: the comparison starts with the bytes of the hash-key, which
    // functions as a trivial rejection mechanism
    
    size_t offset = offsetof(Neighborhood, _hash);

    uint8_t const* ptrA = this->byteStart() + offset;
    uint8_t const* ptrB = other.byteStart() + offset;

    size_t sizeA = this->byteEnd() - ptrA;
    size_t sizeB = other.byteEnd() - ptrB;
    if (sizeA != sizeB)
        return false;

    return std::memcmp(ptrA, ptrB, sizeA) == 0;
}

inline LocalIndex
Neighborhood::FindLocalIndex(Index vertIndex) const {
    // linear search is likely the best compromise here
    ConstIndexArray controlPoints = GetControlPoints();
    assert(controlPoints.size() > 0);
    return controlPoints.FindIndex(vertIndex);
}

inline uint8_t const* 
Neighborhood::byteStart() const {
    return reinterpret_cast<uint8_t const*>(this);
}

inline uint8_t const* 
Neighborhood::byteEnd() const {
    return byteStart() + byteSize(_faceCount, _faceVertsCount, _cornersCount, _creasesCount, 0);
}

inline bool
Neighborhood::IsRegularNeighborhood(Sdc::SchemeType scheme) const
{
    return IsEquivalent(*GetRegularNeighborhood(scheme));
}

inline Neighborhood const*
Neighborhood::GetRegularNeighborhood(Sdc::SchemeType scheme) {

    switch (scheme)
    {
        case Sdc::SchemeType::SCHEME_CATMARK:
            return reinterpret_cast<Neighborhood const*>(_regularCatmarkData.data());
        case Sdc::SchemeType::SCHEME_LOOP:
            return reinterpret_cast<Neighborhood const*>(_regularLoopData.data());
        case Sdc::SchemeType::SCHEME_BILINEAR:
        default:
            break;
    }
    assert(0);
    return nullptr;
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_NEIGHBORHOOD_H */

