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


#ifndef OPENSUBDIV3_TMR_TOPOLOGY_MAP_H
#define OPENSUBDIV3_TMR_TOPOLOGY_MAP_H

#include "../version.h"

#include "../tmr/subdivisionPlan.h"
#include "../tmr/types.h"

#include "../sdc/options.h"
#include "../sdc/types.h"

#include <cstdint>
#include <memory>
#include <shared_mutex>
#include <vector>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

class Neighborhood;

///
/// \brief Hashed dictionary of subdivision plans
///
/// Subdivision plans are indexed based on the binary hash of the 1-ring topological
/// neighborhood around the base face. The hash function is very cheap, but cannot
/// guarantee the absence of collisions. The open-addressing hash table resolves
/// index collisions through systematic bit-wise comparison of the topological data
/// between two Neighborhoods.
/// 
/// When a collision is detected, the open-addressing scheme simply moves forward
/// sequentially until an empty address is found. Assuming the map is not over-
/// populated, the performance impact of probe sequence length (PSL) on searches
/// should remain minimal (see implementation for more details).
/// 
/// note: irregular faces generate 'n' plans, where 'n' corresponds to the number of
/// vertices of the face. Only the first plan (corresponding to the first sub-face)
/// inserts its neighborhood hashes into the map. The remaining plans are expected to 
/// be apended in sequential counter-clockwise order.
/// 
class TopologyMap {

public:

    /// 'Traits' define the unique combination of subdivision rules & settings
    /// that the topology map can represent with a given set of subdivision plans.
    /// Some scheme options and settings are incompatible, so multiple maps may be
    /// necessary to cache all possible topology configurations.
    /// Exclusionary traits include: the subdivision scheme, the Chaikin rule and
    /// use of 'smooth triangle' rule.
    /// Subdivision plan stencil matrices can be generated so that face-varying
    /// boundary interpolation rules are matched to 'vertex' rules. Traits can
    /// automatically select a set of 'vertex' rules that maximizes the re-use of
    /// plans.
    struct Traits {

        /// Selects a set of traits compatible with the input parameters ; set
        /// 'faceVarying' to true to consider face-varying boundary interpolation
        /// rules, otherwise the 'vertex' boundary rules are used for a match.
        void SetCompatible(Sdc::SchemeType schemeType,
            Sdc::Options schemeOptions, EndCapType endcapType, bool faceVarying = false);

        Sdc::SchemeType getSchemeType() const;
        void setSchemeType(Sdc::SchemeType scheme);

        Sdc::Options::VtxBoundaryInterpolation getVtxBoundaryInterpolation() const;
        void setVtxBoundaryInterpolation(Sdc::Options::VtxBoundaryInterpolation b);
        void setVtxBoundaryInterpolation(Sdc::Options::FVarLinearInterpolation b);

        Sdc::Options::CreasingMethod getCreasingMethod() const;
        void setCreasingMethod(Sdc::Options::CreasingMethod c);

        Sdc::Options::TriangleSubdivision getTriangleSubdivision() const;
        void setTriangleSubdivision(Sdc::Options::TriangleSubdivision t);

        EndCapType getEndCapType() const;
        void setEndCapType(EndCapType t);

        uint8_t value() const {
            static_assert(sizeof(Traits) == sizeof(uint8_t));
            return *reinterpret_cast<uint8_t const*>(this); 
        };

        uint8_t scheme : 1;
        uint8_t vtxBoundInterp : 2;
        uint8_t creasingMethod : 1;
        uint8_t triangleSub : 1;
        uint8_t endcapType : 2;
    };

    struct Options {
        
        // XXXX mk : this constructor cannot be defaulted until default 
        // aggregate consutructors are fixed in gcc/clang
        Options(uint32_t id = 0) : uniqueID(id)  {}

        uint32_t const uniqueID : 5 = 0; // optional unique identifier if multiple topology
                                         // maps are being used ; must be less than
                                         // Tmr::SurfaceDescriptor::kMaxMapIndex

        uint32_t initialAddressSpace = 10000; // initial size of the address table

        uint32_t maxAddressSpace = ~uint32_t(0); // limits the maximum size of the address table:
                                                 // any attempt at inserting further plans will
                                                 // return INDEX_INVALID

        float loadFactorThreshold = .75f; // load threshold triggering the expansion of the
                                          // address table (re-indexes all hash entries ;
                                          // see implementation notes for more details)
    };

    /// \brief Constructor
    TopologyMap(Traits traits, Options options = Options());

    Traits const& GetTraits() const { return _traits; }

    //// \brief Returns the configuration options
    Options const& GetOptions() const { return _options; }

    /// \brief Returns the number of subdivision plans stored in the map
    int GetNumSubdivisionPlans() const { return (int)_plansTable.size(); }

    /// \brief Returns the subdivision plan at 'index' (random access in constant O(1) time)
    /// note: can return nullptr (see kRegularPlanAtIndexZero)
    /// note: this function is implictly lock-free as the hash map has no removal function
    SubdivisionPlan const* GetSubdivisionPlan(Index planIndex) const;

    /// \brief Returns the maximum number of patch points of any plan in the map
    int GetNumPatchPointsMax() const { return _numPatchPointsMax; }

public:

    /// \brief Searches the topology dictionary for a configuration matching the neighborhood
    /// and returns a plan index (INDEX_INVALID if none exists) along with the corresponding
    /// 'starting edge' value
    template <bool thread_safe = false> Index FindSubdivisionPlan(Neighborhood const& n, int& startingEdge) const;

    /// \brief Inserts the subdivision plan for a regular face
    template <bool thread_safe = false> Index InsertRegularFace(SubdivisionPlan* plan);

    /// \brief Inserts a sequence of subdivision plans for the sub-faces of an irregular face
    template <bool thread_safe = false> Index InsertIrregularFace(SubdivisionPlan** plans, int numPlans);

public:

    /// \brief Returns the sum of the subdivision plans sizes in bytes
    size_t GetByteSize(bool includeNeighborhood = false) const;

    struct HashmapStats {
        uint32_t pslMin = 0; // probe sequence length
        uint32_t pslMax = 0;
        float pslMean = 0.f;
        size_t hashCount = 0;
        size_t addressCount = 0;
        float loadFactor = 0.f;
    };

    /// \brief Traverse the hash addressing table to compute mean PSL & other statistics
    HashmapStats ComputeHashTableStatistics() const;

    /// \brief Writes a GraphViz 'dot' diagraph of all the subdivision plan trees
    void WriteSubdivisionPlansDigraphs(FILE* fout, char const* title, bool showIndices=true) const;

    typedef uint32_t address_type;
    typedef uint32_t key_type;

    // If true, forces the subdivision plan corresponding to pure regular topology to
    // always be inserted at PlanIndex = 0 (helps implementing dedicated code fast-paths)
    // Note that if no such configuration is hashed, GetSubdivisionPlan(0) will return nullptr
    static constexpr bool const kRegularPlanAtIndexZero = true;

private:

    // note: all functions within the private interface are lock-free (*not* thread-safe)

    static constexpr address_type const kOpenAddress = ~address_type(0);
    static constexpr size_t const kExpansionFactor = 2;

    float loadFactor() const;
    void resizeAddressSpace(size_t size);
    address_type findOpenAddress(key_type hashKey);
    void insertHashKeys(Index planIndex);

    // public interface implementation
    Index findPlan(Neighborhood const& neighborhood, int& startingEdge) const;   
    Index appendPlan(SubdivisionPlan* plan, bool isRegular);    
    Index insertRegularFace(SubdivisionPlan* plan);
    Index insertIrregularFace(SubdivisionPlan** plans, int numPlans);

private:

    Traits _traits = {};
    Options _options;

    mutable std::shared_mutex _mutex;

    int _numPatchPointsMax = 0;

    address_type _hashCount = kRegularPlanAtIndexZero;
    std::vector<address_type> _hashTable;

    std::vector<std::unique_ptr<SubdivisionPlan const>> _plansTable;
};

inline void 
TopologyMap::Traits::SetCompatible(Sdc::SchemeType schemeType,
    Sdc::Options schemeOptions, EndCapType endcapType, bool faceVarying) {

    setSchemeType(schemeType);   
    if (faceVarying)
        setVtxBoundaryInterpolation(schemeOptions.GetFVarLinearInterpolation());
    else
        setVtxBoundaryInterpolation(schemeOptions.GetVtxBoundaryInterpolation());
    setCreasingMethod(schemeOptions.GetCreasingMethod());    
    setTriangleSubdivision(schemeOptions.GetTriangleSubdivision());    
    setEndCapType(endcapType);
}

inline void
TopologyMap::Traits::setSchemeType(Sdc::SchemeType schemeType) {
    static_assert(Sdc::SchemeType::SCHEME_BILINEAR == Sdc::SchemeType(0));
    scheme = uint8_t(schemeType) - 1;
}
inline Sdc::SchemeType 
TopologyMap::Traits::getSchemeType() const {
    return Sdc::SchemeType(scheme + 1);
}

inline void
TopologyMap::Traits::setVtxBoundaryInterpolation(Sdc::Options::VtxBoundaryInterpolation i) {
    vtxBoundInterp = uint8_t(i);
}
inline void
TopologyMap::Traits::setVtxBoundaryInterpolation(Sdc::Options::FVarLinearInterpolation i) {
    using enum Sdc::Options::VtxBoundaryInterpolation;
    using enum Sdc::Options::FVarLinearInterpolation;
    // Attempt to match face-varying boundary rules with a set of vertex boundary rules
    // that maximizes plan re-use. Differences in boundary rules are handled by the builder
    // by sharpening edges and vertices as necessary when hashing face-varying topology.
    // note: Corner+1 and Corner+2 are currently not supported yet as matching their properties
    // is more complex, but not impossible.
    switch (i) {
        case FVAR_LINEAR_NONE:
            vtxBoundInterp = uint8_t(VTX_BOUNDARY_EDGE_ONLY); break;
        case FVAR_LINEAR_CORNERS_ONLY:
        case FVAR_LINEAR_BOUNDARIES: 
            vtxBoundInterp = uint8_t(VTX_BOUNDARY_EDGE_AND_CORNER); break;
        default:
            assert("face-varying boundary interpolation mode not supported"==0);
    }
}
inline Sdc::Options::VtxBoundaryInterpolation 
TopologyMap::Traits::getVtxBoundaryInterpolation() const {
    return Sdc::Options::VtxBoundaryInterpolation(vtxBoundInterp);
}

inline void
TopologyMap::Traits::setCreasingMethod(Sdc::Options::CreasingMethod c) {
    creasingMethod = uint8_t(c);
}
inline Sdc::Options::CreasingMethod 
TopologyMap::Traits::getCreasingMethod() const {
    return Sdc::Options::CreasingMethod(creasingMethod);
}


inline Sdc::Options::TriangleSubdivision 
TopologyMap::Traits::getTriangleSubdivision() const {
    return Sdc::Options::TriangleSubdivision(triangleSub);
}
inline void 
TopologyMap::Traits::setTriangleSubdivision(Sdc::Options::TriangleSubdivision t) {
    triangleSub = Sdc::Options::TriangleSubdivision(t);
}

inline void
TopologyMap::Traits::setEndCapType(EndCapType t) {
    endcapType = uint8_t(t);
}
inline EndCapType
TopologyMap::Traits::getEndCapType() const {
    return EndCapType(endcapType);
}

inline SubdivisionPlan const*
TopologyMap::GetSubdivisionPlan(Index planIndex) const {

    // read-lock here

    return _plansTable[planIndex].get();
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_TOPOLOGY_MAP_H */
