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

#include "../tmr/topologyMap.h"
#include "../tmr/neighborhoodBuilder.h"

#include <cmath>
#include <mutex>
#include <numeric>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

using address_type = TopologyMap::address_type;
using key_type = TopologyMap::key_type;

//
// Given that perfect hashing of Neighborhood topology data is not really practical,
// TopologyMap implements collision resolution through a closed-hashing (or open-
// addressing) method (see https://en.wikipedia.org/wiki/Open_addressing)
// 
// Even though Neighborhoods cache the value of their hash-key, a cheap algorithm is
// preferred, as the performance scales linearly with the number of faces in a base
// mesh. In consequence, much of the search and insertion performance of TopologyMap
// is driven by the general coherence or clustering of hash key values within the
// address table.
// 
// Measuring the effects of hash key clustering on sample meshes seems to indicate
// that the mean PSL (probing sequence length) remains relatively low (less than 10
// iterations), even under table loads above 90%. As expected, mean PSL eventually
// rises sharply as the address space fills up.
// 
// This behavior would seem to confirm the choice of linear probing as the method 
// of choice for resolving address collisions, as opposed to more complex probing
// sequences (ex. quadratic, Hopscotch, Cuckoo, Robin Hood, ...)
// 
// The extremely cheap FNV-1 hash function used to hash the key values also seems
// to generate relatively unclustered values, although no substantial analysis has
// been made so far. 
// 
// Therefore, a reasonable compromise between performance and memory use should
// be relatively easy to achieve by adjusting the load threshold that triggers
// hash-key table resizing.
//

TopologyMap::TopologyMap(Traits traits, Options options) : _traits(traits), _options(options) {

    resizeAddressSpace(_options.initialAddressSpace);

    if constexpr (kRegularPlanAtIndexZero)
        _plansTable.emplace_back(nullptr);
}

//
// lock-free private interface : caller is expected to manage thread-safety
//

inline address_type
computeAddress(key_type hashKey, uint32_t psl, uint32_t tableSize) {
    return  (hashKey + psl) % tableSize;
}

inline address_type 
TopologyMap::findOpenAddress(key_type hashKey) {

    uint32_t tableSize = (uint32_t)_hashTable.size();

    for (uint32_t psl = 0; psl < tableSize; ++psl) {

        address_type address = computeAddress(hashKey, psl, tableSize);

        if (_hashTable[address] == kOpenAddress)
            return address;

        // address was already in use : keep searching
    }

    assert(0); // we cycled over the entire address space : the map is full !

    return kOpenAddress;
}

inline void
TopologyMap::insertHashKeys(Index planIndex) {

    if (SubdivisionPlan const* plan = _plansTable[planIndex].get()) {

        for (int i = 0; i < plan->GetNumNeighborhoods(); ++i) {
        
            address_type address = findOpenAddress(plan->GetNeighborhood(i).GetHash());

            assert(address != kOpenAddress);

            _hashTable[address] = planIndex;
            ++_hashCount;
        }
    }
}

void
TopologyMap::resizeAddressSpace(size_t size) {

    assert(size > _hashTable.size());

    _hashCount = 0;
    _hashTable.clear();
    _hashTable.resize(size, kOpenAddress);

    for (Index planIndex = 0; planIndex < (Index)_plansTable.size(); ) {

        insertHashKeys(planIndex);

        SubdivisionPlan const& plan = *_plansTable[planIndex];

        planIndex += plan.IsRegularFace() ? 1 : plan.GetFaceSize();
    }
}

inline float 
TopologyMap::loadFactor() const {
    return _hashTable.empty() ? 0.f : float(_hashCount) / float(_hashTable.size());
}

inline Index
TopologyMap::appendPlan(SubdivisionPlan* plan, bool indexZero) {

    assert(_plansTable.size() < std::numeric_limits<Index>::max());

    if (loadFactor() > _options.loadFactorThreshold) {
        
        size_t size = _hashTable.size() * kExpansionFactor;
        
        if (size > _options.maxAddressSpace)
            return INDEX_INVALID;

        resizeAddressSpace(size);
    }

    plan->_topologyMap = this;

    // include the control points in case users want to copy them in the same buffer
    int numPatchPoints = (int)plan->GetNumPatchPoints() + (int)plan->GetNumControlPoints();
    _numPatchPointsMax = std::max(_numPatchPointsMax, numPatchPoints);

    if constexpr (kRegularPlanAtIndexZero) {
        if (indexZero) {
            assert(!_plansTable[0]);
            _plansTable[0].reset(plan);
            return Index(0);
        }
    }

    _plansTable.emplace_back().reset(plan);

    return Index(_plansTable.size() - 1);
}

inline Index
TopologyMap::findPlan(Neighborhood const & neighborhood, int& startingEdge) const {

    startingEdge = 0;

    int tableSize = (int)_hashTable.size();

    key_type hashKey = neighborhood.GetHash();

    for (int psl = 0; psl < tableSize; ++psl) {

        address_type address = computeAddress(hashKey, psl, tableSize);

        Index planIndex = _hashTable[address];

        if (planIndex == INDEX_INVALID)
            return INDEX_INVALID;
       
        SubdivisionPlan const* plan = _plansTable[planIndex].get();
        if (plan && plan->IsTopologicalMatch(neighborhood, startingEdge))
            return planIndex;
    }
    return INDEX_INVALID;
}

inline Index 
TopologyMap::insertRegularFace(SubdivisionPlan* plan) {

    bool indexZero = kRegularPlanAtIndexZero ?
        plan->GetNeighborhood(0).IsRegularNeighborhood(_traits.getSchemeType()) : false;

    Index planIndex = appendPlan(plan, indexZero);

    insertHashKeys(planIndex);

    return planIndex;
}

inline Index 
TopologyMap::insertIrregularFace(SubdivisionPlan** plans, int numPlans) {

    assert(numPlans > 0);

    Index planIndex = appendPlan(plans[0], false);

    // only the neighborhoods of the plan for the first sub-face of an irregular
    // face need to be hashed; the rest is just appeneded in sequential (counter-
    // clockwise) order
    insertHashKeys(planIndex);

    for (int i = 1; i < numPlans; ++i)
        appendPlan(plans[i], false);

    return planIndex;
}

//
// public interface specializations (both lock-free and thread-safe)
//

constexpr bool lock_free = false;
constexpr bool thread_safe = true;

template <> Index
TopologyMap::FindSubdivisionPlan<lock_free>(
    Neighborhood const& neighborhood, int& startingEdge) const {
    return findPlan(neighborhood, startingEdge);
}
template <> Index
TopologyMap::FindSubdivisionPlan<thread_safe>(
    Neighborhood const& neighborhood, int& startingEdge) const {
    std::shared_lock lock(_mutex);
    return findPlan(neighborhood, startingEdge);
}


template <> Index
TopologyMap::InsertRegularFace<lock_free>(SubdivisionPlan* plan) {
    return insertRegularFace(plan);
}
template <> Index
TopologyMap::InsertRegularFace<thread_safe>(SubdivisionPlan* plan) {
    std::unique_lock lock(_mutex);
    return insertRegularFace(plan);
}


template <> Index
TopologyMap::InsertIrregularFace<lock_free>(SubdivisionPlan** plans, int numPlans) {
    return insertIrregularFace(plans, numPlans);
}
template <> Index
TopologyMap::InsertIrregularFace<thread_safe>(SubdivisionPlan** plans, int numPlans) {
    std::unique_lock lock(_mutex);
    return insertIrregularFace(plans, numPlans);
}

//
// statistics interface (non thread-safe)
//

size_t
TopologyMap::GetByteSize(bool includeNeighborhoods) const {
    size_t size = 0;
    for (int i = 0; i < GetNumSubdivisionPlans(); ++i)
        if (SubdivisionPlan const* plan = _plansTable[i].get())
            size += plan->GetByteSize(includeNeighborhoods);
    return size;
}

TopologyMap::HashmapStats
TopologyMap::ComputeHashTableStatistics() const {

    uint32_t const tableSize = (uint32_t)_hashTable.size();
    
    uint32_t pslMin = std::numeric_limits<uint32_t>::max();
    uint32_t pslMax = std::numeric_limits<uint32_t>::min();
    size_t pslSum = 0;

    // note: we ignore the regular surface plan if it is forced
    // at index 0, because PSL is always going to be 0 ; this
    // should not really make a significant difference either way
    // for most practical cases
    uint32_t firstPlan = kRegularPlanAtIndexZero;

    for (uint32_t planIndex = firstPlan; planIndex < (uint32_t)_plansTable.size(); ) {

        SubdivisionPlan const* plan = _plansTable[planIndex].get();

        for (int n = 0; n < plan->GetNumNeighborhoods(); ++n) {
            
            key_type hashKey = plan->GetNeighborhood(n).GetHash();

            for (uint32_t psl = 0; psl < tableSize; ++psl) {

                address_type address = computeAddress(hashKey, psl, tableSize);

                if (_hashTable[address] == planIndex) {
                    pslSum += psl;
                    pslMax = std::max(pslMax, psl);
                    pslMin = std::min(pslMin, psl);
                    break;
                }
            }
        }
        planIndex += plan->IsRegularFace() ? 1 : plan->GetFaceSize();
    }

    uint32_t numPlans = (uint32_t)_plansTable.size() - firstPlan;

    if (numPlans > 0)
        return {
            .pslMin = pslMin,
            .pslMax = pslMax,
            .pslMean = float(pslSum) / float(numPlans),
            .hashCount = numPlans,
            .addressCount = _hashTable.size(), 
            .loadFactor = loadFactor(), };
    else
        return {};
}

void
TopologyMap::WriteSubdivisionPlansDigraphs(
    FILE * fout, char const * title, bool showIndices) const {

    // write-lock here

    fprintf(fout, "digraph TopologyMap {\n");

    if (title) {
        fprintf(fout, "  label = \"%s\";\n", title);
        fprintf(fout, "  labelloc = \"t\";\n");
    }

    for (int planIndex = 0; planIndex < GetNumSubdivisionPlans(); ++planIndex) {

        if (SubdivisionPlan const* plan = GetSubdivisionPlan(planIndex))
            plan->WriteTreeDigraph(fout, planIndex, showIndices, /*isSubgraph*/ true);
    }
    fprintf(fout, "}\n");
}


} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv

