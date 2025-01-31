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


#ifndef OPENSUBDIV3_TMR_TYPES_H
#define OPENSUBDIV3_TMR_TYPES_H


#include "../version.h"

#include "../vtr/types.h"

#include <bit>
#include <cstdint>
#include <limits>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

// Endian-agnostic bit-field packing & unpacking

constexpr uint32_t pack(uint32_t value, int width, int offset) {
    return (uint32_t)((value & ((1 << width) - 1)) << offset);
}

constexpr uint32_t unpack(uint32_t value, int width, int offset) {
    return (uint32_t)((value >> offset) & ((1 << width) - 1));
}

constexpr uint8_t countbits(uint32_t value) {
    return std::popcount(value);
}

//
//  Typedefs for indices that are inherited from the Vtr level -- eventually
//  these primitive Vtr types may be declared at a lower, more public level.
//

typedef Vtr::Index       Index;
typedef Vtr::LocalIndex  LocalIndex;

typedef Vtr::IndexArray       IndexArray;
typedef Vtr::ConstIndexArray  ConstIndexArray;

typedef Vtr::LocalIndexArray       LocalIndexArray;
typedef Vtr::ConstLocalIndexArray  ConstLocalIndexArray;

typedef Vtr::Array<float>       FloatArray;
typedef Vtr::ConstArray<float>  ConstFloatArray;


inline bool IndexIsValid(Index index) { return Vtr::IndexIsValid(index); }

static const Index INDEX_INVALID = Vtr::INDEX_INVALID;
static const LocalIndex LOCAL_INDEX_INVALID = ~Vtr::LocalIndex(0);

// maximum adaptive isolation level of a Surface
static constexpr uint8_t const kMaxIsolationLevel = 10u;

// maximum number of SubdivisionPlan in a TopologyMap ; constrained by
// SurfaceDescriptor 'plan index' bit field (24 bits)
static constexpr uint32_t const kMaxNumSubdivisionPlans = (1u << 24) - 1u;

// constrained by Neighborhood::_startingEdge (uint16_t)
static constexpr int const VALENCE_LIMIT = std::numeric_limits<uint16_t>::max();

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif // OPENSUBDIV3_TMR_TYPES_H
