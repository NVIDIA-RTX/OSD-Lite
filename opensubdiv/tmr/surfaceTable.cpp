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

#include "../tmr/surfaceTable.h"
#include "../tmr/neighborhood.h"
#include "../tmr/subdivisionPlan.h"
#include "../tmr/topologyMap.h"

#include <cassert>
#include <cmath>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

template<typename Descriptor> Index const*
    SurfaceTableInterface<Descriptor>::GetControlPointIndices(Index surfaceIndex) const {
    return controlPointIndices.data() + descriptors[surfaceIndex].firstControlPoint;
}
template Index const* SurfaceTableInterface<SurfaceDescriptor>::GetControlPointIndices(Index surfaceIndex) const;
template Index const* SurfaceTableInterface<LinearSurfaceDescriptor>::GetControlPointIndices(Index surfaceIndex) const;


template<typename Descriptor> size_t
    SurfaceTableInterface<Descriptor>::GetByteSize() const {
    size_t size = 0;
    size += sizeof(SurfaceTableInterface);
    size += descriptors.size() * sizeof(typename decltype(descriptors)::value_type);
    size += controlPointIndices.size() * sizeof(typename decltype(controlPointIndices)::value_type);
    return size;
}
template size_t SurfaceTableInterface<SurfaceDescriptor>::GetByteSize() const;
template size_t SurfaceTableInterface<LinearSurfaceDescriptor>::GetByteSize() const;


Domain
LinearSurfaceTable::GetDomain(Index surfaceIndex) const {   
    return GetDescriptor(surfaceIndex).GetDomain();
}

Domain
SurfaceTable::GetDomain(Index surfaceIndex) const {   

    switch (topologyMap.GetTraits().getSchemeType()) {
        case Sdc::SCHEME_CATMARK: {
            SubdivisionPlan const& plan = GetSubdivisionPlan(surfaceIndex);
            return plan.IsRegularFace() ? Domain::Quad : Domain::Quad_Subface;
        }
        case Sdc::SCHEME_LOOP: {
            return Domain::Tri;
        }
        default:
            break;
    }
    assert(0);
    return Domain::Quad;
}

SubdivisionPlan const&
SurfaceTable::GetSubdivisionPlan(Index surfaceIndex) const {
    SurfaceDescriptor desc = GetDescriptor(surfaceIndex);
    // It is safe to dereference the plan pointer because the surface descriptor
    // will never point to a plan that does not exist.
    return *topologyMap.GetSubdivisionPlan(desc.GetSubdivisionPlanIndex()); 
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
