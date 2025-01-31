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

#ifndef OPENSUBDIV3_TMR_SURFACE_TABLE_FACTORY_H
#define OPENSUBDIV3_TMR_SURFACE_TABLE_FACTORY_H

#include "../version.h"

#include "../tmr/surfaceTable.h"
#include "../tmr/neighborhoodBuilder.h"
#include "../tmr/subdivisionPlanBuilder.h"
#include "../tmr/topologyMap.h"

#include <memory>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Far {
    class TopologyRefiner;
}

namespace Tmr {

/// \brief Hashes mesh topology from a TopologyRefiner to build tables of per-face 'Surfaces'

class LinearSurfaceTableFactory {

public:

    // Builds a LinearSurfaceTable for the given control mesh topology.
    // note: bi-linear evaluation does not require a topology map
    std::unique_ptr<LinearSurfaceTable> Create(Far::TopologyRefiner const& refiner,
        int fvarChannel = -1, SurfaceTable const* depTable = nullptr);
};

class SurfaceTableFactory {

public:

    struct Options {

        SubdivisionPlanBuilder::Options planBuilderOptions;

        // note: in order to generate a SurfaceTable for cubic face-varying data,
        // users must specify both the face-varying channel index from the refiner,
        // and the dependent SurfaceTable generated for vertex data interpolation.
        int fvarChannel = -1;

        SurfaceTable const* depTable = nullptr;
    };

    // Builds a (cubic) SurfaceTable for the given control mesh topology
    // The topology map is expected to be shared across many meshes for memory re-use
    template <bool thread_safe = false> std::unique_ptr<SurfaceTable> Create(
        Far::TopologyRefiner const& refiner, TopologyMap& topologyMap, Options const& options);

private:

    NeighborhoodBuilder _neighborhoodBuilder;

    SubdivisionPlanBuilder _planBuilder;
};

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_SURFACE_TABLE_FACTORY_H */

