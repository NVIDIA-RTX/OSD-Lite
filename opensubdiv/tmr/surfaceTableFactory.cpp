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


#include "../tmr/surfaceTableFactory.h"
#include "../tmr/neighborhoodBuilder.h"
#include "../tmr/subdivisionPlanBuilder.h"
#include "../tmr/topologyMap.h"
#include "../tmr/types.h"

#include "../far/topologyRefiner.h"
#include "../vtr/fvarLevel.h"
#include "../vtr/level.h"

#include <limits>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

using Far::TopologyRefiner;
using Vtr::internal::Level;
using Vtr::internal::FVarLevel;

Sdc::Options fromTraits(TopologyMap::Traits const& traits) {
    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(traits.getVtxBoundaryInterpolation());
    options.SetCreasingMethod(traits.getCreasingMethod());
    options.SetTriangleSubdivision(traits.getTriangleSubdivision());   
    return options;
}

static unsigned int 
countSurfaces(Level const& level, int regFaceSize) {
    // note : holes aren't skipped - they are populated with 'no limit' descriptors (see below)
    unsigned int count = 0;
    for (int face = 0; face < level.getNumFaces(); ++face) {
        ConstIndexArray fverts = level.getFaceVertices(face);
        count += fverts.size() == regFaceSize ? 1 : fverts.size();
    }      
    return count;
}

static uint8_t
computeEdgeAdjacencyBits(Level const& level, Index faceIndex, int regFaceSize) {   
    // the quadrangulation of irregular faces introduces T-junctions that tessellation algorithms
    // need to stitch into water-tight meshes.
    // the edge adjcency bit is set to true if one or more faces adjacent to that edge is irregular
    // note: expects faceIndex to point to a regular face
    assert(level.isFaceHole(faceIndex) == false);

    ConstIndexArray edges = level.getFaceEdges(faceIndex);
    assert(edges.size() == regFaceSize);

    uint8_t edgeAdjacency = 0;
    for (int edge = 0; edge < edges.size(); ++edge) {       
        
        ConstIndexArray edgeFaces = level.getEdgeFaces(edges[edge]);
        
        for (int face = 0; face < edgeFaces.size(); ++face) {

            Index adjFaceIndex = edgeFaces[face];

            if (adjFaceIndex == faceIndex || level.isFaceHole(adjFaceIndex))
                continue;

            if (level.getFaceVertices(adjFaceIndex).size() != regFaceSize) {
                edgeAdjacency |= (1u << edge);
                break;
            }
        }
    }
    return edgeAdjacency;
}

template <typename descriptor_type> inline int
setNoLimit(std::vector<descriptor_type>& descriptors, Index surfIndex, int n, int regFaceSize) {
    // holes & other faces that do not render still get a set of 'no limit' descriptors to maintain
    // primtive ordering ; same indexing as ptex, where irregular sized faces require 'n' descriptors
    if (n == regFaceSize) {
        descriptors[surfIndex].SetNoLimit();
        return 1;
    }
    for (int i = 0; i < n; ++i)
        descriptors[surfIndex + i].SetNoLimit();
    return n;
}

uint8_t
computeFvarRot(SurfaceDescriptor vtxDesc, int startingEdge, int regFaceSize) {
    // subdivision plans generate hash keys for each possible neighborhood rotation ; when a surface
    // is matched to a given plan, the implicit parametric space of the basis function will also be
    // rotated.
    // for non-linear face-varying data, there are 2 successive rotations applied: one issued from hashing
    // the vertex neighborhood, and a counter-rotation issued from hashing the face-varying neighborhood. 
    // both must be un-done to retrieve the correct tangent space.
    int rotation = vtxDesc.GetParametricRotation() - startingEdge;
    if (rotation < 0)
        rotation += regFaceSize;
    assert(rotation < 4);
    return rotation;
}


template <bool thread_safe> std::unique_ptr<SurfaceTable>
SurfaceTableFactory::Create(TopologyRefiner const& refiner, 
    TopologyMap& topologyMap, SurfaceTableFactory::Options const& options) {

    TopologyMap::Traits traits = topologyMap.GetTraits();
    
    Index const topologyMapID = topologyMap.GetOptions().uniqueID;

    assert(traits.getSchemeType() == refiner.GetSchemeType());
    assert(traits.getEndCapType() == options.planBuilderOptions.endCapType);
    assert(traits.getCreasingMethod() == refiner.GetSchemeOptions().GetCreasingMethod());
    assert(traits.getTriangleSubdivision() == refiner.GetSchemeOptions().GetTriangleSubdivision());

    Sdc::SchemeType schemeType = traits.getSchemeType();
    Sdc::Options schemeOptions = fromTraits(traits);

    int const regFaceSize = Sdc::SchemeTypeTraits::GetRegularFaceSize(schemeType);

    bool boundaryNone = schemeOptions.GetVtxBoundaryInterpolation() == Sdc::Options::VTX_BOUNDARY_NONE;

    Level const& level = refiner.getLevel(0);
    FVarLevel const* fvlevel = options.fvarChannel >= 0 ? &level.getFVarLevel(options.fvarChannel) : nullptr;      

    int numFaces = level.getNumFaces();
    int numSurfaces = countSurfaces(level, regFaceSize);
    int startingEdge = 0;

    Vtr::internal::StackBuffer<uint8_t, 512> neighborhoodData;

    NeighborhoodBuilder::BuildDescriptor buildDesc = { .refiner = refiner, .fvarChannel = options.fvarChannel, };

    auto populateNeighborhoodData = [&](int faceIndex, int startingEdge = 0) -> Neighborhood const* {
        buildDesc.faceIndex = faceIndex;
        buildDesc.startingEdge = startingEdge;
        return _neighborhoodBuilder.Populate(neighborhoodData, buildDesc) > 0 ?
            reinterpret_cast<Neighborhood*>(&neighborhoodData[0]) : nullptr;
    };
    auto createNeighborhoodData = [&](int faceIndex, int startingEdge = 0) -> std::unique_ptr<uint8_t const[]> {
        buildDesc.faceIndex = faceIndex;
        buildDesc.startingEdge = startingEdge;
        return _neighborhoodBuilder.Create(buildDesc);
    };

    std::unique_ptr<SurfaceTable> surfaceTable = std::make_unique<SurfaceTable>(topologyMap);
    surfaceTable->descriptors.resize(numSurfaces);
    surfaceTable->controlPointIndices.reserve(numFaces * 20);

    for (Index faceIndex = 0, surfIndex = 0; faceIndex < numFaces; ++faceIndex) {

        int n = level.getNumFaceVertices(faceIndex);

        if (level.isFaceHole(faceIndex)) {
            setNoLimit(surfaceTable->descriptors, surfIndex, n, regFaceSize);
            surfIndex += n == regFaceSize ? 1 : n;
            continue;
        }
        
        Neighborhood const* neighborhood = populateNeighborhoodData(faceIndex);

        if (!neighborhood) {
            // degenerate face of some kind that somehow was not detected
            setNoLimit(surfaceTable->descriptors, surfIndex, n, regFaceSize);
            surfIndex += n == regFaceSize ? 1 : n;
            continue;
        }

        assert(neighborhood->GetBaseFaceSize() == n);

        Index planIndex = topologyMap.FindSubdivisionPlan<thread_safe>(*neighborhood, startingEdge);

        if (planIndex == INDEX_INVALID) {

            auto buildPlan = [&](LocalIndex subfaceIndex = 0) -> SubdivisionPlan* {

                std::unique_ptr<SubdivisionPlan> plan = _planBuilder.Create(
                    schemeType, schemeOptions, options.planBuilderOptions, *neighborhood, subfaceIndex);

                plan->reserveNeighborhoods(n);

                for (int subfaceIndex = 0; subfaceIndex < n; ++subfaceIndex) {
                    auto neighborhoodData = createNeighborhoodData(faceIndex, subfaceIndex);
                    plan->addNeighborhood(std::move(neighborhoodData), subfaceIndex);
                }
                return plan.release();
            };

            if (n == regFaceSize) {
                
                planIndex = topologyMap.InsertRegularFace<thread_safe>(buildPlan());

            } else {

                Vtr::internal::StackBuffer<SubdivisionPlan*, 16> plans(n);

                for (int subfaceIndex = 0; subfaceIndex < n; ++subfaceIndex)
                    plans[subfaceIndex] = buildPlan(subfaceIndex);

                planIndex = topologyMap.InsertIrregularFace<thread_safe>(&plans[0], n);
            }
        } else {

            if (startingEdge > 0) {
                // rotated variant of the SubdivisionPlan: rebuild the 1-ring of control points, but with
                // ordering starting from the matching vertex/edge
                //_neighborhoodBuilder.Populate(neighborhoodData, level, faceIndex, startingEdge);      
                neighborhood = populateNeighborhoodData(faceIndex, startingEdge);
            }
        }

        assert(neighborhood->GetNumControlPoints() == topologyMap.GetSubdivisionPlan(planIndex)->GetNumControlPoints());

        // populate the SurfaceTable
 
        unsigned int firstControlPoint = (unsigned int)surfaceTable->controlPointIndices.size();

        ConstIndexArray controlPoints = neighborhood->GetControlPoints(); 
        for (int i = 0; i < controlPoints.size(); ++i)
            surfaceTable->controlPointIndices.push_back(controlPoints[i]);

        if (n == regFaceSize) {
            uint8_t edgeAdjacency = computeEdgeAdjacencyBits(level, faceIndex, regFaceSize);
            uint8_t paramRotation = fvlevel && options.depTable ?
                computeFvarRot(options.depTable->GetDescriptor(surfIndex), startingEdge, regFaceSize) : startingEdge;
            surfaceTable->descriptors[surfIndex].Set(firstControlPoint, planIndex, paramRotation, edgeAdjacency, topologyMapID);
        } else {
            // note: the sub-faces of an irregular face share the same set of control points (the 
            // 'firstControlPoint' offset is the same for all the descriptors)
            uint8_t edgeAdjacency = 0;
            uint8_t paramRotation = 0;
            for (int i = 0; i < n; ++i) {
                // if the surface has been rotated (startingEdge > 0), the sub-faces need to be re-ordered:
                // simply shuffle the order of descriptors
                int descIndex = (i + startingEdge) % n;
                surfaceTable->descriptors[surfIndex + descIndex].Set(firstControlPoint, planIndex + i, paramRotation, edgeAdjacency, topologyMapID);
            }
        }

        surfIndex += n == regFaceSize ? 1 : n;
    }

    surfaceTable->controlPointIndices.shrink_to_fit();

    return surfaceTable;
}

template std::unique_ptr<SurfaceTable>
SurfaceTableFactory::Create<false>(TopologyRefiner const& refiner,
    TopologyMap& topologyMap, SurfaceTableFactory::Options const& options);

template std::unique_ptr<SurfaceTable>
SurfaceTableFactory::Create<true>(TopologyRefiner const& refiner,
    TopologyMap& topologyMap, SurfaceTableFactory::Options const& options);


std::unique_ptr<LinearSurfaceTable>
LinearSurfaceTableFactory::Create(TopologyRefiner const& refiner, int fvarChannel, SurfaceTable const* depTable) {

    std::unique_ptr<LinearSurfaceTable> surfaceTable = std::make_unique<LinearSurfaceTable>();

    Level const& level = refiner.getLevel(0);

    int const regFaceSize = Sdc::SchemeTypeTraits::GetRegularFaceSize(refiner.GetSchemeType());

    int numFaces = level.getNumFaces();
    int numSurfaces = countSurfaces(level, regFaceSize);

    surfaceTable->descriptors.resize(numSurfaces);

    for (Index faceIndex = 0, surfIndex = 0, firstControlPoint = 0; faceIndex < numFaces; ++faceIndex) {

        // face-vertex and face-value arrays have the same layout
        int n = level.getNumFaceVertices(faceIndex);

        if (level.isFaceHole(faceIndex)) {
            surfIndex += setNoLimit(surfaceTable->descriptors, surfIndex, n, regFaceSize);
        } else {
            if (n == regFaceSize) {               
                surfaceTable->descriptors[surfIndex].Set(firstControlPoint, n);
                ++surfIndex;
            } else {
                for (int i = 0 ; i < n; ++i)
                    surfaceTable->descriptors[surfIndex + i].Set(firstControlPoint, n, i);
                surfIndex += n;
            }
            surfaceTable->numControlPointsMax = std::max(surfaceTable->numControlPointsMax, n);
        }
        firstControlPoint += n;
    }

    // populate the control points array

    if (fvarChannel >= 0 && depTable) {

        assert(depTable->GetNumSurfaces() == numSurfaces);

        // the parametric space of a surface is potentially rotated if it has been matched to
        // a rotated version of an existing subdivision plan in the topology map. Bi-linear
        // face-varying interpolation does not use subdivision plans and only require the 0-ring
        // control points: we can fetch the rotation of the subdivision from the dependent 
        // vertex-interpolation surface table & apply it directly here. This makes these
        // rotations completely transparent (no cost) to the run-time evaluation.

        FVarLevel const& fvlevel = level.getFVarLevel(fvarChannel);
        
        surfaceTable->controlPointIndices.resize(fvlevel.getNumFaceValuesTotal());

        Index* controlPoints = surfaceTable->controlPointIndices.data();

        for (Index faceIndex = 0, surfIndex = 0; faceIndex < numFaces; ++faceIndex) {

            ConstIndexArray values = fvlevel.getFaceValues(faceIndex);        
            SurfaceDescriptor desc = depTable->GetDescriptor(surfIndex);
            
            int n = values.size();

            if (int rotation = desc.GetParametricRotation()) {
                for (int i = 0; i < n; ++i)
                    controlPoints[i] = values[(i + rotation) % n];
            } else
                std::memcpy(controlPoints, values.begin(), n * sizeof(Index));

            surfIndex += n == regFaceSize ? 1 : n;
            controlPoints += n;
        }
    } else {
        // general case (not face-varying or no dependent vertex surface table) : copy
        // the face-varying values from the refiner. Sub-faces from faces of irregular
        // size share the same set of control points.
        ConstIndexArray controlPoints = fvarChannel == -1 ? 
            level.getFaceVertices() : level.getFVarLevel(fvarChannel).getFaceValues();
        
        surfaceTable->controlPointIndices.assign(controlPoints.begin(), controlPoints.begin() + controlPoints.size());
    }

    surfaceTable->controlPointIndices.shrink_to_fit();

    return surfaceTable;
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
} // end namespace OpenSubdiv
