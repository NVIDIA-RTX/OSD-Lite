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

#ifndef OPENSUBDIV3_TMR_SURFACE_TABLE_H
#define OPENSUBDIV3_TMR_SURFACE_TABLE_H

#include "../version.h"

#include "../tmr/surfaceDescriptor.h"
#include "../tmr/topologyMap.h"

#include <vector>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Tmr {

class SubdivisionPlan;

///
/// \brief Base interface for tables of 'Surfaces'
/// 
/// note: 'irregular' sized faces require multiple surfaces, so the number
/// of 'Surface' entries in the table is often greater than the number of faces in
/// the base mesh (the indexing of Surface entries is matching Ptex face indices).
/// 

template<typename Descriptor> struct SurfaceTableInterface {

    int GetNumSurfaces() const { return (int)descriptors.size(); }

    Descriptor GetDescriptor(Index surfaceIndex) const { return descriptors[surfaceIndex]; }

    Index const* GetControlPointIndices(Index surfaceIndex) const;

    size_t GetByteSize() const;

    // note: 'regular' mesh faces are associated a single SurfaceDescriptor, while 
    // 'irregular' mesh faces are associated 'n' Surfaces, where 'n' is the number
    // of edges of that  face (ex. a Catmark triangle requires 3 Surfaces). 
    // Ordering of sub-faces is counter-clockwise and matches Ptex face indexing.
    std::vector<Descriptor> descriptors;

    // 'controlPoints' is a dense array of the control point indices for the 1-ring
    // of each base mesh face (the array of control points is not duplicated for the
    // sub-faces of irregular faces)
    std::vector<Index> controlPointIndices;
};


struct SurfaceTable : public SurfaceTableInterface<SurfaceDescriptor> {

    TopologyMap const& topologyMap;

    Domain GetDomain(Index surfaceIndex) const;

    SubdivisionPlan const& GetSubdivisionPlan(Index surfaceIndex) const;

    SurfaceTable(TopologyMap const& topologyMap) : topologyMap(topologyMap) { }

    template <typename T, typename U> int
        EvaluatePatchPoints(Index surfaceIndex, T const* controlPoints, U* patchPoints, int level = kMaxIsolationLevel) const;
};

struct LinearSurfaceTable : public SurfaceTableInterface<LinearSurfaceDescriptor> {

    int numControlPointsMax = 0;

    Domain GetDomain(Index surfaceIndex) const;

    // evaluate only the patch-points relevant to the surface at 'surfaceIndex'
    template <typename T, typename U> int
        EvaluateLocalPatchPoints(Index surfaceIndex, T const* controlPoints, U* patchPoints) const;
    
    // evaluate all the patch-points of the parent face if surfaceIndex points
    // to an irregular face
    template <typename T, typename U> int
        EvaluatePatchPoints(Index surfaceIndex, T const* controlPoints, U* patchPoints) const;
};

//
// Implementation
//

template <typename T, typename U> int
LinearSurfaceTable::EvaluateLocalPatchPoints(Index surfaceIndex, T const* controlPoints, U* patchPoints) const {

    LinearSurfaceDescriptor desc = GetDescriptor(surfaceIndex);

    LocalIndex subFace = desc.GetQuadSubfaceIndex();

    if (subFace == LOCAL_INDEX_INVALID)
        return 0;

    Index const* indices = GetControlPointIndices(surfaceIndex);

    int N = desc.GetFaceSize();

    float invN = 1.f / float(N);

    U& facePoint = patchPoints[0];
    facePoint.Clear();
    for (int i = 0; i < N; ++i) {
        facePoint.AddWithWeight(controlPoints[indices[i]], invN);
    }

    LocalIndex subFaceNext = subFace == (N - 1) ? 0 : subFace + 1;
    patchPoints[1].Set(controlPoints[indices[subFace]], .5f);
    patchPoints[1].AddWithWeight(controlPoints[indices[subFaceNext]], .5f);

    LocalIndex subFacePrev = (subFace > 0 ? subFace : N) - 1 ;
    patchPoints[2].Set(controlPoints[indices[subFacePrev]], .5f);
    patchPoints[2].AddWithWeight(controlPoints[indices[subFace]], .5f);

    return 3;
}

template <typename T, typename U> int
LinearSurfaceTable::EvaluatePatchPoints(Index surfaceIndex, T const* controlPoints, U* patchPoints) const {

    LinearSurfaceDescriptor desc = GetDescriptor(surfaceIndex);

    if (desc.GetQuadSubfaceIndex() == LOCAL_INDEX_INVALID)
        return 0;

    Index const* indices = GetControlPointIndices(surfaceIndex);

    int N = desc.GetFaceSize();
    float invN = 1.f / float(N);

    U& facePoint = patchPoints[0];

    facePoint.Clear();

    for (int i = 0; i < N; ++i) {

        facePoint.AddWithWeight(controlPoints[indices[i]], invN);

        int j = (i < (N - 1)) ? i + 1 : 0;

        U& edgePoint = patchPoints[i + 1];
        edgePoint.Set(controlPoints[indices[i]], .5f);
        edgePoint.AddWithWeight(controlPoints[indices[j]], .5f);
    }
    return N + 1;
}

template<typename T, typename U> inline int
SurfaceTable::EvaluatePatchPoints(Index surfaceIndex, T const* controlPoints, U* patchPoints, int level) const {

    SurfaceDescriptor desc = GetDescriptor(surfaceIndex);

    SubdivisionPlan const& plan = *topologyMap.GetSubdivisionPlan(desc.GetSubdivisionPlanIndex());

    ConstIndexArray controlPointIndices = { GetControlPointIndices(surfaceIndex), plan.GetNumControlPoints() };

    return plan.EvaluatePatchPoints<T, U>(controlPoints, controlPointIndices, patchPoints, level);
}

} // end namespace Tmr

} // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;
} // end namespace OpenSubdiv

#endif /* OPENSUBDIV3_TMR_SURFACE_TABLE_H */

