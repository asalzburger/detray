/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "core/surface_base.hpp"
#include "core/track.hpp"
#include "masks/unmasked.hpp"
#include "utils/unbound.hpp"

#include <cmath>
#include <climits>
#include <optional>

namespace detray
{

    /** This is an intersector struct for planar surface
     */
    struct planar_intersector
    {

        using transform3 = __plugin::transform3;
        using point3 = __plugin::point3;
        using vector3 = __plugin::vector3;
        using point2 = __plugin::point2;

        /** Intersection method for planar surfaces
         * 
         * @tparam track_type The type of the track (which carries the context object)
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param trf the surface to be intersected
         * @param track the track information
         * @param local to the local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * @param tolerance is the mask specific tolerance
         * 
         * @return the intersection with optional parameters
         **/
        template <typename track_type, typename local_type, typename mask_type>
        intersection
        intersect(const transform3 &trf,
                  const track_type &track,
                  const local_type &local,
                  const mask_type &mask,
                  const typename mask_type::mask_tolerance &tolerance = mask_type::within_epsilon) const
        {
            return intersect(trf, track.pos, track.dir, local, mask, tolerance, track.overstep_tolerance);
        }

        /** Intersection method for planar surfaces
         * 
         * @tparam local_type The local frame type to be intersected
         * @tparam mask_type The mask type applied to the local frame
         * 
         * Contextual part:
         * @param trf the transform of the surface to be intersected
         * @param ro the origin of the ray
         * @param rd the direction of the ray
         * @param local transform to the local local frame 
         * 
         * Non-contextual part:
         * @param mask the local mask 
         * @param tolerance is the mask specific tolerance
         * 
         * @return the intersection with optional parameters
         **/
        template <typename local_type = unbound, typename mask_type = unmasked>
        intersection
        intersect(const transform3 &trf,
                  const point3 &ro,
                  const vector3 &rd,
                  const local_type &local = local_type(),
                  const mask_type &mask = mask_type(),
                  const typename mask_type::mask_tolerance &tolerance = mask_type::within_epsilon,
                  scalar overstep_tolerance = 0.) const
        {

            // Retrieve the surface normal & translation (context resolved)
            const auto &sm = trf.matrix();
            auto sn = getter::vector<3>(sm, 0, 2);
            auto st = getter::vector<3>(sm, 0, 3);

            // Intersection code
            scalar denom = vector::dot(rd, sn);
            if (denom != 0.0)
            {
                intersection is;
                is.path = vector::dot(sn, (st - ro)) / (denom);
                is.p3 = ro + is.path * rd;
                is.p2 = local(trf, is.p3);
                is.status = mask.template is_inside<local_type>(is.p2, tolerance);
                is.direction = is.path > overstep_tolerance ? e_along : e_opposite;
                return is;
            }
            return intersection{};
        }
    };

} // namespace detray
