/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection.hpp"

// System include(s)
#include <cmath>
#include <type_traits>

namespace detray {

/// @brief Intersection implementation for cylinder surfaces using helical
/// trajectories.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask. On the @c cylinder3
/// mask, it switches on the check of the radial distance.
template <typename transform3_t>
struct helix_cylinder_intersector {

    using scalar_type = typename transform3_t::scalar_type;
    using matrix_operator = typename transform3_t::matrix_actor;
    using point3 = typename transform3_t::point3;
    using vector3 = typename transform3_t::vector3;
    using helix_type = detail::helix<transform3_t>;

    using intersection_type = line_plane_intersection;
    using output_type = std::array<intersection_type, 2>;

    /// Operator function to find intersections between helix and cylinder mask
    ///
    /// @tparam mask_t is the input mask type
    /// @tparam transform_t is the input transform type
    ///
    /// @param h is the input helix trajectory
    /// @param mask is the input mask
    /// @param trf is the transform
    /// @param mask_tolerance is the tolerance for mask edges
    /// @param overstep_tolerance is the tolerance for track overstepping
    ///
    /// @return the intersection
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const helix_type &h, const mask_t &mask, const transform3_t &trf,
        const scalar_type mask_tolerance = 0) const {

        output_type ret;

        using local_frame = typename mask_t::local_type;

        // Guard against inifinite loops
        constexpr std::size_t max_n_tries{100};
        // Tolerance for convergence
        constexpr scalar_type tol{1e-3};

        // Get the surface placement
        const auto &sm = trf.matrix();
        // Cylinder z axis
        const vector3 sz = getter::vector<3>(sm, 0, 2);
        // Cylinder centre
        const point3 sc = getter::vector<3>(sm, 0, 3);

        // Starting point on the helix for the Newton iteration
        // The mask is a cylinder -> it provides its radius as the first value
        const scalar_type r{mask[0]};
        // Helix path length parameter
        scalar_type s{r * getter::perp(h.dir(tol))};
        // Path length in the previous iteration step
        scalar_type s_prev{s - scalar{0.1}};

        // f(s) = ((h.pos(s) - sc) x sz)^2 - r^2 == 0
        // Run the iteration on s
        std::size_t n_tries{0};
        while (std::abs(s - s_prev) > tol and n_tries < max_n_tries) {

            // f'(s) = 2 * ( (h.pos(s) - sc) x sz) * (h.dir(s) x sz) )
            const vector3 crp = vector::cross(h.pos(s) - sc, sz);
            const scalar_type denom{
                scalar_type{2.} *
                vector::dot(crp, vector::cross(h.dir(s), sz))};
            // No intersection can be found if dividing by zero
            if (denom == scalar_type{0.}) {
                return ret;
            }
            // x_n+1 = x_n - f(s) / f'(s)
            s_prev = s;
            s -= (vector::dot(crp, crp) - r * r) / denom;

            ++n_tries;
        }
        // No intersection found within max number of trials
        if (n_tries == max_n_tries) {
            return ret;
        }

        // Build intersection struct from helix parameter s
        intersection_type &is = ret[0];
        const point3 helix_pos = h.pos(s);

        is.path = getter::norm(helix_pos);
        is.p3 = helix_pos;
        constexpr local_frame local_converter{};
        is.p2 = local_converter.global_to_local(trf, is.p3, h.dir(s));

        auto local3 = trf.point_to_local(is.p3);
        // Explicitly check for radial match
        is.status =
            mask.template is_inside<local_frame, true>(local3, mask_tolerance);
        is.direction = vector::dot(is.p3, h.dir(s)) > scalar_type{0.}
                           ? intersection::direction::e_along
                           : intersection::direction::e_opposite;
        is.link = mask.volume_link();

        return ret;
    }
};

}  // namespace detray