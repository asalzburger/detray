
/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "grids/associator.hpp"
#include "utils/enumerate.hpp"
#include "utils/generators.hpp"

#include <array>
#include <vector>
#include <tuple>

namespace detray
{
    using point2 = __plugin::point2;
    using point3 = __plugin::point3;

    /** Run the bin association of surfaces (via their contour)
     *  - to a given grid.
     * 
     * @param context is the context to win which the association is done
     * @param volume is the volume to which the volume belongs
     * @param grid is the grid that will be filled 
     * @param tolerance is the bin_tolerance in the two local coordinates
     * @param absolute_tolerance is an indicator if the tolerance is to be
     *        taken absolute or relative
     */
    template <typename context_type,
              typename volume_type,
              typename grid_type>
    static inline void bin_association(const context_type &context,
                                       const volume_type &volume,
                                       grid_type &grid,
                                       const std::array<scalar, 2> &bin_tolerance,
                                       bool absolute_tolerance = true)
    {

        // Get surfaces, transforms and masks
        const auto &surfaces = volume.surfaces();
        const auto &surface_transforms = surfaces.transforms(context);
        const auto &surface_masks = surfaces.masks();

        const auto &bounds = volume.bounds();
        bool is_cylinder = std::abs(bounds[1] - bounds[0]) < std::abs(bounds[3] - bounds[2]);

        auto sfm_copy = surface_masks;

        const auto &axis_0 = grid.axis_p0();
        const auto &axis_1 = grid.axis_p1();

        // Disk type bin association
        if (not is_cylinder)
        {
            // Run with two different associators: center of gravity and edge intersection
            center_of_gravity_generic cgs_assoc;
            edges_intersect_generic edges_assoc;

            // Loop over all bins and associate the surfaces
            for (unsigned int bin_0 = 0; bin_0 < axis_0.bins(); ++bin_0)
            {
                for (unsigned int bin_1 = 0; bin_1 < axis_1.bins(); ++bin_1)
                {

                    auto r_borders = axis_0.borders(bin_0);
                    auto phi_borders = axis_1.borders(bin_1);

                    scalar r_add = absolute_tolerance ? bin_tolerance[0] : bin_tolerance[0] * (r_borders[1] - r_borders[0]);
                    scalar phi_add = absolute_tolerance ? bin_tolerance[1] : bin_tolerance[1] * (phi_borders[1] - phi_borders[0]);

                    // Create a contour for the bin
                    std::vector<point2> bin_contour = r_phi_polygon(r_borders[0] - r_add,
                                                                    r_borders[1] + r_add,
                                                                    phi_borders[0] - phi_add,
                                                                    phi_borders[1] + phi_add);

                    // Run through the surfaces and associate them by contour
                    for (auto [is, s] : enumerate(surfaces.objects()))
                    {
                        dvector<point3> vertices = {};
                        const auto &mask_link = s.mask();

                        // Unroll the mask container and generate vertices
                        const auto &transform = surface_transforms[s.transform()];

                        const auto &mask_context = std::get<0>(mask_link);
                        const auto &mask_range = std::get<1>(mask_link);

                        auto vertices_per_masks = unroll_masks_for_vertices(surface_masks, mask_range, mask_context,
                                                                            std::make_integer_sequence<dindex, std::tuple_size_v<decltype(sfm_copy)>>{});

                        // Usually one mask per surface, but design allows - a single association  is sufficient though
                        for (auto &vertices : vertices_per_masks)
                        {
                            if (not vertices.empty())
                            {
                                // Create a surface contour
                                std::vector<point2> surface_contour;
                                surface_contour.reserve(vertices.size());
                                for (const auto &v : vertices)
                                {
                                    auto vg = transform.point_to_global(v);
                                    surface_contour.push_back({vg[0], vg[1]});
                                }
                                // The association has worked
                                if (cgs_assoc(bin_contour, surface_contour) or edges_assoc(bin_contour, surface_contour))
                                {
                                    dindex bin_index = is;
                                    grid.populate(bin_0, bin_1, std::move(bin_index));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {

            center_of_gravity_rectangle cgs_assoc;
            edges_intersect_generic edges_assoc;

            // Loop over all bins and associate the surfaces
            for (unsigned int bin_0 = 0; bin_0 < axis_0.bins(); ++bin_0)
            {
                for (unsigned int bin_1 = 0; bin_1 < axis_1.bins(); ++bin_1)
                {

                    auto z_borders = axis_0.borders(bin_0);
                    auto phi_borders = axis_1.borders(bin_1);

                    scalar z_add = absolute_tolerance ? bin_tolerance[0] : bin_tolerance[0] * (z_borders[1] - z_borders[0]);
                    scalar phi_add = absolute_tolerance ? bin_tolerance[1] : bin_tolerance[1] * (phi_borders[1] - phi_borders[0]);

                    scalar z_min = z_borders[0];
                    scalar z_max = z_borders[1];
                    scalar phi_min = phi_borders[0];
                    scalar phi_max = phi_borders[1];

                    point2 p0_bin = {z_min - z_add, phi_min - phi_add};
                    point2 p1_bin = {z_min - z_add, phi_max + phi_add};
                    point2 p2_bin = {z_max + z_add, phi_max + phi_add};
                    point2 p3_bin = {z_max + z_add, phi_min - phi_add};

                    std::vector<point2> bin_contour = {p0_bin, p1_bin, p2_bin, p3_bin};

                    // Loop over the surfaces within a volume
                    for (auto [is, s] : enumerate(surfaces.objects()))
                    {
                        dvector<point3> vertices = {};
                        const auto &mask_link = s.mask();

                        // Unroll the mask container and generate vertices
                        const auto &transform = surface_transforms[s.transform()];

                        const auto &mask_context = std::get<0>(mask_link);
                        const auto &mask_range = std::get<1>(mask_link);

                        auto vertices_per_masks = unroll_masks_for_vertices(surface_masks, mask_range, mask_context,
                                                                            std::make_integer_sequence<dindex, std::tuple_size_v<decltype(sfm_copy)>>{});

                        for (auto &vertices : vertices_per_masks)
                        {

                            if (not vertices.empty())
                            {
                                // Create a surface contour
                                std::vector<point2> surface_contour;
                                surface_contour.reserve(vertices.size());
                                scalar phi_min = std::numeric_limits<scalar>::max();
                                scalar phi_max = -std::numeric_limits<scalar>::max();
                                // We poentially need the split vertices
                                std::vector<point2> s_c_neg;
                                std::vector<point2> s_c_pos;
                                scalar z_min_neg = std::numeric_limits<scalar>::max();
                                scalar z_max_neg = -std::numeric_limits<scalar>::max();
                                scalar z_min_pos = std::numeric_limits<scalar>::max();
                                scalar z_max_pos = -std::numeric_limits<scalar>::max();

                                for (const auto &v : vertices)
                                {
                                    const point3 vg = transform.point_to_global(v);
                                    scalar phi = std::atan2(vg[1], vg[0]);
                                    phi_min = std::min(phi, phi_min);
                                    phi_max = std::max(phi, phi_max);
                                    surface_contour.push_back({vg[2], phi});
                                    if (phi < 0.)
                                    {
                                        s_c_neg.push_back({vg[2], phi});
                                        z_min_neg = std::min(vg[2], z_min_neg);
                                        z_max_neg = std::max(vg[2], z_max_neg);
                                    }
                                    else
                                    {
                                        s_c_pos.push_back({vg[2], phi});
                                        z_min_pos = std::min(vg[2], z_min_pos);
                                        z_max_pos = std::max(vg[2], z_max_pos);
                                    }
                                }
                                // Check for phi wrapping
                                std::vector<std::vector<point2>> surface_contours;
                                if (phi_max - phi_min > M_PI and phi_max * phi_min < 0.)
                                {
                                    s_c_neg.push_back({z_max_neg, -M_PI});
                                    s_c_neg.push_back({z_min_neg, -M_PI});
                                    s_c_pos.push_back({z_max_pos, M_PI});
                                    s_c_pos.push_back({z_min_pos, M_PI});
                                    surface_contours = {s_c_neg, s_c_pos};
                                }
                                else
                                {
                                    surface_contours = {surface_contour};
                                }

                                // Check the association (with potential splits)
                                bool associated = false;
                                for (const auto &s_c : surface_contours)
                                {
                                    if (cgs_assoc(bin_contour, s_c) or edges_assoc(bin_contour, s_c))
                                    {
                                        associated = true;
                                        break;
                                    }
                                }

                                // Register if associated
                                if (associated)
                                {
                                    dindex bin_index = is;
                                    grid.populate(bin_0, bin_1, std::move(bin_index));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

} // namespace detray
