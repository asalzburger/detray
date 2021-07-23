/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/track.hpp"
#include "core/intersection.hpp"
#include "utils/indexing.hpp"
#include "utils/enumerate.hpp"

#include <utility>
#include <tuple>

namespace detray
{

    /// Transform definition
    using transform3 = __plugin::transform3;

    /** Specialized method to update an intersection when the mask group is resolved
     * 
     * @tparam track_type is the type of the track information including context 
     * @tparam mask_group is the type of the split out mask group from all masks
     * @tparam mask_range is the type of the according mask range object
     * 
     * @param track is the track information 
     * @param trf is the surface transform for the intersection
     * @param masks is the associated (and split out) mask group
     * @param range is the range list of masks to be processed
     * 
     * @note since all masks are on the same surface, only maximally one mask solution at a time
     * is possbile. Once this solution is found, it is returned
     * 
     * @return the surface intersection 
     **/
    template <typename track_type,
              typename mask_group,
              typename mask_range>
    auto intersect_by_group(const track_type &track,
                            const transform3 &trf,
                            const mask_group &masks,
                            const mask_range &range)
    {
        for (auto i : sequence(range))
        {
            const auto &mask = masks[i];
            auto local = mask.local();
            auto sfi = mask.intersector().intersect(trf, track, local, mask);
            if (sfi.status == e_inside)
            {
                auto valid_link = mask.links();
                return std::make_tuple(sfi, valid_link);
            }
        }
        typename mask_group::value_type::mask_links_type invalid_link;
        intersection invalid_intersection;
        return std::make_tuple(invalid_intersection, invalid_link);
    }

    /** Variadic unrolled intersection - last entry 
     * 
     * @tparam intersection_type The type of intersection (to be updated)
     * @tparam links_type The links type (to be updated)
     * @tparam track_type is The type of the surface container
     * @tparam mask_container is the type of the type of the mask container
     * @tparam mask_range is the mask range type 
     * @tparam last_mask_context is the last mask group context 
     * 
     * @param intersection [in, out] the intersection to be updated
     * @param link [in, out] the links to be updated
     * @param tack_type the track information containing the context
     * @param ctf the contextual transform (context resolved)
     * @param masks the masks container 
     * @param range the range within the mask group to be checked
     * @param mask_context the last mask group context 
     *
     * @return a bool for break condition (ignored here as last)
     */
    template <typename intersection_type,
              typename links_type,
              typename track_type,
              typename mask_container,
              typename mask_range,
              dindex last_mask_context>
    bool last_intersect(intersection_type &intersection,
                        links_type &links,
                        const track_type &track,
                        const transform3 &ctf,
                        const mask_container &masks,
                        const mask_range &range,
                        dindex mask_context)
    {
        if (mask_context == last_mask_context)
        {
            auto isg = intersect_by_group(track, ctf, std::get<last_mask_context>(masks), range);
            intersection = std::get<0>(isg);
            links = std::get<1>(isg);
            return true;
        }
        return false;
    }

    /** Variadic unrolled intersection - any integer sequence
     * 
     * @tparam intersection_type The type of intersection (to be updated)
     * @tparam links_type The links type (to be updated)
     * @tparam track_type is The type of the surface container
     * @tparam mask_container is the type of the type of the mask container
     * @tparam mask_range is the mask range type 
     * @tparam first_mask_context is the first mask group context 
     * 
     * @param intersection [in, out] the intersection to be updated
     * @param link [in, out] the links to be updated
     * @param tack_type the track information containing the context
     * @param ctf the contextual transform (context resolved)
     * @param masks the masks container 
     * @param range the range within the mask group to be checked
     * @param mask_context the last mask group context 
     * @param available_contices the mask contices to be checked
     *
     * @return a bool for break condition 
     */
    template <typename intersection_type,
              typename links_type,
              typename track_type,
              typename mask_container,
              typename mask_range,
              dindex first_mask_context,
              dindex... remaining_mask_context>
    bool unroll_intersect(intersection_type &intersection,
                          links_type &links,
                          const track_type &track,
                          const transform3 &ctf,
                          const mask_container &masks,
                          const mask_range &range,
                          dindex mask_context,
                          std::integer_sequence<dindex, first_mask_context, remaining_mask_context...> available_contices)
    {
        // Pick the first one for interseciton
        if (mask_context == first_mask_context)
        {
            auto isg = intersect_by_group(track, ctf, std::get<first_mask_context>(masks), range);
            intersection = std::get<0>(isg);
            links = std::get<1>(isg);
            return true;
        }
        // The reduced integer sequence
        std::integer_sequence<dindex, remaining_mask_context...> remaining;
        // Unroll as long as you have at least 2 entries
        if constexpr (remaining.size() > 1)
        {
            if (unroll_intersect(intersection, links, track, ctf, masks, range, mask_context, remaining))
            {
                return true;
            }
        }
        // Last chance - intersect the last index if possible
        return last_intersect<intersection_type,
                              links_type,
                              track_type,
                              mask_container,
                              mask_range,
                              std::tuple_size_v<mask_container> - 1>(intersection, links, track, ctf, masks, range, mask_context);
    }

    /** Actual kernel method that updates the updates the intersections
     * 
     * @tparam surface_type is The type of the surface container
     * @tparam transform_container is the type of the transform container
     * @tparam mask_container is the type of the type of the mask container
     *
     * @param track the track information including the contexts
     * @param surface the surface type to be intersected
     * @param contextual_transform the transform container
     * @param mask_container the tuple mask container to for the intersection
     * 
     * @return an intersection and the link to the result
    **/
    template <typename surface_type,
              typename transform_container,
              typename mask_container>
    const auto
    intersect(const track<typename transform_container::context> &track,
              const surface_type &surface,
              const transform_container &contextual_transforms,
              const mask_container &masks)
    {
        const auto &ctf = contextual_transforms[surface.transform()];
        auto mask_link = surface.mask();
        const auto &mask_context = std::get<0>(mask_link);
        const auto &mask_range = std::get<1>(mask_link);

        // Create a return intersection and run the variadic unrolling
        const auto &reference_group = std::get<0>(masks);
        typename std::decay_t<decltype(reference_group)>::value_type::mask_links_type result_links;
        intersection result_intersection;

        // Unroll the intersection depending on the mask container size
        unroll_intersect(result_intersection, result_links, track, ctf, masks, mask_range, mask_context,
                         std::make_integer_sequence<dindex, std::tuple_size_v<mask_container>>{});
        // Return the (eventually update) intersection and links
        return std::make_tuple(result_intersection, result_links);
    }
}
