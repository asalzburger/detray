/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/coordinates/coordinates.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/materials/material_map.hpp"
#include "detray/materials/material_rod.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/utils/tuple_helpers.hpp"
#include "detray/utils/type_registry.hpp"
#include "detray/utils/type_traits.hpp"

// System include(s)
#include <type_traits>

namespace detray::io::detail {

/// Infer the IO shape id from the shape type
template <
    typename shape_t,
    std::enable_if_t<std::is_enum_v<typename shape_t::boundaries>, bool> = true>
constexpr io::shape_id get_id() {

    /// Register the mask shapes to the @c shape_id enum
    using shape_registry =
        type_registry<io::shape_id, annulus2D, cuboid3D, cylinder2D, cylinder3D,
                      concentric_cylinder2D, rectangle2D, ring2D, trapezoid2D,
                      wire_cell, straw_tube, single3D<0>, single3D<1>,
                      single3D<2>>;

    // Find the correct shape IO id;
    if constexpr (shape_registry::is_defined(shape_t{})) {
        return shape_registry::get_id(shape_t{});
    } else {
        return io::shape_id::unknown;
    }
}

/// Infer the IO material id from the material type - homogeneous material
template <typename material_t,
          std::enable_if_t<detray::detail::is_hom_material_v<material_t>,
                           bool> = true>
constexpr io::material_id get_id() {
    using scalar_t = typename material_t::scalar_type;

    /// Register the material types to the @c material_id enum
    using mat_registry =
        type_registry<io::material_id, void, void, void, void, void, void,
                      material_slab<scalar_t>, material_rod<scalar_t>,
                      material<scalar_t>>;

    // Find the correct material IO id;
    if constexpr (mat_registry::is_defined(material_t{})) {
        return mat_registry::get_id(material_t{});
    } else {
        return io::material_id::unknown;
    }
}

/// Infer the IO material id from the material type - material maps
template <typename material_t,
          std::enable_if_t<detray::detail::is_material_map_v<material_t> ||
                               detray::detail::is_volume_material_v<material_t>,
                           bool> = true>
constexpr io::material_id get_id() {

    using map_frame_t = typename material_t::local_frame_type;
    using algebra_t = typename map_frame_t::transform3_type;

    /// Register the material types to the @c material_id enum
    using mat_registry =
        type_registry<io::material_id, polar2<algebra_t>, cartesian2<algebra_t>,
                      cartesian3<algebra_t>, concentric_cylindrical2<algebra_t>,
                      cylindrical2<algebra_t>, cylindrical3<algebra_t>, void,
                      void>;

    // Find the correct material IO id;
    if constexpr (mat_registry::is_defined(map_frame_t{})) {
        return mat_registry::get_id(map_frame_t{});
    } else {
        return io::material_id::unknown;
    }
}

/// Infer the grid id from its coordinate system
template <
    typename grid_t,
    std::enable_if_t<detray::detail::is_surface_grid_v<grid_t>, bool> = true>
constexpr io::accel_id get_id() {

    using frame_t = typename grid_t::local_frame_type;
    using algebra_t = typename frame_t::transform3_type;

    /// Register the grid shapes to the @c accel_id enum
    /// @note the first type corresponds to a non-grid type in the enum
    /// (brute force)
    using frame_registry =
        type_registry<io::accel_id, void, cartesian2<algebra_t>,
                      cartesian3<algebra_t>, polar2<algebra_t>,
                      concentric_cylindrical2<algebra_t>,
                      cylindrical2<algebra_t>, cylindrical3<algebra_t>>;

    // Find the correct grid shape IO id;
    if constexpr (frame_registry::is_defined(frame_t{})) {
        return frame_registry::get_id(frame_t{});
    } else {
        return io::accel_id::unknown;
    }
}

/// Determine the type and id of a shape of a mask without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Mask shape unknown by detector
template <io::shape_id shape, typename detector_t, typename = void>
struct mask_info {
    using shape_id = typename detector_t::masks::id;
    using type = void;
    static constexpr shape_id value{detray::detail::invalid_value<shape_id>()};
};

/// Check for a stereo annulus shape
template <typename detector_t>
struct mask_info<io::shape_id::annulus2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<annulus2D, std::uint_least16_t>>(),
                                  void>> {
    using type = annulus2D;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_annulus2};
};

/// Check for a 2D cylinder shape
template <typename detector_t>
struct mask_info<io::shape_id::cylinder2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<cylinder2D, std::uint_least16_t>>(),
                                  void>> {
    using type = cylinder2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cylinder2};
};

/// Check for a 2D cylinder portal shape
template <typename detector_t>
struct mask_info<
    io::shape_id::portal_cylinder2, detector_t,
    std::enable_if_t<detector_t::masks::template is_defined<
                         mask<concentric_cylinder2D, std::uint_least16_t>>(),
                     void>> {
    using type = concentric_cylinder2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_portal_cylinder2};
};

/// Check for a cell wire line shape
template <typename detector_t>
struct mask_info<io::shape_id::cell_wire, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<wire_cell, std::uint_least16_t>>(),
                                  void>> {
    using type = wire_cell;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_cell_wire};
};

/// Check for a straw wire line shape
template <typename detector_t>
struct mask_info<io::shape_id::straw_wire, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<straw_tube, std::uint_least16_t>>(),
                                  void>> {
    using type = straw_tube;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_straw_wire};
};

/// Check for a rectangle shape
template <typename detector_t>
struct mask_info<io::shape_id::rectangle2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<rectangle2D, std::uint_least16_t>>(),
                                  void>> {
    using type = rectangle2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_rectangle2};
};

/// Check for a ring/disc shape
template <typename detector_t>
struct mask_info<io::shape_id::ring2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<ring2D, std::uint_least16_t>>(),
                                  void>> {
    using type = ring2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_portal_ring2};
};

/// Check for a single masked value (1st value is checked)
template <typename detector_t>
struct mask_info<io::shape_id::single1, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<0>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<0>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single1};
};

/// Check for a single masked value (2nd value is checked)
template <typename detector_t>
struct mask_info<io::shape_id::single2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<1>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<1>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single2};
};

/// Check for a single masked value (3rd value is checked)
template <typename detector_t>
struct mask_info<io::shape_id::single3, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<single3D<2>, std::uint_least16_t>>(),
                                  void>> {
    using type = single3D<2>;
    static constexpr
        typename detector_t::masks::id value{detector_t::masks::id::e_single3};
};

/// Check for a trapezoid shape
template <typename detector_t>
struct mask_info<io::shape_id::trapezoid2, detector_t,
                 std::enable_if_t<detector_t::masks::template is_defined<
                                      mask<trapezoid2D, std::uint_least16_t>>(),
                                  void>> {
    using type = trapezoid2D;
    static constexpr typename detector_t::masks::id value{
        detector_t::masks::id::e_trapezoid2};
};
/// @}

/// Determine the type and id of a material map without triggering a compiler
/// error (sfinae) if the detector does not know the type / enum entry
/// @{
/// Material map unknown by detector
template <io::material_id mat, typename detector_t, typename = void>
struct mat_map_info {
    using material_id = typename detector_t::materials::id;
    using type = void;
    static constexpr material_id value{
        detray::detail::invalid_value<material_id>()};
};

/// Check for a 2D disc material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::ring2_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         ring2D, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<ring2D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_disc2_map};
};

/// Check for a 2D cartesian material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::rectangle2_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         rectangle2D, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<rectangle2D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_rectangle2_map};
};

/// Check for a 3D cuboid volume material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::cuboid3_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         cuboid3D, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<cuboid3D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cuboid3_map};
};

/// Check for a 2D cylindrical material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::cylinder2_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         cylinder2D, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<cylinder2D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cylinder2_map};
};

/// Check for a 2D concentric cylindrical material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::concentric_cylinder2_map, detector_t,
    std::enable_if_t<
        detector_t::materials::template is_defined<material_map<
            concentric_cylinder2D, typename detector_t::scalar_type>>(),
        void>> {
    using type =
        material_map<concentric_cylinder2D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_concentric_cylinder2_map};
};

/// Check for a 3D cylindrical volume material map
template <typename detector_t>
struct mat_map_info<
    io::material_id::cylinder3_map, detector_t,
    std::enable_if_t<detector_t::materials::template is_defined<material_map<
                         cylinder3D, typename detector_t::scalar_type>>(),
                     void>> {
    using type = material_map<cylinder3D, typename detector_t::scalar_type>;
    static constexpr typename detector_t::materials::id value{
        detector_t::materials::id::e_cylinder3_map};
};
/// @}

}  // namespace detray::io::detail
