/** Detray library, part of the ACTS project
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/coordinates/coordinate_base.hpp"
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Frame projection into a cartesian coordinate frame
 */
template <typename transform3_t>
struct cartesian2 final : public coordinate_base<cartesian2, transform3_t> {

    /// @name Type definitions for the struct
    /// @{

    /// Base type
    using base_type = coordinate_base<cartesian2, transform3_t>;
    /// Sclar type
    using scalar_type = typename base_type::scalar_type;
    /// Point in 2D space
    using point2 = typename base_type::point2;
    /// Point in 3D space
    using point3 = typename base_type::point3;
    /// Vector in 3D space
    using vector3 = typename base_type::vector3;
    /// Matrix actor
    using matrix_operator = typename base_type::matrix_operator;
    /// Matrix size type
    using size_type = typename base_type::size_type;
    /// 2D matrix type
    template <size_type ROWS, size_type COLS>
    using matrix_type = typename base_type::template matrix_type<ROWS, COLS>;

    // Trigonometrics
    using trigonometrics = typename base_type::trigonometrics;
    // Vector types
    using bound_vector = typename base_type::bound_vector;
    using free_vector = typename base_type::free_vector;

    /// @}

    /** This method transform from a point from 2D cartesian frame to a 2D
     * cartesian point */
    DETRAY_HOST_DEVICE
    inline point2 operator()(const point2 &local2) const { return local2; }

    /** This method transform from a point from 3D cartesian frame to a 2D
     * cartesian point */
    DETRAY_HOST_DEVICE
    inline point2 operator()(const point3 &local3) const {

        return {local3[0], local3[1]};
    }

    /** This method transform from a point from global cartesian 3D frame to a
     * local 2D cartesian point */
    DETRAY_HOST_DEVICE
    inline point2 global_to_local(const transform3_t &trf, const point3 &p,
                                  const vector3 & /*d*/) const {
        const auto local3 = trf.point_to_local(p);
        return this->operator()(local3);
    }

    /** This method transform from a local 2D cartesian point to a point global
     * cartesian 3D frame*/
    template <typename mask_t>
    DETRAY_HOST_DEVICE inline point3 local_to_global(
        const transform3_t &trf, const mask_t & /*mask*/, const point2 &p,
        const vector3 & /*d*/) const {
        return trf.point_to_global(point3{p[0], p[1], 0.});
    }

    DETRAY_HOST_DEVICE
    inline matrix_type<3, 2> bound_to_free_rotation(
        const transform3_t &trf3, const trigonometrics & /*t*/) {

        // Get d(x,y,z)/d(loc0, loc1)
        return matrix_operator().template block<3, 2>(trf3.matrix(), 0, 0);
    }

    DETRAY_HOST_DEVICE
    inline matrix_type<2, 3> free_to_bound_rotation(
        const transform3_t &trf3, const trigonometrics & /*t*/) {

        // Get transpose of transform3 matrix
        const auto trf3T = matrix_operator().transpose(trf3);

        // Get d(loc0, loc1)/d(x,y,z)
        return matrix_operator().template block<2, 3>(trf3T.matrix(), 0, 0);
    }

};  // struct cartesian2

}  // namespace detray