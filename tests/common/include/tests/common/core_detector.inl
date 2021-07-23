/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"

#include <gtest/gtest.h>

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, detector)
{
    using namespace detray;
    using namespace __plugin;

    using detector = detector<>;

    static_transform_store<>::storage static_storage;
    static_transform_store<>::context ctx0;

    detector::surface_mask_container masks;

    /// Surface 0
    point3 t0{0., 0., 0.};
    transform3 tf0{t0};
    static_storage.push_back(std::move(tf0));
    detector::surface_rectangle rect = {-3., 3.};
    std::get<detector::surface_rectangle::mask_context>(masks).push_back(rect);

    /// Surface 1
    point3 t1{1., 0., 0.};
    transform3 tf1{t1};
    static_storage.push_back(std::move(tf1));
    detector::surface_annulus anns = {1., 2., 3., 4., 5., 6., 7.};
    std::get<detector::surface_annulus::mask_context>(masks).push_back(anns);

    /// Surface 2
    point3 t2{2., 0., 0.};
    transform3 tf2{t2};
    static_storage.push_back(std::move(tf2));
    detector::surface_trapezoid trap = {1., 2., 3.};
    std::get<detector::surface_trapezoid::mask_context>(masks).push_back(trap);

    detector d("test_detector");
    auto &v = d.new_volume("test_volume", {0., 10., -5., 5., -M_PI, M_PI});
    d.add_surface_transforms(ctx0, v, std::move(static_storage));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
