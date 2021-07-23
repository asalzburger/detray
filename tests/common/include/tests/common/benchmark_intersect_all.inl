/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/detector.hpp"
#include "core/transform_store.hpp"
#include "io/csv_io.hpp"
#include "tools/intersection_kernel.hpp"

#include <iostream>
#include <fstream>

#include <benchmark/benchmark.h>

using namespace detray;

using transform3 = __plugin::transform3;
using point3 = point3;
using vector3 = vector3;
using surface = surface_base<transform3>;

__plugin::cartesian2 cartesian2;
using point2 = __plugin::point2;

#ifdef DETRAY_BENCHMARKS_REP
unsigned int gbench_repetitions = DETRAY_BENCHMARKS_REP;
#else
unsigned int gbench_repetitions = 0;
#endif

unsigned int theta_steps = 100;
unsigned int phi_steps = 100;
bool stream_file = false;

/** Read the detector from file */
auto read_detector()
{
    auto env_d_d = std::getenv("DETRAY_TEST_DATA_DIR");
    if (env_d_d == nullptr)
    {
        throw std::ios_base::failure("Test data directory not found. Please set DETRAY_TEST_DATA_DIR.");
    }
    auto data_directory = std::string(env_d_d);

    std::string name = "tml";
    std::string surfaces = data_directory + "odd.csv";
    std::string volumes = data_directory + "odd-layer-volumes.csv";
    std::string grids = data_directory + "odd-surface-grids.csv";
    std::string grid_entries = "";
    return detray::detector_from_csv<>(name, surfaces, volumes, grids, grid_entries);
};

auto d = read_detector();

namespace __plugin
{
    // This test runs intersection with all surfaces of the TrackML detector
    static void BM_INTERSECT_ALL(benchmark::State &state)
    {

        /*std::ofstream hit_out;
        if (stream_file)
        {
            hit_out.open("tml_hits.csv");
        }*/
        unsigned int hits = 0;
        unsigned int missed = 0;

        point3 ori = {0., 0., 0.};

        using detray_context = decltype(d)::transform_store::context;
        detray_context default_context;

        for (auto _ : state)
        {
            track<detray_context> track;
            track.pos = point3{0., 0., 0.};

            // Loops of theta values
            for (unsigned int itheta = 0; itheta < theta_steps; ++itheta)
            {
                scalar theta = 0.1 + itheta * (M_PI - 0.1) / theta_steps;
                scalar sin_theta = std::sin(theta);
                scalar cos_theta = std::cos(theta);

                // Loops of phi values
                for (unsigned int iphi = 0; iphi < phi_steps; ++iphi)
                {
                    // The direction
                    scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
                    scalar sin_phi = std::sin(phi);
                    scalar cos_phi = std::cos(phi);
                    track.dir = {cos_phi * sin_theta, sin_phi * sin_theta, cos_theta};

                    // Loop over volumes
                    for (const auto &v : d.volumes())
                    {
                        const auto &surfaces = v.surfaces();

                        // Loop over surfaces
                        for (const auto &s : surfaces.objects())
                        {
                            auto sfi_surface = intersect(track, s, surfaces.transforms(default_context), surfaces.masks());

                            const auto &sfi = std::get<0>(sfi_surface);
                            
                            benchmark::DoNotOptimize(hits);
                            benchmark::DoNotOptimize(missed);
                            if (sfi.status == intersection_status::e_inside)
                            {
                                /* state.PauseTiming();
                                if (stream_file)
                                {
                                    hit_out << sfi.p3[0] << "," << sfi.p3[1] << "," << sfi.p3[2] << "\n";
                                }
                                state.ResumeTiming();*/
                                ++hits;
                            }
                            else
                            {
                                ++missed;
                            }
                            benchmark::ClobberMemory();
                        }
                    }
                }
            }
        }

        #ifndef DETRAY_BENCHMARKS_MULTITHREAD
        std::cout << "[detray] hits / missed / total = " << hits << " / " << missed << " / " << hits + missed << std::endl;
        #endif
        /**if (stream_file)
        {
            hit_out.close();
        }*/
    }

    BENCHMARK(BM_INTERSECT_ALL)
    #ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadPerCpu()
    #endif
    ->Unit(benchmark::kMillisecond)
    ->Repetitions(gbench_repetitions)
    ->DisplayAggregatesOnly(true);

} // namespace __plugin

BENCHMARK_MAIN();
