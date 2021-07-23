/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "core/detector.hpp"
#include "core/volume_connector.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "io/csv_io_types.hpp"
#include "tools/bin_association.hpp"

#include <climits>
#include <map>
#include <vector>

namespace detray
{

  /// Function to read the detector from the CSV file
  ///
  /// @tparam array_type is the type of the array used for the detector class
  /// @tparam tuple_type is the type of the tuple used for the detector class
  /// @tparam vector_type is the type of the tuple used for the detector class
  ///
  /// All other container types are STL containers, as they do not leave the function
  ///
  /// @param detector_name is the name of the detector
  /// @param surface_file_name is the name of the detector surface file
  /// @param grid_file_name is the name of the surface grid file
  /// @param layer_volume_file_name is the name of the file containing layer/volume information
  /// @param r_sync_tolerance is a tolerance to be for synching volumes in r
  /// @param z_sync_tolerance is a toleranced to be used for synchinng volumes in z
  ///
  /// @return a detector object
  template <template <typename, unsigned int> class array_type = darray,
            template <typename...> class tuple_type = dtuple,
            template <typename> class vector_type = dvector,
            typename alignable_store = static_transform_store<vecm_vector>,
            typename surface_source_link = dindex,
            typename bounds_source_link = dindex>
  detector<array_type,
           tuple_type,
           vector_type,
           alignable_store,
           surface_source_link,
           bounds_source_link>
  detector_from_csv(const std::string &detector_name,
                    const std::string &surface_file_name,
                    const std::string &layer_volume_file_name,
                    const std::string &grid_file_name,
                    const std::string &grid_entries_file_name,
                    scalar r_sync_tolerance = 0.,
                    scalar z_sync_tolerance = 0.)
  {
    using typed_detector = detector<array_type, tuple_type, vector_type, alignable_store, surface_source_link, bounds_source_link>;

    typed_detector d(detector_name);

    // Surface reading
    surface_reader s_reader(surface_file_name);
    csv_surface io_surface;

    // Layer volume reading
    layer_volume_reader lv_reader(layer_volume_file_name);
    csv_layer_volume io_layer_volume;

    // Surface grid reading
    surface_grid_reader sg_reader(grid_file_name);
    csv_surface_grid io_surface_grid;

    using volume_layer_index = std::pair<uint32_t, uint32_t>;
    std::map<volume_layer_index, typename typed_detector::volume *> volumes;

    // Read in with a default context
    typename alignable_store::storage surface_transform_storage;
    typename alignable_store::context surface_default_context;

    // Flushable containers
    typename typed_detector::volume *c_volume = nullptr;
    typename typed_detector::surface_container c_surfaces;
    typename typed_detector::surface_mask_container c_masks;

    std::map<volume_layer_index, array_type<scalar, 6>> volume_bounds;

    // Remember the r/z attachements
    std::map<scalar, std::vector<dindex>> r_min_attachments;
    std::map<scalar, std::vector<dindex>> z_min_attachments;
    scalar r_max = 0;
    scalar z_max = -std::numeric_limits<scalar>::max();

    /** Helper method to attach volumes to bins
     * 
     * @param attachments The attachnment map
     * @param value the (eventually new) value for insertion
     * @param volume_index the index of the attached volume
     */
    auto attach_volume = [](std::map<scalar, std::vector<dindex>> &attachments, scalar value, dindex volume_index) -> void
    {
      if (attachments.find(value) == attachments.end())
      {
        attachments[value] = {volume_index};
      }
      else
      {
        attachments[value].push_back(volume_index);
      }
    };

    // (A) Pre-read the bounds values
    std::map<scalar, scalar> z_min_layer_volumes;
    std::map<scalar, scalar> z_max_layer_volumes;
    std::map<scalar, scalar> r_min_layer_volumes;
    std::map<scalar, scalar> r_max_layer_volumes;

    while (lv_reader.read(io_layer_volume))
    {
      volume_layer_index c_index = {io_layer_volume.volume_id, io_layer_volume.layer_id};
      array_type<scalar, 6> c_bounds = {
          io_layer_volume.min_v0,
          io_layer_volume.max_v0,
          io_layer_volume.min_v1,
          io_layer_volume.max_v1,
          io_layer_volume.min_v2,
          io_layer_volume.max_v2,
      };
      volume_bounds[c_index] = c_bounds;

      // Register low, high volume
      if (io_layer_volume.layer_id % 2 == 0)
      {
        z_min_layer_volumes[io_layer_volume.min_v1] = io_layer_volume.min_v1;
        z_max_layer_volumes[-1. * io_layer_volume.max_v1] = io_layer_volume.max_v1;
        r_min_layer_volumes[io_layer_volume.min_v0] = io_layer_volume.min_v0;
        r_max_layer_volumes[-1. * io_layer_volume.max_v0] = io_layer_volume.max_v0;
      }
    }
    /** Helper function to cluster boundaries 
     * 
     * @param boundaries is the unclustered boundaries map
     * @param tolerance is the tolerance parameter in which a cluster can lie
     * @param flip is a reverse flag for upper boundaries
     * 
     */
    auto cluster_boundaries = [&](std::map<scalar, scalar> &boundaries, scalar tolerance, int flip = 1) -> void
    {
      scalar last_ = std::numeric_limits<scalar>::max();
      for (auto &[key, boundary] : boundaries)
      {
        // Do not adust the last one for max values
        if (flip < 0 and key == boundaries.begin()->first)
        {
          continue;
        }

        if (std::abs(last_ - flip * key) < tolerance)
        {
          boundary = last_;
        }
        last_ = boundary;
      }
    };

    // Cluster boundary synchronization
    cluster_boundaries(z_min_layer_volumes, 5.);
    cluster_boundaries(z_max_layer_volumes, 5., -1);
    cluster_boundaries(r_min_layer_volumes, 5.);
    cluster_boundaries(r_max_layer_volumes, 5., -1);

    /** Helper functions to find and replace in case
     * 
     * @param value is the value in question
     * @param value_map is the map for the replacement value
     * 
     **/
    auto find_and_replace = [](scalar &value, const std::map<scalar, scalar> &value_map, int flip = 1) -> void
    {
      // synchronize lower bound
      if (value_map.find(flip * value) != value_map.end())
      {
        value = value_map.find(flip * value)->second;
      }
    };

    /** Helper function to return syncrhonized boundary objects
     * 
     * @param bounds the unsynchronized bounds object
     * @param gap_volume an indicator if its a gap volume 
     **/
    auto synchronize_bounds = [&](const array_type<scalar, 6> &bounds, bool gap_volume) -> array_type<scalar, 6>
    {
      scalar r_min = bounds[0];
      scalar r_max = bounds[1];
      scalar z_min = bounds[2];
      scalar z_max = bounds[3];
      scalar phi_min = bounds[4];
      scalar phi_max = bounds[5];

      if (not gap_volume)
      {
        find_and_replace(r_min, r_min_layer_volumes);
        find_and_replace(r_max, r_max_layer_volumes, -1);
        find_and_replace(z_min, z_min_layer_volumes);
        find_and_replace(z_max, z_max_layer_volumes, -1);
      }
      else
      {
        if (std::abs(z_max + z_min) < 0.1)
        {
          find_and_replace(r_min, r_max_layer_volumes, -1);
          find_and_replace(r_max, r_min_layer_volumes);
          find_and_replace(z_min, z_min_layer_volumes);
          find_and_replace(z_max, z_max_layer_volumes, -1);
        }
        else
        {
          find_and_replace(r_min, r_min_layer_volumes);
          find_and_replace(r_max, r_max_layer_volumes, -1);
          find_and_replace(z_min, z_max_layer_volumes, -1);
          find_and_replace(z_max, z_min_layer_volumes);
        }
      }

      return {r_min, r_max, z_min, z_max, phi_min, phi_max};
    };

    // Create the surface finders & reserve
    std::map<volume_layer_index, dindex> surface_finder_entries;

    using surfaces_r_axis = typename typed_detector::surfaces_regular_axis;
    using surfaces_z_axis = typename typed_detector::surfaces_regular_axis;
    using surfaces_phi_axis = typename typed_detector::surfaces_circular_axis;

    using surfaces_r_phi_grid = typename typed_detector::surfaces_regular_circular_grid;
    using surfaces_z_phi_grid = typename typed_detector::surfaces_regular_circular_grid;
    using surfaces_finder = typename typed_detector::surfaces_finder;

    vector_type<surfaces_finder> detector_surfaces_finders;

    // (B) Pre-read the grids & create local object finders
    while (sg_reader.read(io_surface_grid))
    {
      volume_layer_index c_index = {io_surface_grid.volume_id, io_surface_grid.layer_id};

      surface_finder_entries[c_index] = detector_surfaces_finders.size();

      bool is_disk = (io_surface_grid.type_loc0 == 3);

      // Prepare z axis parameters
      scalar z_min = is_disk ? std::numeric_limits<scalar>::min() : io_surface_grid.min_loc1;
      scalar z_max = is_disk ? std::numeric_limits<scalar>::max() : io_surface_grid.max_loc1;
      dindex z_bins = is_disk ? 1u : static_cast<dindex>(io_surface_grid.nbins_loc1);

      // Prepare r axis parameters
      scalar r_min = is_disk ? io_surface_grid.min_loc0 : 0.;
      scalar r_max = is_disk ? io_surface_grid.max_loc0 : std::numeric_limits<scalar>::max();
      dindex r_bins = is_disk ? static_cast<dindex>(io_surface_grid.nbins_loc0) : 1u;

      // Prepare phi axis parameters
      scalar phi_min = is_disk ? io_surface_grid.min_loc1 : io_surface_grid.min_loc0;
      scalar phi_max = is_disk ? io_surface_grid.max_loc1 : io_surface_grid.max_loc0;
      dindex phi_bins = is_disk ? static_cast<dindex>(io_surface_grid.nbins_loc1) : static_cast<dindex>(io_surface_grid.nbins_loc0);

      surfaces_z_axis z_axis{z_bins, z_min, z_max};
      surfaces_r_axis r_axis{r_bins, r_min, r_max};
      surfaces_phi_axis phi_axis{phi_bins, phi_min, phi_max};

      // negative / positive / inner / outer
      surfaces_r_phi_grid rphi_grid_n(r_axis, phi_axis);
      surfaces_r_phi_grid rphi_grid_p(r_axis, phi_axis);
      surfaces_z_phi_grid zphi_grid_i{z_axis, phi_axis};
      surfaces_z_phi_grid zphi_grid_o{z_axis, phi_axis};

      detector_surfaces_finders.push_back(rphi_grid_n);
      detector_surfaces_finders.push_back(rphi_grid_p);
      detector_surfaces_finders.push_back(zphi_grid_i);
      detector_surfaces_finders.push_back(zphi_grid_o);
    }

    // (C) Read the surfaces and fill it
    while (s_reader.read(io_surface))
    {
      volume_layer_index c_index = {io_surface.volume_id, io_surface.layer_id};
      auto c_volume_itr = volumes.find(c_index);
      if (c_volume_itr == volumes.end())
      {
        // Flush the former information / c_volume still points to the prior volume
        if (c_volume != nullptr and not surface_transform_storage.empty())
        {
          // Construction with move semantics
          //c_volume->add_surface_transforms(surface_default_context, std::move(surface_transform_storage));
          d.add_surface_transforms(surface_default_context, *c_volume, std::move(surface_transform_storage));
          c_volume->add_surface_components(std::move(c_surfaces), std::move(c_masks));

          // Get new clean containers
          surface_transform_storage = typename alignable_store::storage();
          c_surfaces = typename typed_detector::surface_container();
          c_masks = typename typed_detector::surface_mask_container();
        }

        // Create a new volume & assign
        std::string volume_name = detector_name;
        volume_name += std::string("_vol_") + std::to_string(io_surface.volume_id);
        volume_name += std::string("_lay_") + std::to_string(io_surface.layer_id);
        // Find and fill the bounds
        auto new_bounds = volume_bounds.find(c_index);
        if (new_bounds == volume_bounds.end())
        {
          // Bounds not found, do not build the volume
          continue;
        }

        const auto &unsynchronized_volume_bounds = new_bounds->second;
        // Check if you need to synchronize
        bool is_gap = (io_surface.layer_id % 2 != 0);
        auto volume_bounds = synchronize_bounds(unsynchronized_volume_bounds, is_gap);

        // Check if this volume has a surface finder entry associated
        dindex surfaces_finder_entry = dindex_invalid;
        auto surface_finder_itr = surface_finder_entries.find(c_index);
        if (surface_finder_itr != surface_finder_entries.end())
        {
          surfaces_finder_entry = surface_finder_itr->second;
        }

        auto &new_volume = d.new_volume(volume_name, volume_bounds, surfaces_finder_entry);

        // RZ attachment storage
        attach_volume(r_min_attachments, volume_bounds[0], new_volume.index());
        attach_volume(z_min_attachments, volume_bounds[2], new_volume.index());

        r_max = std::max(r_max, volume_bounds[1]);
        z_max = std::max(z_max, volume_bounds[3]);

        c_volume = &new_volume;
        // Insert to volume map
        volumes[c_index] = c_volume;
      }
      else
      {
        c_volume = c_volume_itr->second;
      }

      // Do not fill navigation layers
      if (io_surface.layer_id % 2 == 0)
      {

        // Read the transform
        vector3 t{io_surface.cx, io_surface.cy, io_surface.cz};
        vector3 x{io_surface.rot_xu, io_surface.rot_yu, io_surface.rot_zu};
        vector3 z{io_surface.rot_xw, io_surface.rot_yw, io_surface.rot_zw};

        // Translate the mask & add it to the mask container
        unsigned int bounds_type = io_surface.bounds_type;
        std::vector<scalar> bounds;
        bounds.push_back(io_surface.bound_param0);
        bounds.push_back(io_surface.bound_param1);
        bounds.push_back(io_surface.bound_param2);
        bounds.push_back(io_surface.bound_param3);
        bounds.push_back(io_surface.bound_param4);
        bounds.push_back(io_surface.bound_param5);
        bounds.push_back(io_surface.bound_param6);

        // Acts naming convention for bounds
        typename typed_detector::surface_mask_index mask_index = {dindex_invalid, dindex_invalid};

        if (bounds_type == 1)
        {
          // Cylinder bounds
          constexpr auto cylinder_context = typed_detector::surface_cylinder::mask_context;

          // Get the cylinder mask container
          auto &cylinder_masks = std::get<cylinder_context>(c_masks);
          dindex cylinder_index = cylinder_masks.size();
          cylinder_masks.push_back({io_surface.bound_param0, io_surface.cz - io_surface.bound_param1, io_surface.cz + io_surface.bound_param1});
          // The read is valid: set the index
          mask_index = {cylinder_context, cylinder_index};
        }
        else if (bounds_type == 3)
        {
          // Disc bounds
        }
        else if (bounds_type == 6)
        {
          // Rectangle bounds
          constexpr auto rectangle_context = typed_detector::surface_rectangle::mask_context;

          // Get the rectangle mask container
          auto &rectangle_masks = std::get<rectangle_context>(c_masks);
          dindex rectangle_index = rectangle_masks.size();
          scalar half_x = 0.5 * (io_surface.bound_param2 - io_surface.bound_param0);
          scalar half_y = 0.5 * (io_surface.bound_param3 - io_surface.bound_param1);
          rectangle_masks.push_back({half_x, half_y});
          // The read is valid: set the index
          mask_index = {rectangle_context, rectangle_index};
        }
        else if (bounds_type == 7)
        {
          // Trapezoid bounds
          constexpr auto trapezoid_context = typed_detector::surface_trapezoid::mask_context;
          // Get the trapezoid mask container
          auto &trapezoid_masks = std::get<trapezoid_context>(c_masks);
          dindex trapezoid_index = trapezoid_masks.size();
          trapezoid_masks.push_back({io_surface.bound_param0, io_surface.bound_param1, io_surface.bound_param2});
          // The read is valid: set the index
          mask_index = {trapezoid_context, trapezoid_index};
        }
        else if (bounds_type == 11)
        {
          // Annulus bounds
          constexpr auto annulus_context = typed_detector::surface_annulus::mask_context;
          // Get the trapezoid mask container
          auto &annulus_masks = std::get<annulus_context>(c_masks);
          dindex annulus_index = annulus_masks.size();
          annulus_masks.push_back({io_surface.bound_param0,
                                   io_surface.bound_param1,
                                   io_surface.bound_param2,
                                   io_surface.bound_param3,
                                   io_surface.bound_param4,
                                   io_surface.bound_param5,
                                   io_surface.bound_param6});
          // The read is valid: set the index
          mask_index = {annulus_context, annulus_index};
        }

        // Fill the surface into the temporary container
        if (mask_index[0] != dindex_invalid)
        {
          dindex transform_index = surface_transform_storage.size();
          surface_transform_storage.push_back(transform3{t, z, x});
          c_surfaces.push_back({transform_index, mask_index, c_volume->index(), io_surface.geometry_id});
        }

      } // end of exclusion for navigation layers
    }

    /** Helper method to sort and remove duplicates
     * 
     * @param att attribute vector for sorting and duplicate removal
     * 
     * @return the key values
     */
    auto
        sort_and_remove_duplicates = [](std::map<scalar, std::vector<dindex>> &att) -> dvector<scalar>
    {
      dvector<scalar> keys;
      keys.reserve(att.size());
      for (auto [key, value] : att)
      {
        keys.push_back(key);
        std::sort(value.begin(), value.end());
        value.erase(std::unique(value.begin(), value.end()), value.end());
      }
      return keys;
    };

    // Drawing the lines for the grid search
    auto rs = sort_and_remove_duplicates(r_min_attachments);
    rs.push_back(r_max);
    auto zs = sort_and_remove_duplicates(z_min_attachments);
    zs.push_back(z_max);

    // Create axes and volume grid
    axis::irregular raxis{{rs}};
    axis::irregular zaxis{{zs}};

    typename typed_detector::volume_grid v_grid(std::move(raxis), std::move(zaxis));

    // A step into the volume (stepsilon), can be read in from the smallest difference
    scalar stepsilon = 1.;

    // Run the bin association and write out
    surface_grid_entries_writer sge_writer("grid-entries.csv");
    bool write_grid_entries = (grid_entries_file_name.find("write") != std::string::npos);
    bool read_grid_entries = not grid_entries_file_name.empty() and not write_grid_entries
          and not (grid_entries_file_name.find("none") != std::string::npos);

    // Loop over the volumes
    // - fill the volume grid
    // - run the bin association
    for (auto [iv, v] : enumerate(d.volumes()))
    {
      // Get the volume bounds for filling
      const auto &v_bounds = v.bounds();

      dindex irl = v_grid.axis_p0().bin(v_bounds[0] + stepsilon);
      dindex irh = v_grid.axis_p0().bin(v_bounds[1] - stepsilon);
      dindex izl = v_grid.axis_p1().bin(v_bounds[2] + stepsilon);
      dindex izh = v_grid.axis_p1().bin(v_bounds[3] - stepsilon);
      dindex volume_index = v.index();

      auto r_low = v_grid.axis_p0().borders(irl)[0];
      auto r_high = v_grid.axis_p0().borders(irh)[1];
      auto z_low = v_grid.axis_p1().borders(izl)[0];
      auto z_high = v_grid.axis_p1().borders(izh)[1];

      bool is_cylinder = std::abs(v_bounds[1] - v_bounds[0]) < std::abs(v_bounds[3] - v_bounds[2]);

      for (dindex ir = irl; ir <= irh; ++ir)
      {
        for (dindex iz = izl; iz <= izh; ++iz)
        {
          v_grid.populate(ir, iz, std::move(volume_index));
        }
      }

      dindex sfi = v.surfaces_finder_entry();
      if (sfi != dindex_invalid and write_grid_entries)
      {
        auto &grid = is_cylinder ? detector_surfaces_finders[sfi + 2] : detector_surfaces_finders[sfi];
        bin_association(surface_default_context, v, grid, {0.1, 0.1}, false);

        csv_surface_grid_entry csv_ge;
        csv_ge.detray_volume_id = static_cast<int>(iv);
        size_t nbins0 = grid.axis_p0().bins();
        size_t nbins1 = grid.axis_p1().bins();
        for (size_t b0 = 0; b0 < nbins0; ++b0)
        {
          for (size_t b1 = 0; b1 < nbins1; ++b1)
          {
            csv_ge.detray_bin0 = b0;
            csv_ge.detray_bin1 = b1;
            for (auto e : grid.bin(b0, b1))
            {
              csv_ge.detray_entry = e;
              sge_writer.append(csv_ge);
            }
          }
        }
      }
    }

    // Fast option, read the grid entries back in
    if (read_grid_entries)
    {

      surface_grid_entries_reader sge_reader(grid_entries_file_name);
      csv_surface_grid_entry surface_grid_entry;
      while (sge_reader.read(surface_grid_entry))
      {
        // Get the volume bounds for fillind
        const auto &v = d.indexed_volume(surface_grid_entry.detray_volume_id);
        const auto &v_bounds = v.bounds();
        dindex sfi = v.surfaces_finder_entry();
        if (sfi != dindex_invalid)
        {
          bool is_cylinder = std::abs(v_bounds[1] - v_bounds[0]) < std::abs(v_bounds[3] - v_bounds[2]);
          auto &grid = is_cylinder ? detector_surfaces_finders[sfi + 2] : detector_surfaces_finders[sfi];
          // Fill the entry
          grid.populate(static_cast<dindex>(surface_grid_entry.detray_bin0),
                        static_cast<dindex>(surface_grid_entry.detray_bin1),
                        static_cast<dindex>(surface_grid_entry.detray_entry));
        }
      }
    }

    // Connect the cylindrical volumes
    connect_cylindrical_volumes<typed_detector, array_type, tuple_type, vector_type>(d, v_grid);

    // Add the surface finders to the detector
    d.add_surfaces_finders(std::move(detector_surfaces_finders));

    // Add the volume grid to the detector
    d.add_volume_grid(std::move(v_grid));

    return d;
  }

} // namespace detray
