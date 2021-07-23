/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include <dfe/dfe_namedtuple.hpp>
#include <dfe/dfe_io_dsv.hpp>
#include <vecmem/containers/vector.hpp>

#include <fstream>

namespace detray
{

  template<typename T>
  using vecm_vector = vecmem::vector<T>;

  struct csv_surface
  {
    /// Surface identifier. Not available in the TrackML datasets.
    uint64_t geometry_id;
    /// Partially decoded surface identifier components.
    uint32_t volume_id, boundary_id, layer_id, module_id;
    /// Center position components in mm.
    scalar cx, cy, cz;
    /// Rotation matrix components.
    scalar rot_xu, rot_xv, rot_xw;
    scalar rot_yu, rot_yv, rot_yw;
    scalar rot_zu, rot_zv, rot_zw;
    /// The type of the surface bpounds object, determines the parameters filled
    int bounds_type;
    scalar bound_param0 = -1.f;
    scalar bound_param1 = -1.f;
    scalar bound_param2 = -1.f;
    scalar bound_param3 = -1.f;
    scalar bound_param4 = -1.f;
    scalar bound_param5 = -1.f;
    scalar bound_param6 = -1.f;

    DFE_NAMEDTUPLE(csv_surface, geometry_id, volume_id, boundary_id, layer_id,
                   module_id, cx, cy, cz, rot_xu, rot_xv, rot_xw, rot_yu, rot_yv,
                   rot_yw, rot_zu, rot_zv, rot_zw, bounds_type, bound_param0,
                   bound_param1, bound_param2, bound_param3, bound_param4,
                   bound_param5, bound_param6);
  };

  using surface_reader = dfe::NamedTupleCsvReader<csv_surface>;

  struct csv_surface_grid
  {
    /// Surface identifier. Not available in the TrackML datasets.
    uint64_t geometry_id;
    /// Partially decoded surface identifier components.
    uint32_t volume_id, layer_id, surface_id;
    /// The number of bins in loc 0 / 1
    int type_loc0 = -1;
    int nbins_loc0 = -1;
    scalar min_loc0, max_loc0;
    int type_loc1 = -1;
    int nbins_loc1 = -1;
    scalar min_loc1, max_loc1;

    DFE_NAMEDTUPLE(csv_surface_grid, geometry_id, volume_id, layer_id, surface_id, type_loc0, nbins_loc0, min_loc0,
                   max_loc0, type_loc1, nbins_loc1, min_loc1, max_loc1);
  };

  using surface_grid_reader = dfe::NamedTupleCsvReader<csv_surface_grid>;


  struct csv_surface_grid_entry
  {
    /// detray volume identifier
    int detray_volume_id = -1;
    int detray_bin0 = -1;
    int detray_bin1 = -1;
    int detray_entry = -1;

    DFE_NAMEDTUPLE(csv_surface_grid_entry, detray_volume_id, detray_bin0, detray_bin1, detray_entry);
  };

  using surface_grid_entries_writer = dfe::NamedTupleCsvWriter<csv_surface_grid_entry>;
  using surface_grid_entries_reader = dfe::NamedTupleCsvReader<csv_surface_grid_entry>;


  struct csv_layer_volume
  {
    /// Surface identifier. Not available in the TrackML datasets.
    uint64_t geometry_id;
    /// Partially decoded surface identifier components.
    uint32_t volume_id, layer_id;
    /// The type of the surface bpounds object, determines the parameters filled
    int volume_type;
    float min_v0 = -1.f;
    float max_v0 = -1.f;
    float min_v1 = -1.f;
    float max_v1 = -1.f;
    float min_v2 = -1.f;
    float max_v2 = -1.f;

    DFE_NAMEDTUPLE(csv_layer_volume, geometry_id, volume_id, layer_id, min_v0,
                   max_v0, min_v1, max_v1, min_v2, max_v2);
  };

  using layer_volume_reader = dfe::NamedTupleCsvReader<csv_layer_volume>;

} // namespace detray
