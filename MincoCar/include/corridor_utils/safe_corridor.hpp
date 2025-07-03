#ifndef _SAFE_CORRIDOR_HPP_
#define _SAFE_CORRIDOR_HPP_

#include <vector>
#include <Eigen/Eigen>

#include "map_utils/voxel_map.hpp"
#include "vehicle_utils/vehicle_geometry.hpp"

namespace corridor_utils
{
    /* Calculate corridors using rectangles */
    void calcRectangleCorridor(std::vector<Eigen::MatrixXd> &hPolys,
                               const std::vector<Eigen::Vector3d> &state_list,
                               const map_utils::VoxelMap::Ptr &map_ptr,
                               const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                               double limitBound = 10.0);

    /* Calculate auxiliary corridors using rectangles */
    void calcAuxRectangleCorridor(std::vector<std::vector<Eigen::MatrixXd>> &hPolys_list,
                                  const std::vector<Eigen::Vector3d> &state_list,
                                  const map_utils::VoxelMap::Ptr &map_ptr,
                                  const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                                  double limitBound = 10.0);

    /* Generate corridors using FIRI */
    void calcFIRICorridor(std::vector<Eigen::MatrixXd> &hPolys,
                          const std::vector<Eigen::Vector3d> &state_list,
                          const map_utils::VoxelMap::Ptr &map_ptr,
                          const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                          double range = 10.0);

    /* Generate auxiliary corridors using FIRI */
    void calcAuxFIRICorridor(std::vector<std::vector<Eigen::MatrixXd>> &hPolys_list,
                             const std::vector<Eigen::Vector3d> &state_list,
                             const map_utils::VoxelMap::Ptr &map_ptr,
                             const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                             double range = 10.0);

} // namespace corridor_utils

#endif // _SAFE_CORRIDOR_HPP_