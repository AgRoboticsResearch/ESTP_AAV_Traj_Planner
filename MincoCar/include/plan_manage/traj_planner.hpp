#ifndef _TRAJ_MANAGER_HPP_
#define _TRAJ_MANAGER_HPP_

#include <memory>
#include <set>
#include <string>
#include <iostream>
#include <cmath>
#include <chrono>
#include <yaml-cpp/yaml.h>

#include "map_utils/voxel_map.hpp"
#include "vehicle_utils/vehicle_geometry.hpp"
#include "corridor_utils/safe_corridor.hpp"
#include "path_searching/kino_astar.hpp"
#include "plan_manage/traj_optimizer.hpp"

namespace plan_manage
{

    using namespace map_utils;
    using namespace vehicle_utils;
    using namespace corridor_utils;
    using namespace path_searching;

    class TrajPlanner
    {
    public:
        TrajPlanner(){};
        TrajPlanner(const std::string &config_string, VoxelMap::Ptr &map_ptr, VehicleGeometry::Ptr &vehicle_geo_ptr);
        ~TrajPlanner(){};
        typedef std::shared_ptr<TrajPlanner> Ptr;

        /* main functions */
        void setMap(const VoxelMap::Ptr &ptr) { map_ptr_ = ptr; };
        void setVehicleGeometry(const VehicleGeometry::Ptr &ptr) { vehicle_geo_ptr_ = ptr; };
        bool generateKinoPath(const Eigen::Vector4d &start, const Eigen::Vector4d &end);
        bool RunMINCO();
        bool checkCollisionUsingPosAndYaw(const Eigen::Vector3d &state);
        
        /* helper functions */
        KinoAstar::Ptr getKinoAstarPathFinder() const { return kino_astar_path_finder_; };
        plan_utils::KinoTrajData getFlatTraj() const { return kino_trajs_; };
        std::vector<std::vector<Eigen::Vector3d>> getCorridorState() const { return sfc_state_container_; };
        std::vector<std::vector<Eigen::MatrixXd>> getCorridor() const { return sfc_container_; };
        std::vector<std::vector<std::vector<Eigen::MatrixXd>>> getAuxCorridor() const { return aux_sfc_container_; };
        plan_utils::SingulTrajData getOptTraj() const { return traj_container_.singul_traj; };
        std::vector<Eigen::Vector3d> getConstraintPts() const { return poly_traj_opt_->getConstraintPts(); };
        std::vector<Eigen::VectorXd> getCostsHistory() const { return poly_traj_opt_->getCostsHistory(); }; 

        /* debug */
        std::vector<Eigen::MatrixXd> getIniStates() const { return iniStates_; };
        std::vector<Eigen::MatrixXd> getFinStates() const { return finStates_; };
        std::vector<Eigen::MatrixXd> getInitInnerPts() const { return initInnerPts_; };
        Eigen::VectorXd getInitTs() const { return initTs_; };
        std::vector<int> getSinguls() const { return singuls_; };

    private:
        // map
        VoxelMap::Ptr map_ptr_;

        // vehicle geometry
        VehicleGeometry::Ptr vehicle_geo_ptr_;

        // front-end path planning
        KinoAstar::Ptr kino_astar_path_finder_;
        plan_utils::KinoTrajData kino_trajs_;

        // safe corridor
        std::vector<std::vector<Eigen::Vector3d>> sfc_state_container_;
        std::vector<std::vector<Eigen::MatrixXd>> sfc_container_;
        std::vector<std::vector<std::vector<Eigen::MatrixXd>>> aux_sfc_container_;
        bool use_firi_;
        double sfc_range_;

        // back-end trajectory planning
        plan_utils::TrajContainer traj_container_;
        plan_manage::PolyTrajOptimizer::Ptr poly_traj_opt_;
        double traj_piece_duration_;
        int traj_res_, dense_traj_res_;
        bool gear_opt_;
        double non_siguav_;
        bool collision_check_type_;

        // debug
        bool verbose_;
        std::vector<Eigen::MatrixXd> iniStates_;
        std::vector<Eigen::MatrixXd> finStates_;
        std::vector<Eigen::MatrixXd> initInnerPts_;
        Eigen::VectorXd initTs_;
        std::vector<int> singuls_;

    };

} // namespace plan_manage

#endif // _TRAJ_MANAGER_HPP_
