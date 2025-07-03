#include "plan_manage/traj_planner.hpp"

#define USE_SEARCHTIME false // there is a bug with searchTime

namespace plan_manage
{

    TrajPlanner::TrajPlanner(const std::string &config_string, VoxelMap::Ptr &map_ptr, VehicleGeometry::Ptr &vehicle_geo_ptr)
    {

        YAML::Node config = YAML::Load(config_string);
        auto search_config = config["search"];
        auto optimization_config = config["optimization"];

        collision_check_type_ = search_config["collision_check_type"].as<int>(0);
        non_siguav_ = search_config["non_siguav"].as<double>(0.05);

        traj_piece_duration_ = optimization_config["traj_piece_duration"].as<double>(1.0);
        traj_res_ = optimization_config["traj_resolution"].as<int>(8);
        dense_traj_res_ = optimization_config["dense_traj_resolution"].as<int>(20);
        gear_opt_ = optimization_config["gear_opt"].as<bool>(false);
        verbose_ = optimization_config["verbose"].as<bool>(false);
        use_firi_ = optimization_config["use_firi"].as<bool>(false);
        sfc_range_ = optimization_config["sfc_range_"].as<double>(10.0);

        // Voxel Map
        setMap(map_ptr);

        // Vehicle Geometry
        setVehicleGeometry(vehicle_geo_ptr);

        // Kinodynamic a* path finder
        kino_astar_path_finder_ = std::make_shared<KinoAstar>(config_string, map_ptr, vehicle_geo_ptr);

        // Polynomial trajectory optimizer
        poly_traj_opt_ = std::make_shared<PolyTrajOptimizer>(config_string, vehicle_geo_ptr);
    }

    // use kinodynamic a* to generate a path
    bool TrajPlanner::generateKinoPath(const Eigen::Vector4d &start, const Eigen::Vector4d &end)
    {
        Eigen::Vector4d start_state, end_state; // x, y, yaw, vel
        start_state = start;
        end_state = end;
        if (start_state(3) == 0.0)
        {
            start_state(3) = 0.05; // if velocity is 0, assume it's forward
        }
        else
        {
            if (fabs(start_state(3)) < non_siguav_)
                start_state(3) = non_siguav_ * ((start_state(3) < 0) ? -1.0 : 1.0);
        }
        if (verbose_)
            printf("[Traj Planner]: collision_check_type : %d \n", collision_check_type_);

        Eigen::Vector2d init_ctrl(0.0, 0.0); // steer and acc
        int status;
        Eigen::Vector2d start_pos = start_state.head(2);
        kino_astar_path_finder_->findNearestNode(start_pos, true);
        kino_astar_path_finder_->reset();
        chrono::time_point<std::chrono::high_resolution_clock> search_t1, search_t2;
        if (collision_check_type_ == 0)
        {
            search_t1 = chrono::high_resolution_clock::now();
            status = kino_astar_path_finder_->search(start_state, init_ctrl, end_state, true);
            search_t2 = chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> search_ms = search_t2 - search_t1;
            if (verbose_)
                printf("[Traj Planner]: kino astar search time: %f ms\n", search_ms.count());
        }
        else if (collision_check_type_ == 1)
        {
            search_t1 = chrono::high_resolution_clock::now();
            status = kino_astar_path_finder_->searchUseCircle(start_state, init_ctrl, end_state, true);
            search_t2 = chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> search_ms = search_t2 - search_t1;
            if (verbose_)
                printf("[Traj Planner]: kino astar circle_search time: %f ms\n", search_ms.count());
        }
        else
        {
            printf("[Traj Planner]: kino astar search error:");
        }

        if (status == KinoAstar::NO_PATH)
        {
            if (verbose_)
                printf("[Traj Planner]: kino astar search failed!\n");
            // retry searching with discontinuous initial state
            kino_astar_path_finder_->reset();
            status = kino_astar_path_finder_->search(start_state, init_ctrl, end_state, false);
            if (status == KinoAstar::NO_PATH)
            {
                if (verbose_)
                    printf("[Traj Planner]: kino astar couldn't find a path.\n");
                return false;
            }
            else
            {
                if (verbose_)
                    printf("[Traj Planner]: kino astar retry search succeeded.\n");
            }
        }
        else
        {
            if (verbose_)
                printf("[Traj Planner]: kino astar search success.\n");
        }

        auto time_searcht1 = chrono::high_resolution_clock::now();
        kino_astar_path_finder_->getFullposeLists();
        kino_astar_path_finder_->getSingulNodes();
        /* Use searchTime */
        if (USE_SEARCHTIME)
        {
            bool middle_node_success = kino_astar_path_finder_->searchTime(kino_trajs_);
        }
        /* Use getKinoNode */
        else
        {
            kino_astar_path_finder_->getKinoNode(kino_trajs_);
        }
        bool middle_node_success = true;
        auto time_searcht2 = chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> time_searchms = time_searcht2 - time_searcht1;
        // if (verbose_)
        //     printf("[Traj Planner]: kino astar speed planning search time: %f ms\n", time_searchms.count());

        if (!middle_node_success)
        {
            if (verbose_)
                printf("[Traj Planner]: kino astar speed planning failed!\n");
            return false;
        }
        else
        {
            if (verbose_)
                printf("[Traj Planner]: kino astar speed planning succeeded!\n");
            return true;
        }
    }

    bool TrajPlanner::RunMINCO()
    {
        Eigen::MatrixXd flat_finalState(2, 3), flat_headState(2, 3);
        Eigen::VectorXd ego_piece_dur_vec;
        Eigen::MatrixXd ego_innerPs;

        /* try to merge optimization process */
        sfc_state_container_.clear();
        sfc_container_.clear();
        aux_sfc_container_.clear();
        std::vector<int> singul_container;
        Eigen::VectorXd duration_container;
        std::vector<Eigen::MatrixXd> waypoints_container;
        std::vector<Eigen::MatrixXd> iniState_container, finState_container;
        duration_container.resize(kino_trajs_.size());
        double base_time = 0.0;
        double total_corridor_ms = 0.0;

        for (unsigned int seg_idx = 0; seg_idx < kino_trajs_.size(); seg_idx++)
        {
            double timePerPiece = traj_piece_duration_;
            plan_utils::FlatTrajData kino_traj = kino_trajs_.at(seg_idx);
            singul_container.push_back(kino_traj.singul);
            std::vector<Eigen::Vector3d> pts = kino_traj.traj_pts;
            int piece_nums = 0;
            double initTotalduration = 0.0;
            for (const auto &pt : pts)
            {
                initTotalduration += pt[2];
            }
            piece_nums = std::max(static_cast<int>(std::ceil(initTotalduration / timePerPiece)), 2);
            timePerPiece = initTotalduration / piece_nums;
            ego_piece_dur_vec.resize(piece_nums);
            ego_piece_dur_vec.setConstant(timePerPiece);
            duration_container[seg_idx] = initTotalduration;
            ego_innerPs.resize(2, piece_nums - 1);
            std::vector<Eigen::Vector3d> state_list;
            double piece_time = 0;
            for (int i = 0; i < piece_nums; i++)
            {
                int resolution;
                if (i == 0 || i == piece_nums - 1)
                {
                    resolution = dense_traj_res_;
                }
                else
                {
                    resolution = traj_res_;
                }
                for (int k = 0; k <= resolution; k++)
                {
                    double t;
                    Eigen::Vector3d pos;
                    if (USE_SEARCHTIME)
                    {
                        /* Use searchTime */
                        t = piece_time + 1.0 * k / resolution * ego_piece_dur_vec[i];
                        pos = kino_astar_path_finder_->CalculateInitPos(t, kino_traj.singul);
                    }
                    else
                    {
                        /* Use getKinoNode */
                        t = base_time + piece_time + 1.0 * k / resolution * ego_piece_dur_vec[i];
                        pos = kino_astar_path_finder_->evaluatePos(t);
                    }
                    state_list.push_back(pos);
                    if (k == resolution && i != piece_nums - 1)
                    {
                        ego_innerPs.col(i) = pos.head(2);
                    }
                }
                piece_time += ego_piece_dur_vec[i];
            }
            sfc_state_container_.push_back(state_list);

            // calculate safe corridors
            auto corridor_t1 = chrono::high_resolution_clock::now();
            std::vector<Eigen::MatrixXd> hPolys;                  // for car
            std::vector<std::vector<Eigen::MatrixXd>> aux_hPolys; // for auxiliary
            /** Option 1: FIRI corridor **/
            if (use_firi_)
            {
                calcFIRICorridor(hPolys, state_list, map_ptr_, vehicle_geo_ptr_, sfc_range_);
                calcAuxFIRICorridor(aux_hPolys, state_list, map_ptr_, vehicle_geo_ptr_, sfc_range_);
            }
            /** Option 2: Rectangle corridor **/
            else
            {
                calcRectangleCorridor(hPolys, state_list, map_ptr_, vehicle_geo_ptr_, sfc_range_);
                calcAuxRectangleCorridor(aux_hPolys, state_list, map_ptr_, vehicle_geo_ptr_, sfc_range_);
            }
            sfc_container_.push_back(hPolys);
            aux_sfc_container_.push_back(aux_hPolys);
            auto corridor_t2 = chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> corridor_ms = corridor_t2 - corridor_t1;
            total_corridor_ms += corridor_ms.count();

            waypoints_container.push_back(ego_innerPs);
            iniState_container.push_back(kino_traj.start_state);
            finState_container.push_back(kino_traj.final_state);
            base_time += initTotalduration;
        }

        if (verbose_)
            printf("[Traj Planner]: corridor computing time: %f ms\n", total_corridor_ms);

        // debug
        iniStates_ = iniState_container;
        finStates_ = finState_container;
        initInnerPts_ = waypoints_container;
        initTs_ = duration_container;
        singuls_ = singul_container;

        // start optimization
        bool flag_success;
        flag_success = poly_traj_opt_->OptimizeTrajectory(iniState_container, finState_container,
                                                          waypoints_container, duration_container,
                                                          sfc_container_,
                                                          aux_sfc_container_,
                                                          singul_container,
                                                          verbose_);

        if (flag_success)
        {
            traj_container_.clearSingul();
            double baseTime = 0.0;
            Eigen::VectorXd optimizedDuration(poly_traj_opt_->getTrajNum());
            for (unsigned int i = 0; i < poly_traj_opt_->getTrajNum(); i++)
            {
                traj_container_.addSingulTraj((*poly_traj_opt_->getMinJerkOptPtr())[i].getTraj(poly_traj_opt_->getSinguls()[i]), baseTime);
                optimizedDuration[i] = (*poly_traj_opt_->getMinJerkOptPtr())[i].getTotalDuration();
                // if (verbose_)
                // {
                //     printf("[Traj Planner]: seg%d init duration: %f, optimized duration: %f\n", i, duration_container[i], (*poly_traj_opt_->getMinJerkOptPtr())[i].getTotalDuration());
                //     printf("[Traj Planner]: seg%d pieceNum: %d\n", i, waypoints_container[i].cols() + 1);
                //     printf("[Traj Planner]: seg%d optimized jerk cost: %f\n", i, (*poly_traj_opt_->getMinJerkOptPtr())[i].getTrajJerkCost());
                // }

                baseTime = traj_container_.singul_traj.back().end_time;
            }
            if (verbose_)
            {
                printf("[Traj Planner]: init traj duration: %f s, optimal traj duration: %f s\n", duration_container.sum(), optimizedDuration.sum());
                printf("[Traj Planner]: Planning success!\n");
            }
            return true;
        }
        else
        {
            if (verbose_)
                printf("[Traj Planner] Planning failed!\n");
            return false;
        }
    }

    bool TrajPlanner::checkCollisionUsingPosAndYaw(const Eigen::Vector3d &state)
    {
        bool res;
        kino_astar_path_finder_->checkCollisionUsingPosAndYaw(state, res);

        return res;
    }

} // namespace plan_manage
