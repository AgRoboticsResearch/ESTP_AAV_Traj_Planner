#include "plan_manage/traj_optimizer.hpp"

namespace plan_manage
{
    int PolyTrajOptimizer::earlyExitCallback(void *func_data, const double *x, const double *g, const double fx, const double xnorm, const double gnorm, const double step, int n, int k, int ls)
    {
        PolyTrajOptimizer *opt = reinterpret_cast<PolyTrajOptimizer *>(func_data);

        return (opt->force_stop_type_ == STOP_FOR_ERROR || opt->force_stop_type_ == STOP_FOR_REBOUND);
    }

    /* Mappings between real world time and unconstrained virtual time */
    template <typename EIGENVEC>
    void PolyTrajOptimizer::RealT2VirtualT(const Eigen::VectorXd &RT, EIGENVEC &VT)
    {
        for (int i = 0; i < RT.size(); ++i)
        {
            VT(i) = RT(i) > 1.0 ? (sqrt(2.0 * RT(i) - 1.0) - 1.0)
                                : (1.0 - sqrt(2.0 / RT(i) - 1.0));
        }
    }

    void PolyTrajOptimizer::RealT2VirtualT(double &RT, double &VT)
    {
        VT = RT > 1.0 ? (sqrt(2.0 * RT - 1.0) - 1.0)
                      : (1.0 - sqrt(2.0 / RT - 1.0));
    }

    template <typename EIGENVEC>
    void PolyTrajOptimizer::VirtualT2RealT(const EIGENVEC &VT, Eigen::VectorXd &RT)
    {
        for (int i = 0; i < VT.size(); ++i)
        {
            RT(i) = VT(i) > 0.0 ? ((0.5 * VT(i) + 1.0) * VT(i) + 1.0)
                                : 1.0 / ((0.5 * VT(i) - 1.0) * VT(i) + 1.0);
        }
    }

    void PolyTrajOptimizer::VirtualT2RealT(double &VT, double &RT)
    {
        RT = VT > 0.0 ? ((0.5 * VT + 1.0) * VT + 1.0)
                      : 1.0 / ((0.5 * VT - 1.0) * VT + 1.0);
    }

    template <typename EIGENVEC, typename EIGENVECGD>
    void PolyTrajOptimizer::VirtualTGradCost(
        const Eigen::VectorXd &RT, const EIGENVEC &VT,
        const Eigen::VectorXd &gdRT, EIGENVECGD &gdVT)
    {
        for (int i = 0; i < VT.size(); ++i)
        {
            double gdVT2Rt;
            if (VT(i) > 0)
            {
                gdVT2Rt = VT(i) + 1.0;
            }
            else
            {
                double denSqrt = (0.5 * VT(i) - 1.0) * VT(i) + 1.0;
                gdVT2Rt = (1.0 - VT(i)) / (denSqrt * denSqrt);
            }

            gdVT(i) = gdRT(i) * gdVT2Rt;
        }
    }

    void PolyTrajOptimizer::VirtualTGradCost(const double &RT, const double &VT, const double &gdRT, double &gdVT)
    {
        double gdVT2Rt;
        if (VT > 0)
        {
            gdVT2Rt = VT + 1.0;
        }
        else
        {
            double denSqrt = (0.5 * VT - 1.0) * VT + 1.0;
            gdVT2Rt = (1.0 - VT) / (denSqrt * denSqrt);
        }

        gdVT = gdRT * gdVT2Rt;
    }

    /* First-order postive smoothing approximation */
    void PolyTrajOptimizer::positiveSmoothedL1(const double &x, double &f, double &df)
    {
        // viola*, viola*Pena, viola*PenaD
        double pe = 1.0e-4;
        double f3c = 1.0 / (pe * pe);
        double f4c = -0.5 * f3c / pe;

        if (x <= 0)
        {
            f = 0.0;
            df = 0.0;
        }
        else if (x <= pe)
        {
            f = (f4c * x + f3c) * x * x * x;
            df = (4.0 * f4c * x + 3.0 * f3c) * x * x;
        }
        else
        {
            f = x - 0.5 * pe;
            df = 1.0;
        }
        return;
    }

    /* Third-order postive smoothing approximation */
    void PolyTrajOptimizer::positiveSmoothedL3(const double &x, double &f, double &df)
    {
        if (x <= 0)
        {
            f = 0.0;
            df = 0.0;
        }
        else
        {
            f = x * x * x;
            df = 3.0 * x * x;
        }
    }

    // double PolyTrajOptimizer::log_sum_exp(double alpha, Eigen::VectorXd &all_dists, double &exp_sum)
    // {
    //     // all_dists will be std::exp(alpha * (all_dists(j) - d_max));
    //     double d_0;
    //     if (alpha > 0)
    //     {
    //         d_0 = all_dists.maxCoeff();
    //     }
    //     else
    //     {
    //         d_0 = all_dists.minCoeff();
    //     }

    //     exp_sum = 0;
    //     for (unsigned int j = 0; j < all_dists.size(); j++)
    //     {
    //         all_dists(j) = std::exp(alpha * (all_dists(j) - d_0));
    //         exp_sum += all_dists(j);
    //     }

    //     return std::log(exp_sum) / alpha + d_0;
    // }

    PolyTrajOptimizer::PolyTrajOptimizer(const std::string &config_string, VehicleGeometry::Ptr &vehicle_geo_ptr)
    {
        YAML::Node config = YAML::Load(config_string);
        auto optimization_config = config["optimization"];
        auto search_config = config["search"];

        traj_res_ = optimization_config["traj_resolution"].as<int>(8);
        dense_traj_res_ = optimization_config["dense_traj_resolution"].as<int>(20);
        wei_obs_ = optimization_config["wei_sta_obs"].as<double>(7000.0);
        wei_feas_ = optimization_config["wei_feas"].as<double>(1000.0);
        wei_time_ = optimization_config["wei_time"].as<double>(500.0);
        max_vel_ = optimization_config["max_vel"].as<double>(5.0);
        max_acc_ = optimization_config["max_acc"].as<double>(5.0);
        max_cur_ = optimization_config["max_cur"].as<double>(0.350877);
        max_phi_dot_ = optimization_config["max_phi_dot"].as<double>(5.0);
        max_latacc_ = optimization_config["max_latacc"].as<double>(5.0);
        max_omg_ = optimization_config["max_omg"].as<double>(10000.0);
        non_siguav_ = search_config["non_siguav"].as<double>(0.05);
        epis_ = optimization_config["epis"].as<double>(1.0e-4);
        gear_opt_ = optimization_config["gear_opt"].as<bool>(true);
        remove_gear_ = optimization_config["remove_gear"].as<bool>(false);
        log_cost_ = optimization_config["log_cost"].as<bool>(false);
        // L-BFGS parameters
        lbfgs_mem_size = optimization_config["lbfgs_mem_size"].as<int>(128);
        lbfgs_g_epsilon = optimization_config["lbfgs_g_epsilon"].as<double>(1.0e-6);
        lbfgs_past = optimization_config["lbfgs_past"].as<int>(3);
        lbfgs_delta = optimization_config["lbfgs_delta"].as<double>(1.0e-6);
        lbfgs_max_linesearch = optimization_config["lbfgs_max_linesearch"].as<int>(100);
        lbfgs_max_iterations = optimization_config["lbfgs_max_iterations"].as<int>(0);
        lbfgs_f_dec_coeff = optimization_config["lbfgs_f_dec_coeff"].as<double>(1.0e-4);
        lbfgs_s_curv_coeff = optimization_config["lbfgs_s_curv_coeff"].as<double>(0.9);
        lbfgs_cautious_factor = optimization_config["lbfgs_cautious_factor"].as<double>(1.0e-6);

        car_wheelbase_ = vehicle_geo_ptr->getWheelbase();

        B_h << 0, -1,
            1, 0;

        // vehicle vertices
        ego_vertices_.clear();
        for (const auto &vertex : vehicle_geo_ptr->getCarVertices())
        {
            ego_vertices_.push_back(vertex);
        }

        if (ego_vertices_.front() != ego_vertices_.back())
        {
            ego_vertices_.push_back(ego_vertices_.front());
        }

        // aux vertices
        aux_vertices_list_.clear();
        aux_vertices_list_.reserve(vehicle_geo_ptr->getAuxCount());
        for (int i = 0; i < vehicle_geo_ptr->getAuxCount(); i++)
        {
            std::vector<Eigen::Vector2d> aux_vertices = vehicle_geo_ptr->getAuxVertices(i);
            std::vector<Eigen::Vector2d> vertices;
            for (const auto &vertex : aux_vertices)
            {
                vertices.push_back(vertex);
            }
            if (vertices.front() != vertices.back())
            {
                vertices.push_back(vertices.front());
            }
            aux_vertices_list_.push_back(vertices);
        }
    }

    bool PolyTrajOptimizer::OptimizeTrajectory(
        const std::vector<Eigen::MatrixXd> &iniStates, const std::vector<Eigen::MatrixXd> &finStates,
        const std::vector<Eigen::MatrixXd> &initInnerPts, const Eigen::VectorXd &initTs,
        const std::vector<std::vector<Eigen::MatrixXd>> &hPoly_container,
        const std::vector<std::vector<std::vector<Eigen::MatrixXd>>> &aux_hPoly_container,
        const std::vector<int> &singuls,
        const bool verbose)
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        iniState_container = iniStates;
        finState_container = finStates;
        cfgHs_container = hPoly_container;
        aux_cfgHs_container = aux_hPoly_container;
        singul_container = singuls;

        trajnum_ = initInnerPts.size();
        jerkOpt_container.clear();
        jerkOpt_container.resize(trajnum_);
        piece_num_container.clear();
        piece_num_container.resize(trajnum_);

        int variable_num = 0; // number of optimization variables

        // check trajectory number
        if (initTs.size() != trajnum_)
        {
            if (verbose)
                printf("[Poly Traj Optimizer]: Error! initTs.size() != initInnerPts.size() \n");
            return false;
        }

        for (int i = 0; i < trajnum_; i++)
        {
            // check piece number
            if (initInnerPts[i].cols() == 0)
            {
                if (verbose)
                    printf("[Poly Traj Optimizer]: Error! There is only one piece in segment %i\n", i);
                return false;
            }

            int piece_num = initInnerPts[i].cols() + 1;
            piece_num_container[i] = piece_num;

            // check corridor dimension
            // TODO: remove corridor dimension check
            if (cfgHs_container[i].size() != (piece_num - 2) * (traj_res_ + 1) + 2 * (dense_traj_res_ + 1))
            {
                if (verbose)
                    printf("[Poly Traj Optimizer]: Error! cfgHs size = %zu\n", cfgHs_container[i].size());
                return false;
            }

            // bound initial and final states velocity and acceleration
            if (iniState_container[i].col(1).norm() >= max_vel_)
                iniState_container[i].col(1) = iniState_container[i].col(1).normalized() * (max_vel_ - 1.0e-2);
            if (iniState_container[i].col(2).norm() >= max_acc_)
                iniState_container[i].col(2) = iniState_container[i].col(2).normalized() * (max_acc_ - 1.0e-2);
            if (finState_container[i].col(1).norm() >= max_vel_)
                finState_container[i].col(1) = finState_container[i].col(1).normalized() * (max_vel_ - 1.0e-2);
            if (finState_container[i].col(2).norm() >= max_acc_)
                finState_container[i].col(2) = finState_container[i].col(2).normalized() * (max_acc_ - 1.0e-2);

            // reset optimizer
            jerkOpt_container[i].reset(piece_num);

            /* variables */
            variable_num += 2 * (piece_num - 1); // waypoints
        }
        variable_num += trajnum_; // times

        if (gear_opt_)
        {
            variable_num += 2 * (trajnum_ - 1); // gear positions
            variable_num += 1 * (trajnum_ - 1); // angles
        }
        if (verbose)
        {
            if (gear_opt_)
                printf("[Poly Traj Optimizer]: gear pos opt is enabled!\n");
            else
                printf("[Poly Traj Optimizer]: gear pos opt is disabled!\n");
        }

        /* ---------- L-BFGS parameters ---------- */
        lbfgs::lbfgs_parameter_t lbfgs_params;
        lbfgs_params.mem_size = lbfgs_mem_size;               // 128
        lbfgs_params.g_epsilon = lbfgs_g_epsilon;             // 1.0e-6
        lbfgs_params.past = lbfgs_past;                       // 3
        lbfgs_params.delta = lbfgs_delta;                     // 1.0e-6
        lbfgs_params.max_linesearch = lbfgs_max_linesearch;   // 100
        lbfgs_params.max_iterations = lbfgs_max_iterations;   // 0
        lbfgs_params.f_dec_coeff = lbfgs_f_dec_coeff;         // 1.0e-4
        lbfgs_params.s_curv_coeff = lbfgs_s_curv_coeff;       // 0.9
        lbfgs_params.cautious_factor = lbfgs_cautious_factor; // 1.0e-6

        /* ---------- Prepare ---------- */
        double final_cost = 0;
        iter_num_ = 0;
        costs_history_container.clear();
        bool flag_success = false;
        force_stop_type_ = DONT_STOP;

        /* ---------- Initial guess ---------- */
        Eigen::VectorXd x(variable_num);
        int offset = 0;
        // waypoints
        for (int i = 0; i < trajnum_; i++)
        {
            memcpy(x.data() + offset, initInnerPts[i].data(), initInnerPts[i].size() * sizeof(x[0]));
            offset += initInnerPts[i].size();
        }
        // times
        Eigen::Map<Eigen::VectorXd> Vt(x.data() + offset, initTs.size());
        RealT2VirtualT(initTs, Vt);
        offset += initTs.size();

        if (gear_opt_)
        {
            // gear positions
            for (int i = 0; i < trajnum_ - 1; i++)
            {
                memcpy(x.data() + offset, finState_container[i].col(0).data(), 2 * sizeof(x[0]));
                offset += 2;
            }
            // angles
            Eigen::Map<Eigen::VectorXd> angles(x.data() + offset, trajnum_ - 1);
            for (int i = 0; i < trajnum_ - 1; i++)
            {
                Eigen::Vector2d gearv = finState_container[i].col(1);
                angles[i] = std::atan2(gearv[1], gearv[0]);
            }
        }

        /* ---------- Optimize ---------- */
        // if (verbose)
        //     printf("[Poly Traj Optimizer]: Begin to optimize!\n");
        auto t1 = std::chrono::high_resolution_clock::now();
        int result = lbfgs::lbfgs_optimize(
            x,
            final_cost,
            PolyTrajOptimizer::costFunctionCallback,
            NULL,
            NULL,
            this,
            lbfgs_params);
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> opt_time_ms = t2 - t1;
        std::chrono::duration<double, std::milli> total_time_ms = t2 - t0;

        /* ---------- get result and check collision ---------- */
        if (result == lbfgs::LBFGS_CONVERGENCE ||
            result == lbfgs::LBFGS_CANCELED ||
            result == lbfgs::LBFGS_STOP ||
            result == lbfgs::LBFGSERR_MAXIMUMITERATION ||
            result == lbfgs::LBFGSERR_MAXIMUMLINESEARCH)
        {
            flag_success = true;
            if (verbose)
            {
                if (result == lbfgs::LBFGSERR_MAXIMUMITERATION)
                    printf("[Poly Traj Optimizer]: Lbfgs - Maximum number of iterations has been reached.\n");
                if (result == lbfgs::LBFGSERR_MAXIMUMLINESEARCH)
                    printf("[Poly Traj Optimizer]: Lbfgs - The line-search routine reaches the maximum number of evaluations.\n");

                printf("[Poly Traj Optimizer]: elapsed time: %5.3f ms, iter num: %d, cost: %5.3f\n", total_time_ms.count(), iter_num_, final_cost);
            }

            /* Remove the gear switch if possible */
            if (remove_gear_)
            {
                // /* Print out */
                // offset = 0;
                // // waypoints
                // std::cout << "waypoints: " << std::endl;
                // for (int trajid = 0; trajid < trajnum_; trajid++)
                // {
                //     std::cout << "piece number: " << piece_num_container[trajid] << std::endl;
                //     Eigen::Map<const Eigen::MatrixXd> P(x.data() + offset, 2, piece_num_container[trajid] - 1);
                //     offset += 2 * (piece_num_container[trajid] - 1);
                //     std::cout << "(" << P.transpose() << ")," << std::endl;
                // }

                // // times
                // std::cout << "times: " << std::endl;
                // Eigen::Map<const Eigen::VectorXd> t(x.data() + offset, trajnum_);
                // Eigen::VectorXd T(trajnum_);
                // VirtualT2RealT(t, T);
                // std::cout << T.transpose() << std::endl;
                // offset += trajnum_;

                // if (gear_opt_)
                // {
                //     // gear positions
                //     std::cout << "gear position: " << std::endl;
                //     for (int trajid = 0; trajid < trajnum_ - 1; trajid++)
                //     {
                //         Eigen::Map<const Eigen::MatrixXd> Gear(x.data() + offset, 2, 1);
                //         offset += 2;
                //         std::cout << Gear.transpose() << std::endl;
                //     }
                //     // angles
                //     std::cout << "angles: " << std::endl;
                //     Eigen::Map<const Eigen::VectorXd> Angles(x.data() + offset, trajnum_ - 1);
                //     std::cout << Angles.transpose() << std::endl;
                // }

                // check trajectory length
                if (trajnum_ > 1)
                {
                    double traj_duration = 0;
                    for (int i = 0; i < trajnum_; i++)
                    {
                        traj_duration += jerkOpt_container[i].getTotalDuration();
                    }
                    printf("[Poly Traj Optimizer]: traj duration: %f s\n", traj_duration);
                    
                    std::vector<std::vector<int>> trajid_container;
                    std::vector<int> trajid_list;
                    std::vector<int> skipped_index;
                    std::vector<Eigen::MatrixXd> new_merge_point;

                    int trajid = 0;
                    int offset = 0;
                    int last_singul = singul_container[0]; // !! Note: we assume the first trajectory segment is not short
                    int waypoint_variable_num = 0;
                    for (const auto piece_num : piece_num_container)
                        waypoint_variable_num += 2 * (piece_num - 1);

                    while (trajid < trajnum_)
                    {
                        double traj_length = jerkOpt_container[trajid].getTraj(singul_container[trajid]).getTotalLength();
                        if (traj_length < 1 && trajid != 0) // if the trajectory is too short
                        {
                            std::cout << "[Poly Traj Optimizer]: Skip short trajectory segment: " << trajid << std::endl;
                            for (int i = 0; i < 2 * (piece_num_container[trajid] - 1); i++)
                            {
                                skipped_index.push_back(offset + i); // skipped waypoint index
                            }
                            skipped_index.push_back(waypoint_variable_num + trajid); // skipped time index

                            if (gear_opt_)
                            {
                                int gear_base_index = waypoint_variable_num + trajnum_;
                                skipped_index.push_back(gear_base_index + 2 * (trajid - 1)); // skipped gear position index
                                skipped_index.push_back(gear_base_index + 2 * (trajid - 1) + 1);

                                if (trajid != trajnum_ - 1)
                                {
                                    skipped_index.push_back(gear_base_index + 2 * (trajid - 1) + 2); // skipped gear position index
                                    skipped_index.push_back(gear_base_index + 2 * (trajid - 1) + 3);
                                    // Calculate new position
                                    // Eigen::Map<const Eigen::MatrixXd> gear_before(x.data() + waypoint_variable_num + trajnum_ + 2 * (trajid - 1), 2, 1);
                                    // Eigen::Map<const Eigen::MatrixXd> gear_after(x.data() + waypoint_variable_num + trajnum_ + 2 * trajid, 2, 1);
                                    // Eigen::MatrixXd gear_average = (gear_before + gear_after) / 2;
                                    int k = 0;
                                    for (int i = 0; i < trajid; i++)
                                    {
                                        k += 2 * (piece_num_container[i] - 1);
                                    }
                                    Eigen::Map<const Eigen::MatrixXd> gear_before(x.data() + k - 2, 2, 1);
                                    k = 0;
                                    for (int i = 0; i < trajid + 1; i++)
                                    {
                                        k += 2 * (piece_num_container[i] - 1);
                                    }
                                    Eigen::Map<const Eigen::MatrixXd> gear_after(x.data() + k, 2, 1);
                                    Eigen::MatrixXd gear_average = (gear_before + gear_after) / 2;
                                    new_merge_point.push_back(gear_average);
                                }

                                int angle_base_index = waypoint_variable_num + trajnum_ + 2 * (trajnum_ - 1);
                                skipped_index.push_back(angle_base_index + trajid - 1); // skipped angle index
                                if (trajid != trajnum_ - 1)
                                {
                                    skipped_index.push_back(waypoint_variable_num + trajnum_ + 2 * (trajnum_ - 1) + trajid); // skipped angle index
                                }
                            }
                        }
                        else
                        {
                            if (singul_container[trajid] != last_singul)
                            {
                                trajid_container.push_back(trajid_list);
                                trajid_list.clear();
                                last_singul = singul_container[trajid];
                            }
                            trajid_list.push_back(trajid);
                        }

                        if (trajid == trajnum_ - 1)
                        {
                            trajid_container.push_back(trajid_list);
                        }

                        offset += 2 * (piece_num_container[trajid] - 1);
                        trajid++;
                    }

                    // If there are skipped trajectories
                    if (trajid_container.size() != trajnum_)
                    {
                        int trajnum_new = trajid_container.size();

                        std::vector<int> piece_num_container_new(trajnum_new);
                        Eigen::VectorXd traj_time_new(trajnum_new);
                        std::vector<Eigen::MatrixXd> iniState_container_new(trajnum_new);
                        std::vector<Eigen::MatrixXd> finState_container_new(trajnum_new);
                        std::vector<int> singul_container_new(trajnum_new);
                        std::vector<std::vector<Eigen::MatrixXd>> cfgHs_container_new(trajnum_new);
                        std::vector<std::vector<std::vector<Eigen::MatrixXd>>> aux_cfgHs_container_new(trajnum_new);
                        std::vector<int> waypoint_insert_index_container;

                        int variable_num_new_so_far = 0;
                        for (int seg = 0; seg < trajnum_new; seg++)
                        {
                            std::vector<int> &trajid_list = trajid_container[seg];
                            int piece_num = 0;
                            double traj_time = 0;
                            std::vector<Eigen::MatrixXd> cfgHs_new;
                            std::vector<std::vector<Eigen::MatrixXd>> aux_cfgHs_new(aux_vertices_list_.size());

                            for (int i = 0; i < trajid_list.size(); i++)
                            {
                                int trajid = trajid_list[i];
                                piece_num += piece_num_container[trajid];
                                traj_time += jerkOpt_container[trajid].getTraj(singul_container[trajid]).getTotalDuration();
                                cfgHs_new.insert(cfgHs_new.end(), cfgHs_container[trajid].begin(), cfgHs_container[trajid].end());
                                for (int m = 0; m < aux_vertices_list_.size(); m++)
                                {
                                    aux_cfgHs_new[m].insert(aux_cfgHs_new[m].end(), aux_cfgHs_container[trajid][m].begin(), aux_cfgHs_container[trajid][m].end());
                                }

                                if (trajid_list.size() > 1 && i != trajid_list.size() - 1)
                                {
                                    waypoint_insert_index_container.push_back(variable_num_new_so_far + 2 * (piece_num - 1));
                                }
                            }

                            piece_num_container_new[seg] = piece_num;
                            traj_time_new[seg] = traj_time * 0.9;
                            singul_container_new[seg] = singul_container[trajid_list[0]];
                            cfgHs_container_new[seg] = std::move(cfgHs_new);
                            aux_cfgHs_container_new[seg] = std::move(aux_cfgHs_new);

                            variable_num_new_so_far += 2 * (piece_num_container_new[seg]);
                        }

                        Eigen::VectorXd virtual_traj_time_new(trajnum_new);
                        RealT2VirtualT(traj_time_new, virtual_traj_time_new);

                        // Assign x_new
                        int waypoint_variable_num_new = 0;
                        for (const auto piece_num : piece_num_container_new)
                            waypoint_variable_num_new += 2 * (piece_num - 1);

                        Eigen::VectorXd x_new_waypoint(waypoint_variable_num_new);
                        Eigen::VectorXd x_new_time(trajnum_new);
                        Eigen::VectorXd x_new_gear(2 * (trajnum_new - 1));
                        Eigen::VectorXd x_new_angle(trajnum_new - 1);
                        Eigen::VectorXd x_new(waypoint_variable_num_new + trajnum_new + 3 * (trajnum_new - 1));

                        // New waypoint
                        int j = 0, k = 0;
                        for (int i = 0; i < waypoint_variable_num && j < waypoint_variable_num_new; i++)
                        {
                            if (k < waypoint_insert_index_container.size() && j == waypoint_insert_index_container[k])
                            {
                                x_new_waypoint.segment<2>(j) = new_merge_point[k];
                                j += 2;
                                k++;
                            }

                            if (std::find(skipped_index.begin(), skipped_index.end(), i) == skipped_index.end())
                            {
                                x_new_waypoint[j] = x[i];
                                j++;
                            }
                        }

                        // New time
                        j = 0;
                        for (int i = waypoint_variable_num; i < waypoint_variable_num + trajnum_ && j < trajnum_new; i++)
                        {
                            if (std::find(skipped_index.begin(), skipped_index.end(), i) == skipped_index.end())
                            {
                                x_new_time[j] = virtual_traj_time_new[j];
                                j++;
                            }
                        }

                        // New gear position
                        j = 0;
                        for (int i = waypoint_variable_num + trajnum_; i < waypoint_variable_num + trajnum_ + 2 * (trajnum_ - 1) && j < 2 * (trajnum_new - 1); i++)
                        {
                            if (std::find(skipped_index.begin(), skipped_index.end(), i) == skipped_index.end())
                            {
                                x_new_gear[j] = x[i];
                                j++;
                            }
                        }

                        // New angle
                        j = 0;
                        for (int i = waypoint_variable_num + trajnum_ + 2 * (trajnum_ - 1); i < variable_num && j < trajnum_new - 1; i++)
                        {
                            if (std::find(skipped_index.begin(), skipped_index.end(), i) == skipped_index.end())
                            {
                                x_new_angle[j] = x[i];
                                j++;
                            }
                        }

                        // Merge x_new
                        offset = 0;
                        x_new.segment(offset, x_new_waypoint.size()) = x_new_waypoint;
                        offset += x_new_waypoint.size();
                        x_new.segment(offset, x_new_time.size()) = x_new_time;
                        offset += x_new_time.size();
                        x_new.segment(offset, x_new_gear.size()) = x_new_gear;
                        offset += x_new_gear.size();
                        x_new.segment(offset, x_new_angle.size()) = x_new_angle;

                        // waypoints
                        offset = 0;
                        // std::cout << "waypoints: " << std::endl;
                        for (int trajid = 0; trajid < trajnum_new; trajid++)
                        {
                            // std::cout << "piece number: " << piece_num_container_new[trajid] << std::endl;
                            Eigen::Map<const Eigen::MatrixXd> P_new(x_new.data() + offset, 2, piece_num_container_new[trajid] - 1);
                            offset += 2 * (piece_num_container_new[trajid] - 1);
                            // std::cout << "(" << P_new.transpose() << ")," << std::endl;
                        }

                        // times
                        // std::cout << "times: " << std::endl;
                        Eigen::Map<const Eigen::VectorXd> t_new(x_new.data() + offset, trajnum_new);
                        Eigen::VectorXd T_new(trajnum_new);
                        VirtualT2RealT(t_new, T_new);
                        // std::cout << T_new.transpose() << std::endl;
                        offset += trajnum_new;

                        std::vector<Eigen::Map<const Eigen::MatrixXd>> Gear_new_container;
                        Eigen::Map<const Eigen::VectorXd> Angles_new(nullptr, 0);
                        if (gear_opt_)
                        {
                            // gear positions
                            // std::cout << "gear position: " << std::endl;
                            for (int trajid = 0; trajid < trajnum_new - 1; trajid++)
                            {
                                Eigen::Map<const Eigen::MatrixXd> Gear_new(x_new.data() + offset, 2, 1);
                                offset += 2;
                                Gear_new_container.push_back(Gear_new);
                                // std::cout << Gear_new.transpose() << std::endl;
                            }
                            // angles
                            // std::cout << "angles: " << std::endl;
                            new (&Angles_new) Eigen::Map<const Eigen::VectorXd>(x_new.data() + offset, trajnum_new - 1);
                            // std::cout << Angles_new.transpose() << std::endl;
                        }

                        // update containers
                        iniState_container_new[0] = iniState_container[0];
                        finState_container_new[trajnum_new - 1] = finState_container[trajnum_ - 1];
                        for (int trajid = 0; trajid < trajnum_new; trajid++)
                        {
                            Eigen::MatrixXd iniS(2, 3), finS(2, 3);
                            iniS.setZero();
                            finS.setZero();

                            if (trajid > 0)
                            {
                                double theta = Angles_new[trajid - 1];
                                iniS.col(0) = Gear_new_container[trajid - 1];
                                iniS.col(1) = Eigen::Vector2d(-non_siguav_ * cos(theta), -non_siguav_ * sin(theta));
                                iniState_container_new[trajid] = iniS;
                            }
                            if (trajid < trajnum_new - 1)
                            {
                                double theta = Angles_new[trajid];
                                finS.col(0) = Gear_new_container[trajid];
                                finS.col(1) = Eigen::Vector2d(non_siguav_ * cos(theta), non_siguav_ * sin(theta));
                                finState_container_new[trajid] = finS;
                            }
                        }
                        Angles_new.~Map();

                        jerkOpt_container.clear();
                        jerkOpt_container.resize(trajnum_new);
                        for (int i = 0; i < trajnum_new; i++)
                        {
                            jerkOpt_container[i].reset(piece_num_container_new[i]);
                        }

                        trajnum_ = trajnum_new;
                        piece_num_container = std::move(piece_num_container_new);
                        iniState_container = std::move(iniState_container_new);
                        finState_container = std::move(finState_container_new);
                        singul_container = std::move(singul_container_new);
                        cfgHs_container = std::move(cfgHs_container_new);
                        aux_cfgHs_container = std::move(aux_cfgHs_container_new);

                        lbfgs_params.delta = std::min(lbfgs_params.delta, 1.0e-10);
                        int result = lbfgs::lbfgs_optimize(
                            x_new,
                            final_cost,
                            PolyTrajOptimizer::costFunctionCallback,
                            NULL,
                            NULL,
                            this,
                            lbfgs_params);

                        /* ---------- get result and check collision ---------- */
                        if (result == lbfgs::LBFGS_CONVERGENCE ||
                            result == lbfgs::LBFGS_CANCELED ||
                            result == lbfgs::LBFGS_STOP ||
                            result == lbfgs::LBFGSERR_MAXIMUMITERATION ||
                            result == lbfgs::LBFGSERR_MAXIMUMLINESEARCH)
                        {
                            flag_success = true;
                            if (verbose)
                            {
                                if (result == lbfgs::LBFGSERR_MAXIMUMITERATION)
                                    printf("[Poly Traj Optimizer]: Lbfgs - Maximum number of iterations has been reached.\n");
                                if (result == lbfgs::LBFGSERR_MAXIMUMLINESEARCH)
                                    printf("[Poly Traj Optimizer]: Lbfgs - The line-search routine reaches the maximum number of evaluations.\n");

                                t2 = std::chrono::high_resolution_clock::now();
                                total_time_ms = t2 - t0;
                                printf("[Poly Traj Optimizer]: elapsed time: %5.3f ms, iter num: %d, cost: %5.3f\n", total_time_ms.count(), iter_num_, final_cost);
                            }
                        }
                        else
                        {
                            if (verbose)
                                printf("[Poly Traj Optimizer]: Solver error! Return = %d, %s. Skip this planning.\n", result, lbfgs::lbfgs_strerror(result));
                        }
                    }
                }
            }
        }
        else
        {
            if (verbose)
                printf("[Poly Traj Optimizer]: Solver error! Return = %d, %s. Skip this planning.\n", result, lbfgs::lbfgs_strerror(result));
        }

        return flag_success;
    }

    /* Callbacks by the L-BFGS optimizer */
    double PolyTrajOptimizer::costFunctionCallback(void *func_data, const Eigen::VectorXd &x, Eigen::VectorXd &grad)
    {
        double total_smooth_cost = 0, total_time_cost = 0;
        double total_obs_penalty_cost = 0, total_feas_penalty_cost = 0;
        PolyTrajOptimizer *opt = reinterpret_cast<PolyTrajOptimizer *>(func_data);
        int offset = 0;

        // Waypoints
        std::vector<Eigen::Map<const Eigen::MatrixXd>> P_container;
        std::vector<Eigen::Map<Eigen::MatrixXd>> gradP_container;
        for (int trajid = 0; trajid < opt->trajnum_; trajid++)
        {
            Eigen::Map<const Eigen::MatrixXd> P(x.data() + offset, 2, opt->piece_num_container[trajid] - 1);
            Eigen::Map<Eigen::MatrixXd> gradP(grad.data() + offset, 2, opt->piece_num_container[trajid] - 1);
            offset += 2 * (opt->piece_num_container[trajid] - 1);
            gradP.setZero();
            P_container.push_back(P);
            gradP_container.push_back(gradP);
        }

        // Times
        Eigen::Map<const Eigen::VectorXd> t(x.data() + offset, opt->trajnum_);
        Eigen::Map<Eigen::VectorXd> gradt(grad.data() + offset, opt->trajnum_);
        gradt.setZero();
        Eigen::VectorXd T(opt->trajnum_);
        opt->VirtualT2RealT(t, T);
        // Eigen::VectorXd gradT(opt->trajnum_);
        // gradT.setZero();

        std::vector<double> trajtimes;
        trajtimes.push_back(0.0);
        for (int trajid = 0; trajid < opt->trajnum_; trajid++)
        {
            // T(i) is sum time of i-segment traj
            trajtimes.push_back(T[trajid]);
        }
        offset += opt->trajnum_;

        // if (T.sum() > 1000 || T.sum() < 0.1)
        // {
        //     for (int trajid = 0; trajid < opt->trajnum_; trajid++)
        //     {
        //         gradP_container[trajid].setZero();
        //     }
        //     gradt.setZero();
        //     return 999999.0;
        // }

        std::vector<Eigen::Map<const Eigen::MatrixXd>> Gear_container;
        std::vector<Eigen::Map<Eigen::MatrixXd>> gradGear_container;
        Eigen::Map<const Eigen::VectorXd> Angles(nullptr, 0);
        Eigen::Map<Eigen::VectorXd> gradAngles(nullptr, 0);
        if (opt->gear_opt_)
        {
            // Gear Positions
            for (int trajid = 0; trajid < opt->trajnum_ - 1; trajid++)
            {
                Eigen::Map<const Eigen::MatrixXd> Gear(x.data() + offset, 2, 1);
                Eigen::Map<Eigen::MatrixXd> gradGear(grad.data() + offset, 2, 1);
                offset += 2;
                gradGear.setZero();
                Gear_container.push_back(Gear);
                gradGear_container.push_back(gradGear);
            }

            // Angles
            new (&Angles) Eigen::Map<const Eigen::VectorXd>(x.data() + offset, opt->trajnum_ - 1);
            new (&gradAngles) Eigen::Map<Eigen::VectorXd>(grad.data() + offset, opt->trajnum_ - 1);
            gradAngles.setZero();
        }

        // Calculate costs
        opt->constraint_pts_container.clear();
        for (int trajid = 0; trajid < opt->trajnum_; trajid++)
        {
            Eigen::VectorXd penalty_cost_array(2); // collision, feasibility
            penalty_cost_array.setZero();

            Eigen::MatrixXd IniS, FinS;
            IniS = opt->iniState_container[trajid];
            FinS = opt->finState_container[trajid];

            if (opt->gear_opt_)
            {
                if (trajid > 0)
                {
                    double theta = Angles[trajid - 1];
                    IniS.col(0) = Gear_container[trajid - 1];
                    IniS.col(1) = Eigen::Vector2d(-opt->non_siguav_ * cos(theta), -opt->non_siguav_ * sin(theta));
                }
                if (trajid < opt->trajnum_ - 1)
                {
                    double theta = Angles[trajid];
                    FinS.col(0) = Gear_container[trajid];
                    FinS.col(1) = Eigen::Vector2d(opt->non_siguav_ * cos(theta), opt->non_siguav_ * sin(theta));
                }
            }

            opt->jerkOpt_container[trajid].generate(P_container[trajid], T[trajid] / opt->piece_num_container[trajid], IniS, FinS);
            opt->jerkOpt_container[trajid].initSmGradCost();
            double smooth_cost = opt->jerkOpt_container[trajid].getTrajJerkCost(); // Smoothness cost
            opt->addPVAGradCost2CT(penalty_cost_array, trajid, trajtimes[trajid]); // Penalty cost

            total_smooth_cost += smooth_cost;
            total_obs_penalty_cost += penalty_cost_array(0);
            total_feas_penalty_cost += penalty_cost_array(1);
        }

        // Calculate gradients
        for (int trajid = 0; trajid < opt->trajnum_; trajid++)
        {
            opt->jerkOpt_container[trajid].calGrads_PT(); // gradC->gradT gradC->gradP
            gradP_container[trajid] = opt->jerkOpt_container[trajid].get_gdP();

            if (opt->gear_opt_)
            {
                // gdhead gdtail
                Eigen::Matrix<double, 2, 3> gradIni, gradFin;
                gradIni = opt->jerkOpt_container[trajid].get_gdHead();
                gradFin = opt->jerkOpt_container[trajid].get_gdTail();

                if (trajid > 0)
                {
                    double theta = Angles[trajid - 1];
                    gradGear_container[trajid - 1] += gradIni.col(0);
                    gradAngles[trajid - 1] += gradIni.col(1).transpose() * Eigen::Vector2d(opt->non_siguav_ * sin(theta), -opt->non_siguav_ * cos(theta));
                }
                if (trajid < opt->trajnum_ - 1)
                {
                    double theta = Angles[trajid];
                    gradGear_container[trajid] += gradFin.col(0);
                    gradAngles[trajid] += gradFin.col(1).transpose() * Eigen::Vector2d(-opt->non_siguav_ * sin(theta), opt->non_siguav_ * cos(theta));
                }
            }

            opt->VirtualTGradCost(T[trajid], t[trajid], opt->jerkOpt_container[trajid].get_gdT() / opt->piece_num_container[trajid] + opt->wei_time_, gradt[trajid]);
            double time_cost = T[trajid] * opt->wei_time_;
            total_time_cost += time_cost;
        }

        opt->iter_num_ += 1;

        // log cost history
        if (opt->log_cost_)
        {
            Eigen::Matrix<double, 5, 1> costs;
            costs << static_cast<double>(opt->iter_num_), total_smooth_cost, total_time_cost, total_obs_penalty_cost, total_feas_penalty_cost;
            opt->costs_history_container.push_back(costs);
        }

        if (opt->gear_opt_)
        {
            Angles.~Map();
            gradAngles.~Map();
        }

        return total_smooth_cost + total_time_cost + total_obs_penalty_cost + total_feas_penalty_cost;
    }

    void PolyTrajOptimizer::addPVAGradCost2CT(Eigen::VectorXd &costs, const int trajid, const double trajtime)
    {
        costs.setZero(); // collision, feasibility

        int singul = singul_container[trajid];
        std::vector<Eigen::MatrixXd> cfgHs = cfgHs_container[trajid];
        std::vector<std::vector<Eigen::MatrixXd>> aux_cfgHs_list = aux_cfgHs_container[trajid];

        int M = piece_num_container[trajid];

        Eigen::Vector2d sigma, dsigma, ddsigma, dddsigma, ddddsigma;
        double vel2, vel2_reci, vel3_reci, vel4_reci, vel5_reci, vel6_reci;
        double acc2, cur, omg2;
        // double latacc2;
        // double phi_dot;
        Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
        double t1, t2, t3, t4, t5;
        double z_h0, z_h1, z_h2, z_h3, z_h4;
        Eigen::Matrix2d ego_R;
        double step, alpha, omg;

        Eigen::Matrix<double, 6, 2> gradViolaPc, gradViolaVc, gradViolaAc, gradViolaKLc, gradViolaKRc, gradViolaOc;
        Eigen::Matrix<double, 6, 2> gradViolaLatAc, gradViolaPhidotLc, gradViolaPhidotRc;
        double gradViolaPt, gradViolaVt, gradViolaAt, gradViolaKLt, gradViolaKRt, gradViolaOt;
        // double gradViolaLatAt;
        // double gradViolaPhidotLt, gradViolaPhidotRt;
        double violaPos, violaVel, violaAcc, violaCurL, violaCurR, violaOmg;
        // double violaLatAcc;
        // double violaPhidotL, violaPhidotR;
        double violaPosPenaD, violaVelPenaD, violaAccPenaD, violaCurPenaDL, violaCurPenaDR, violaOmgPenaD;
        // double violaLatAccPenaD;
        // double violaPhidotPenaDL, violaPhidotPenaDR;
        double violaPosPena, violaVelPena, violaAccPena, violaCurPenaL, violaCurPenaR, violaOmgPena;
        // double violaLatAccPena;
        // double violaPhidotPenaL, violaPhidotPenaR;
        Eigen::Vector2d outerNormal;
        int pointid = -1;

        for (int i = 0; i < M; ++i)
        {
            int K = (i == 0 || i == M - 1) ? dense_traj_res_ : traj_res_;
            const Eigen::Matrix<double, 6, 2> &c = jerkOpt_container[trajid].getCoeffs().block<6, 2>(i * 6, 0);
            step = jerkOpt_container[trajid].getDt() / K; // delta_T_i / K
            t1 = 0.0;

            for (int j = 0; j <= K; ++j)
            {
                t2 = t1 * t1;
                t3 = t2 * t1;
                t4 = t3 * t1;
                t5 = t4 * t1;
                beta0 << 1.0, t1, t2, t3, t4, t5;
                beta1 << 0.0, 1.0, 2.0 * t1, 3.0 * t2, 4.0 * t3, 5.0 * t4;
                beta2 << 0.0, 0.0, 2.0, 6.0 * t1, 12.0 * t2, 20.0 * t3;
                beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * t1, 60.0 * t2;
                beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120 * t1;
                alpha = 1.0 / K * j;

                // update t1 for the next iteration
                t1 += step;
                pointid++;

                sigma = c.transpose() * beta0;
                dsigma = c.transpose() * beta1;
                ddsigma = c.transpose() * beta2;
                dddsigma = c.transpose() * beta3;
                ddddsigma = c.transpose() * beta4;

                omg = (j == 0 || j == K) ? 0.5 : 1.0;

                // some help values
                z_h0 = dsigma.norm();
                z_h1 = ddsigma.transpose() * dsigma;
                z_h2 = dddsigma.transpose() * dsigma;
                z_h3 = ddsigma.transpose() * B_h * dsigma;

                // save constraint points
                Eigen::Vector3d constraint_pt;
                constraint_pt << sigma[0], sigma[1], std::atan2(singul * dsigma[1], singul * dsigma[0]);
                constraint_pts_container.push_back(constraint_pt);

                // // avoid singularity
                // if (z_h0 < 1e-4 || (j == 0 && i == 0) || (i == M - 1 && j == K))
                // {
                //     continue;
                // }

                vel2_reci = 1.0 / (z_h0 * z_h0);
                vel3_reci = vel2_reci * (1.0 / z_h0);
                vel4_reci = vel2_reci * vel2_reci;
                vel5_reci = vel4_reci * (1.0 / z_h0);
                vel6_reci = vel2_reci * vel4_reci;

                double vel2_inv, vel3_inv, vel4_inv, vel5_inv, vel6_inv;
                if (z_h0 < epis_)
                {
                    double z_h0_epis = z_h0 + epis_;
                    double vel2_reci_e, vel3_reci_e, vel4_reci_e, vel5_reci_e, vel6_reci_e;
                    vel2_reci_e = 1.0 / (z_h0_epis * z_h0_epis);
                    vel3_reci_e = 1.0 / (z_h0_epis * z_h0_epis * z_h0_epis);
                    vel4_reci_e = 1.0 / (z_h0_epis * z_h0_epis * z_h0_epis * z_h0_epis);
                    vel5_reci_e = 1.0 / (z_h0_epis * z_h0_epis * z_h0_epis * z_h0_epis * z_h0_epis);
                    vel6_reci_e = 1.0 / (z_h0_epis * z_h0_epis * z_h0_epis * z_h0_epis * z_h0_epis * z_h0_epis);

                    vel2_inv = vel2_reci_e;
                    vel3_inv = vel3_reci_e;
                    vel4_inv = vel4_reci_e;
                    vel5_inv = vel5_reci_e;
                    vel6_inv = vel6_reci_e;
                }
                else
                {
                    vel2_inv = vel2_reci;
                    vel3_inv = vel3_reci;
                    vel4_inv = vel4_reci;
                    vel5_inv = vel5_reci;
                    vel6_inv = vel6_reci;
                }

                /* feasibility with velocity */
                vel2 = z_h0 * z_h0;
                violaVel = vel2 - max_vel_ * max_vel_;

                if (violaVel > 0.0)
                {
                    positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);
                    gradViolaVc = 2.0 * beta1 * dsigma.transpose(); // 6*2
                    gradViolaVt = 2.0 * alpha * z_h1;               // 1*1

                    jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaVelPenaD * gradViolaVc;
                    jerkOpt_container[trajid].get_gdT() += omg * wei_feas_ * (violaVelPenaD * gradViolaVt * step + violaVelPena / K);
                    costs(1) += omg * step * wei_feas_ * violaVelPena;
                }

                /* feasibility with acceleration */
                z_h4 = z_h1 * vel2_inv;
                acc2 = z_h1 * z_h1 * vel2_inv;
                violaAcc = acc2 - max_acc_ * max_acc_;

                if (violaAcc > 0.0)
                {
                    positiveSmoothedL1(violaAcc, violaAccPena, violaAccPenaD);
                    gradViolaAc = 2.0 * beta1 * (z_h4 * ddsigma - z_h4 * z_h4 * dsigma).transpose() +
                                  2.0 * beta2 * (z_h4 * dsigma).transpose();                         // 6*2
                    gradViolaAt = 2.0 * alpha * z_h4 * (ddsigma.squaredNorm() - z_h1 * z_h4 + z_h2); // 1*1

                    jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaAccPenaD * gradViolaAc;
                    jerkOpt_container[trajid].get_gdT() += omg * wei_feas_ * (violaAccPenaD * gradViolaAt * step + violaAccPena / K);
                    costs(1) += omg * step * wei_feas_ * violaAccPena;
                }

                // // feasibility with lateral acceleration
                // latacc2 = z_h3 * z_h3 * vel2_inv;
                // violaLatAcc = latacc2 - max_latacc_ * max_latacc_;

                // if (violaLatAcc > 0.0)
                // {
                //     positiveSmoothedL1(violaLatAcc, violaLatAccPena, violaLatAccPenaD);
                //     gradViolaLatAc = 2.0 * beta1 * (z_h3 * vel2_inv * B_h.transpose() * ddsigma - z_h3 * z_h3 * vel2_inv * vel2_inv * dsigma).transpose() +
                //                      2.0 * beta2 * (z_h3 * vel2_inv * B_h * dsigma).transpose();                                                                                                                                           // 6*2
                //     gradViolaLatAt = 2.0 * alpha * (z_h3 * vel2_inv * (ddsigma.transpose() * B_h.transpose() * ddsigma)(0, 0) - z_h3 * z_h3 * vel2_inv * vel2_inv * z_h1 + z_h3 * vel2_inv * (dddsigma.transpose() * B_h * dsigma)(0, 0)); // 1*1

                //     jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaLatAccPenaD * gradViolaLatAc;
                //     jerkOpt_container[trajid].get_gdT() += omg * wei_feas_ * (violaLatAccPenaD * gradViolaLatAt * step + violaLatAccPena / K);
                //     costs(1) += omg * step * wei_feas_ * violaLatAccPena;
                // }

                /* feasibility with curvature */
                cur = singul * z_h3 * vel3_inv;
                violaCurL = cur - max_cur_;
                violaCurR = -cur - max_cur_;
                double factor = 50;

                if (violaCurL > 0.0)
                {
                    positiveSmoothedL1(violaCurL, violaCurPenaL, violaCurPenaDL);
                    gradViolaKLc = singul * beta1 * (vel3_inv * B_h.transpose() * ddsigma - 3 * z_h3 * vel5_inv * dsigma).transpose() +
                                   singul * beta2 * (vel3_inv * B_h * dsigma).transpose();                                                                                           // 6*2
                    gradViolaKLt = alpha * singul * vel3_inv * (ddsigma.transpose() * B_h.transpose() * ddsigma - 3 * z_h1 * z_h3 * vel2_inv + dddsigma.transpose() * B_h * dsigma); // 1*1

                    jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * (wei_feas_ * factor) * violaCurPenaDL * gradViolaKLc;
                    jerkOpt_container[trajid].get_gdT() += omg * (wei_feas_ * factor) * (violaCurPenaDL * gradViolaKLt * step + violaCurPenaL / K);
                    costs(1) += omg * step * (wei_feas_ * factor) * violaCurPenaL;
                }

                if (violaCurR > 0.0)
                {
                    positiveSmoothedL1(violaCurR, violaCurPenaR, violaCurPenaDR);
                    gradViolaKRc = -(singul * beta1 * (vel3_inv * B_h.transpose() * ddsigma - 3 * z_h3 * vel5_inv * dsigma).transpose() +
                                     singul * beta2 * (vel3_inv * B_h * dsigma).transpose());                                                                                         // 6*2
                    gradViolaKRt = -alpha * singul * vel3_inv * (ddsigma.transpose() * B_h.transpose() * ddsigma - 3 * z_h1 * z_h3 * vel2_inv + dddsigma.transpose() * B_h * dsigma); // 1*1

                    jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * (wei_feas_ * factor) * violaCurPenaDR * gradViolaKRc;
                    jerkOpt_container[trajid].get_gdT() += omg * (wei_feas_ * factor) * (violaCurPenaDR * gradViolaKRt * step + violaCurPenaR / K);
                    costs(1) += omg * step * (wei_feas_ * factor) * violaCurPenaR;
                }

                /* feasibility with omega */
                omg2 = z_h3 * z_h3 * vel2_inv * vel2_inv;
                violaOmg = omg2 - max_omg_ * max_omg_;

                if (violaOmg > 0.0)
                {
                    positiveSmoothedL1(violaOmg, violaOmgPena, violaOmgPenaD);
                    gradViolaOc = 2.0 * beta1 * (z_h3 * vel4_inv * B_h.transpose() * ddsigma - 2.0 * z_h3 * z_h3 * vel6_inv * dsigma).transpose() +
                                  2.0 * beta2 * (z_h3 * vel4_inv * B_h * dsigma).transpose();                                                                                                                      // 6*2
                    gradViolaOt = 2.0 * alpha * (z_h3 * vel4_inv * ddsigma.transpose() * B_h.transpose() * ddsigma - 2.0 * z_h3 * z_h3 * vel6_inv * z_h1 + z_h3 * vel4_inv * dddsigma.transpose() * B_h * dsigma); // 1*1

                    jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_feas_ * violaOmgPenaD * gradViolaOc;
                    jerkOpt_container[trajid].get_gdT() += omg * wei_feas_ * (violaOmgPenaD * gradViolaOt * step + violaOmgPena / K);
                    costs(1) += omg * step * wei_feas_ * violaOmgPena;
                }

                // /* feasibility with steer rate */
                // double L = car_wheelbase_;
                // double n1 = z_h0;
                // double n2 = n1 * n1;
                // double n3 = n2 * n1;
                // double n4 = n3 * n1;
                // double n5 = n4 * n1;
                // double n6 = n5 * n1;

                // double z_h5 = dddsigma.transpose() * B_h * dsigma;
                // double z_h6 = ddsigma.transpose() * B_h * ddsigma;

                // double phidot_numerator = singul * L * ((z_h5 + z_h6) * n3 - 3 * z_h3 * z_h1 * n1);
                // double phidot_denominator = n6 + L * L * z_h3 * z_h3;
                // phi_dot = phidot_numerator / (phidot_denominator + epis_);
                // violaPhidotL = phi_dot - max_phi_dot_;
                // violaPhidotR = -phi_dot - max_phi_dot_;
                // factor = 1.0;

                // if (violaPhidotL > 0.0)
                // {
                //     positiveSmoothedL1(violaPhidotL, violaPhidotPenaL, violaPhidotPenaDL);
                //     Eigen::Vector2d partial_N_over_partial_dsigma = singul * L * (3 * z_h5 * n1 * dsigma + n3 * B_h.transpose() * dddsigma + 3 * z_h6 * n1 * dsigma - 3 * B_h.transpose() * ddsigma * z_h1 * n1 - 3 * z_h3 * ddsigma * n1 - 3 * z_h3 * z_h1 * dsigma / n1);
                //     Eigen::Vector2d partial_D_over_partial_dsigma = 6 * n4 * dsigma + 2 * L * L * z_h3 * B_h.transpose() * ddsigma;
                //     Eigen::Vector2d partial_N_over_partial_ddsigma = singul * L * (2 * B_h * ddsigma * n3 - 3 * B_h * dsigma * z_h1 * n1 - 3 * z_h3 * dsigma * n1);
                //     Eigen::Vector2d partial_D_over_partial_ddsigma = 2 * L * L * z_h3 * B_h * dsigma;
                //     Eigen::Vector2d partial_N_over_partial_dddsigma = singul * L * B_h * dsigma * n3;
                //     Eigen::Vector2d partial_D_over_partial_dddsigma = Eigen::Vector2d::Zero();

                //     Eigen::Vector2d partial_phi_dot_over_partial_dsigma = (partial_N_over_partial_dsigma * phidot_denominator - partial_D_over_partial_dsigma * phidot_numerator) / pow(phidot_denominator, 2);
                //     Eigen::Vector2d partial_phi_dot_over_partial_ddsigma = (partial_N_over_partial_ddsigma * phidot_denominator - partial_D_over_partial_ddsigma * phidot_numerator) / pow(phidot_denominator, 2);
                //     Eigen::Vector2d partial_phi_dot_over_partial_dddsigma = (partial_N_over_partial_dddsigma * phidot_denominator - partial_D_over_partial_dddsigma * phidot_numerator) / pow(phidot_denominator, 2);

                //     gradViolaPhidotLc = beta1 * partial_phi_dot_over_partial_dsigma.transpose() + beta2 * partial_phi_dot_over_partial_ddsigma.transpose() + beta3 * partial_phi_dot_over_partial_dddsigma.transpose();
                //     gradViolaPhidotLt = alpha * (ddsigma.transpose() * partial_phi_dot_over_partial_dsigma + dddsigma.transpose() * partial_phi_dot_over_partial_ddsigma + ddddsigma.transpose() * partial_phi_dot_over_partial_dddsigma)(0,0);

                //     jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * (wei_feas_ * factor) * violaPhidotPenaDL * gradViolaPhidotLc;
                //     jerkOpt_container[trajid].get_gdT() += omg * (wei_feas_ * factor) * (violaPhidotPenaDL * gradViolaPhidotLt * step + violaPhidotPenaL / K);
                //     costs(1) += omg * step * (wei_feas_ * factor) * violaPhidotPenaL;
                // }

                // if (violaPhidotR > 0.0)
                // {
                //     positiveSmoothedL1(violaPhidotR, violaPhidotPenaR, violaPhidotPenaDR);
                //     Eigen::Vector2d partial_N_over_partial_dsigma = singul * L * (3 * z_h5 * n1 * dsigma + n3 * B_h.transpose() * dddsigma + 3 * z_h6 * n1 * dsigma - 3 * B_h.transpose() * ddsigma * z_h1 * n1 - 3 * z_h3 * ddsigma * n1 - 3 * z_h3 * z_h1 * dsigma / n1);
                //     Eigen::Vector2d partial_D_over_partial_dsigma = 6 * n4 * dsigma + 2 * L * L * z_h3 * B_h.transpose() * ddsigma;
                //     Eigen::Vector2d partial_N_over_partial_ddsigma = singul * L * (2 * B_h * ddsigma * n3 - 3 * B_h * dsigma * z_h1 * n1 - 3 * z_h3 * dsigma * n1);
                //     Eigen::Vector2d partial_D_over_partial_ddsigma = 2 * L * L * z_h3 * B_h * dsigma;
                //     Eigen::Vector2d partial_N_over_partial_dddsigma = singul * L * B_h * dsigma * n3;
                //     Eigen::Vector2d partial_D_over_partial_dddsigma = Eigen::Vector2d::Zero();

                //     Eigen::Vector2d partial_phi_dot_over_partial_dsigma = (partial_N_over_partial_dsigma * phidot_denominator - partial_D_over_partial_dsigma * phidot_numerator) / pow(phidot_denominator, 2);
                //     Eigen::Vector2d partial_phi_dot_over_partial_ddsigma = (partial_N_over_partial_ddsigma * phidot_denominator - partial_D_over_partial_ddsigma * phidot_numerator) / pow(phidot_denominator, 2);
                //     Eigen::Vector2d partial_phi_dot_over_partial_dddsigma = (partial_N_over_partial_dddsigma * phidot_denominator - partial_D_over_partial_dddsigma * phidot_numerator) / pow(phidot_denominator, 2);

                //     gradViolaPhidotRc = -(beta1 * partial_phi_dot_over_partial_dsigma.transpose() + beta2 * partial_phi_dot_over_partial_ddsigma.transpose() + beta3 * partial_phi_dot_over_partial_dddsigma.transpose());
                //     gradViolaPhidotRt = -alpha * (ddsigma.transpose() * partial_phi_dot_over_partial_dsigma + dddsigma.transpose() * partial_phi_dot_over_partial_ddsigma + ddddsigma.transpose() * partial_phi_dot_over_partial_dddsigma)(0,0);

                //     jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * (wei_feas_ * factor) * violaPhidotPenaDR * gradViolaPhidotRc;
                //     jerkOpt_container[trajid].get_gdT() += omg * (wei_feas_ * factor) * (violaPhidotPenaDR * gradViolaPhidotRt * step + violaPhidotPenaR / K);
                //     costs(1) += omg * step * (wei_feas_ * factor) * violaPhidotPenaR;
                // }

                /* static obstacle collision */
                ego_R << dsigma(0), -dsigma(1),
                    dsigma(1), dsigma(0);
                ego_R = singul * ego_R / z_h0;

                // iterate over ego vertices
                int corr_k = cfgHs[pointid].cols();
                
                for (const auto &le : ego_vertices_)
                {
                    Eigen::Vector2d bpt = sigma + ego_R * le;
                    Eigen::Matrix2d temp_l_Bl;
                    temp_l_Bl << le(0), -le(1),
                        le(1), le(0);
                    Eigen::Matrix2d F_le = singul * temp_l_Bl.transpose() / z_h0 - vel2_inv * dsigma * (ego_R * le).transpose();

                    for (int k = 0; k < corr_k; k++)
                    {
                        outerNormal = cfgHs[pointid].col(k).head<2>();
                        violaPos = outerNormal.dot(bpt - cfgHs[pointid].col(k).tail<2>());

                        if (violaPos > 0)
                        {
                            positiveSmoothedL1(violaPos, violaPosPena, violaPosPenaD);
                            gradViolaPc = beta0 * outerNormal.transpose() +
                                          beta1 * (F_le * outerNormal).transpose();                                                      // 6*2
                            gradViolaPt = alpha * (dsigma.transpose() * outerNormal + ddsigma.transpose() * (F_le * outerNormal))(0, 0); // 1*1

                            jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_obs_ * violaPosPenaD * gradViolaPc;
                            jerkOpt_container[trajid].get_gdT() += omg * wei_obs_ * (violaPosPenaD * gradViolaPt * step + violaPosPena / K);
                            costs(0) += omg * step * wei_obs_ * violaPosPena; // cost is the same
                        }
                    }
                }

                // iterate over aux vertices
                for (int m = 0; m < aux_vertices_list_.size(); m++)
                {
                    int corr_k = aux_cfgHs_list[m][pointid].cols();
                    
                    for (const auto &le : aux_vertices_list_[m])
                    {
                        Eigen::Vector2d bpt = sigma + ego_R * le;
                        Eigen::Matrix2d temp_l_Bl;
                        temp_l_Bl << le(0), -le(1),
                            le(1), le(0);
                        Eigen::Matrix2d F_le = singul * temp_l_Bl.transpose() / z_h0 - vel2_inv * dsigma * (ego_R * le).transpose();

                        for (int k = 0; k < corr_k; k++)
                        {
                            outerNormal = aux_cfgHs_list[m][pointid].col(k).head<2>();
                            violaPos = outerNormal.dot(bpt - aux_cfgHs_list[m][pointid].col(k).tail<2>());

                            if (violaPos > 0)
                            {
                                positiveSmoothedL1(violaPos, violaPosPena, violaPosPenaD);
                                gradViolaPc = beta0 * outerNormal.transpose() +
                                              beta1 * (F_le * outerNormal).transpose();                                                      // 6*2
                                gradViolaPt = alpha * (dsigma.transpose() * outerNormal + ddsigma.transpose() * (F_le * outerNormal))(0, 0); // 1*1

                                jerkOpt_container[trajid].get_gdC().block<6, 2>(i * 6, 0) += omg * step * wei_obs_ * violaPosPenaD * gradViolaPc;
                                jerkOpt_container[trajid].get_gdT() += omg * wei_obs_ * (violaPosPenaD * gradViolaPt * step + violaPosPena / K);
                                costs(0) += omg * step * wei_obs_ * violaPosPena; // cost is the same
                            }
                        }
                    }
                }
            }
        }
    }

} // namespace plan_manage