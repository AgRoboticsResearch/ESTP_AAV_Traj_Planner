#ifndef _TRAJ_OPTIMIZER_HPP_
#define _TRAJ_OPTIMIZER_HPP_

#include <Eigen/Eigen>
#include <yaml-cpp/yaml.h>
#include <chrono>
#include <cmath>

#include "vehicle_utils/vehicle_geometry.hpp"
#include "plan_utils/traj_container.hpp"
#include "optimization/lbfgs.hpp"

namespace plan_manage
{

    using namespace std;
    using namespace vehicle_utils;

    class PolyTrajOptimizer
    {

    private:
        /* trajectory */
        int traj_res_;       // number of distinctive constraint points each piece
        int dense_traj_res_; // number of distinctive constraint points of the first and last piece (should be more dense)!
        int trajnum_;        // number of trajectory pieces

        /* optimization parameters */
        double wei_obs_;                     // obstacle weight
        double wei_feas_;                    // feasibility weight
        double wei_time_;                    // time weight
        double max_vel_, max_acc_, max_cur_; // dynamic limits
        double max_phi_dot_, max_latacc_;
        double max_omg_;
        double non_siguav_; // small positive velocity
        double epis_;       // helper epis

        enum FORCE_STOP_OPTIMIZE_TYPE
        {
            DONT_STOP,
            STOP_FOR_REBOUND,
            STOP_FOR_ERROR
        } force_stop_type_;

        int iter_num_; // iteration of the solver
        Eigen::Matrix<double, 2, 2> B_h;

        /* L-BFGS parameters */
        int lbfgs_mem_size;
        int lbfgs_past;
        double lbfgs_delta;
        double lbfgs_g_epsilon;
        int lbfgs_max_linesearch;
        int lbfgs_max_iterations;
        double lbfgs_f_dec_coeff;
        double lbfgs_s_curv_coeff;
        double lbfgs_cautious_factor;

        /* vehicle */
        double car_wheelbase_;
        std::vector<Eigen::Vector2d> ego_vertices_;
        std::vector<std::vector<Eigen::Vector2d>> aux_vertices_list_;

        /* containers */
        std::vector<Eigen::MatrixXd> iniState_container;
        std::vector<Eigen::MatrixXd> finState_container;
        std::vector<int> piece_num_container;
        std::vector<Eigen::Vector3d> constraint_pts_container;
        std::vector<int> singul_container;
        // Each col of cfgHs denotes a facet (outter_normal^T,point^T)^T
        std::vector<std::vector<Eigen::MatrixXd>> cfgHs_container;
        std::vector<std::vector<std::vector<Eigen::MatrixXd>>> aux_cfgHs_container;
        
        std::vector<plan_utils::MinJerkOpt> jerkOpt_container;
        std::vector<Eigen::VectorXd> costs_history_container;

        /* flags */
        bool gear_opt_;    // gear switch optimization flag
        bool remove_gear_; // remove gear switch if possible
        bool log_cost_;    // log cost history flag

    public:
        PolyTrajOptimizer() {};
        PolyTrajOptimizer(const std::string &config_string, VehicleGeometry::Ptr &vehicle_geo_ptr);
        ~PolyTrajOptimizer() {};
        typedef std::shared_ptr<PolyTrajOptimizer> Ptr;

        /* helper functions */
        inline const std::vector<plan_utils::MinJerkOpt> *getMinJerkOptPtr(void) { return &jerkOpt_container; }
        inline int get_traj_res_() { return traj_res_; };
        inline int get_dense_traj_res_() { return dense_traj_res_; };
        const int getTrajNum() const { return trajnum_; }
        const std::vector<int> getSinguls() const { return singul_container; }
        const std::vector<Eigen::Vector3d> getConstraintPts() const { return constraint_pts_container; }
        const std::vector<Eigen::VectorXd> getCostsHistory() const { return costs_history_container; }

        /* main planning API */
        bool OptimizeTrajectory(const std::vector<Eigen::MatrixXd> &iniStates, const std::vector<Eigen::MatrixXd> &finStates,
                                const std::vector<Eigen::MatrixXd> &initInnerPts, const Eigen::VectorXd &initTs,
                                const std::vector<std::vector<Eigen::MatrixXd>> &hPoly_container,
                                const std::vector<std::vector<std::vector<Eigen::MatrixXd>>> &aux_hPoly_container,
                                const std::vector<int> &singuls,
                                const bool verbose = true);

        // double log_sum_exp(double alpha, Eigen::VectorXd &all_dists, double &exp_sum);

    private:
        /* callbacks by the L-BFGS optimizer */
        static double costFunctionCallback(void *func_data, const Eigen::VectorXd &x, Eigen::VectorXd &grad);

        static int earlyExitCallback(void *func_data, const double *x, const double *g,
                                     const double fx, const double xnorm, const double gnorm,
                                     const double step, int n, int k, int ls);

        /* mappings between real world time and unconstrained virtual time */
        template <typename EIGENVEC>
        void RealT2VirtualT(const Eigen::VectorXd &RT, EIGENVEC &VT);
        void RealT2VirtualT(double &RT, double &VT);

        template <typename EIGENVEC>
        void VirtualT2RealT(const EIGENVEC &VT, Eigen::VectorXd &RT);
        void VirtualT2RealT(double &VT, double &RT);

        template <typename EIGENVEC, typename EIGENVECGD>
        void VirtualTGradCost(const Eigen::VectorXd &RT, const EIGENVEC &VT,
                              const Eigen::VectorXd &gdRT, EIGENVECGD &gdVT);
        void VirtualTGradCost(const double &RT, const double &VT,
                              const double &gdRT, double &gdVT);

        /* gradient and cost evaluation functions */
        void addPVAGradCost2CT(Eigen::VectorXd &costs, const int trajid, const double trajtime);

        void positiveSmoothedL1(const double &x, double &f, double &df);
        void positiveSmoothedL3(const double &x, double &f, double &df);
    };

} // namespace plan_manage

#endif // _TRAJ_OPTIMIZER_HPP_