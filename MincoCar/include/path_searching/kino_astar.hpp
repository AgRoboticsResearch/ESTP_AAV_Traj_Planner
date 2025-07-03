#ifndef _KINODYNAMIC_ASTAR_HPP_
#define _KINODYNAMIC_ASTAR_HPP_

#include <Eigen/Eigen>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <queue>
#include <math.h>
#include <utility>
#include <chrono>
#include <yaml-cpp/yaml.h>

#include "circle_cover/circle_cover.hpp"
#include "map_utils/voxel_map.hpp"
#include "map_utils/raycast.hpp"
#include "vehicle_utils/vehicle_geometry.hpp"
#include "plan_utils/traj_container.hpp"
#include <ompl/base/spaces/ReedsSheppStateSpace.h>
#include <ompl/base/spaces/DubinsStateSpace.h>
#include <ompl/geometric/SimpleSetup.h>

namespace path_searching
{
  using namespace std;
  using namespace map_utils;
  using namespace vehicle_utils;

  constexpr char IN_CLOSE_SET = 'a';
  constexpr char IN_OPEN_SET = 'b';
  constexpr char NOT_EXPAND = 'c';
  constexpr double inf = 1 >> 30;

  class PathNode
  {
  public:
    /* -------------------- */
    Eigen::Vector2i index;
    int yaw_idx;
    /* --- the state is x y theta(orientation) */
    Eigen::Vector3d state;
    double g_score, f_score;
    /* control input should be steer and arc */
    Eigen::Vector2d input;
    std::shared_ptr<PathNode> parent;
    char node_state;
    int singul = 0;
    int number;
    /* -------------------- */
    PathNode()
    {
      parent = NULL;
      node_state = NOT_EXPAND;
    }
    ~PathNode(){};
  };
  typedef std::shared_ptr<PathNode> PathNodePtr;

  class MiddleNode
  {
  public:
    double start_s;
    double start_t;
    double start_vel;
    double end_s;
    double end_t;
    double end_vel;
    Eigen::Vector3i s_t_v_idx;
    double length;
    double distance2end;
    double acc;
    double g_score, f_score;
    char node_state;

    std::shared_ptr<MiddleNode> parent;

    MiddleNode()
    {
      parent = NULL;
      node_state = NOT_EXPAND;
    }
    ~MiddleNode(){};
  };
  typedef std::shared_ptr<MiddleNode> MiddleNodePtr;

  class NodeComparator
  {
  public:
    template <class NodePtr>
    bool operator()(NodePtr node1, NodePtr node2)
    {
      return node1->f_score > node2->f_score;
    }
  };

  template <typename T>
  struct matrix_hash : std::unary_function<T, size_t>
  {
    std::size_t operator()(T const &matrix) const
    {
      size_t seed = 0;
      for (long int i = 0; i < matrix.size(); ++i)
      {
        auto elem = *(matrix.data() + i);
        seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  template <class NodePtr>
  class NodeHashTable
  {
  private:
    /* data */

    std::unordered_map<Eigen::Vector2i, NodePtr, matrix_hash<Eigen::Vector2i>> data_2d_;
    std::unordered_map<Eigen::Vector3i, NodePtr, matrix_hash<Eigen::Vector3i>> data_3d_;

  public:
    NodeHashTable(/* args */)
    {
    }
    ~NodeHashTable()
    {
    }
    // : for 2d vehicle planning
    void insert(Eigen::Vector2i idx, NodePtr node)
    {
      data_2d_.insert(std::make_pair(idx, node));
    }
    // for 3d vehicle planning
    void insert(Eigen::Vector2i idx, int yaw_idx, NodePtr node)
    {
      data_3d_.insert(std::make_pair(Eigen::Vector3i(idx(0), idx(1), yaw_idx), node));
    }
    void insert(Eigen::Vector3i idx, NodePtr node)
    {
      data_3d_.insert(std::make_pair(idx, node));
    }

    NodePtr find(Eigen::Vector2i idx)
    {
      auto iter = data_2d_.find(idx);
      return iter == data_2d_.end() ? NULL : iter->second;
    }
    NodePtr find(Eigen::Vector2i idx, int yaw_idx)
    {
      auto iter = data_3d_.find(Eigen::Vector3i(idx(0), idx(1), yaw_idx));
      return iter == data_3d_.end() ? NULL : iter->second;
    }

    void clear()
    {
      data_2d_.clear();
      data_3d_.clear();
    }
  };

  template <class MiddleNode>
  class MiddleNodeHashTable
  {
  private:
    /* data */
    std::unordered_map<Eigen::Vector3i, MiddleNodePtr, matrix_hash<Eigen::Vector3i>> middle_data_3d_;

  public:
    MiddleNodeHashTable(/* args */)
    {
    }
    ~MiddleNodeHashTable()
    {
    }

    // For middle node searching
    void insert(Eigen::Vector3i idx, MiddleNodePtr node)
    {
      middle_data_3d_.insert(std::make_pair(idx, node));
    }

    MiddleNodePtr find(Eigen::Vector3i idx)
    {
      auto iter = middle_data_3d_.find(idx);
      return iter == middle_data_3d_.end() ? NULL : iter->second;
    }

    void clear()
    {
      middle_data_3d_.clear();
    }
  };

  class KinoAstar
  {

  private:
    /* ---------- main data structure ---------- */
    std::vector<PathNodePtr> path_node_pool_;
    int use_node_num_, iter_num_;
    NodeHashTable<PathNodePtr> expanded_nodes_;
    std::priority_queue<PathNodePtr, std::vector<PathNodePtr>, NodeComparator> open_set_;
    std::vector<PathNodePtr> path_nodes_;
    std::vector<Eigen::Vector2d> last_path_pos_, last_path_pos_temp_;

    /*--------------speed planning-----------------*/
    MiddleNodeHashTable<MiddleNodePtr> expanded_middle_nodes_;
    double calculateMiniTime2singulPoint(double initspeed, double distance2end);
    int s2idx(double s);
    int t2idx(double t);
    int v2idx(double v);

    std::vector<MiddleNodePtr> middle_node_pool_;
    std::vector<MiddleNodePtr> middle_nodes_;
    int use_time_node_num_;
    double total_s_;
    std::multimap<double, MiddleNodePtr> openSet_;
    double gettimeHeu(int &idx, double &vel);
    double calculateCurNode2nearestSingul(int idx);

    /* ---------- record data ---------- */
    Eigen::Vector4d start_state_, end_state_;
    Eigen::Vector2d start_ctrl_;
    bool is_shot_succ_ = false;
    bool has_path_ = false;
    /* ---------- parameter ---------- */
    /* search */
    double step_arc_ = 1.0;
    double step_arc_res_ = 0.5;
    double step_steer_res_ = 0.5;
    double max_vel_ = 10.0;
    double max_acc_ = 3.0;
    double min_vel_ = -1.0;
    double min_acc_ = -1.0;
    double max_cur_ = 0.35;
    double max_steer_ = 0.78539815;
    double max_seach_time = 0.1;
    double traj_forward_penalty = 1.0;
    double traj_back_penalty = 1.0;
    double traj_gear_switch_penalty = 10.0;
    double traj_steer_penalty = 0.0;
    double traj_steer_change_penalty = 0.0;
    double max_horizon_;
    double end_horizon_;
    double lambda_heu_;
    int collision_check_type_;

    double time_resolution_ = 0.1;
    double distance_resolution_ = 0.5;
    double velocity_resolution_ = 0.5;

    // double margin_;
    int allocate_num_;
    int check_num_;
    double tie_breaker_ = 1.0 + 1.0 / 10000;
    int nearest_idx_;
    double nearest_init_distance_;

    /* vehicle parameters*/
    double car_wheelbase_;
    double max_sagitta;

    /* map */
    double resolution_, inv_resolution_, yaw_resolution_, inv_yaw_resolution_;
    Eigen::Vector2d origin_, map_size_3d_;
    double yaw_origin_ = -M_PI;
    /*shot hzc*/
    std::vector<double> shot_timeList;
    std::vector<double> shot_lengthList;
    std::vector<int> shotindex;
    std::vector<int> shot_SList;
    std::vector<Eigen::Vector3d> SampleTraj;
    std::vector<Eigen::Vector2d> positionList_;
    std::vector<int> direction_change_idx_;

    /* helper */
    Eigen::Vector2i posToIndex(Eigen::Vector2d pt);
    int yawToIndex(double &yaw);
    double normalize_angle(const double &theta);
    void retrievePath(PathNodePtr end_node);

    /* shot trajectory */
    vector<double> cubic(double a, double b, double c, double d);
    vector<double> quartic(double a, double b, double c, double d, double e);
    double estimateHeuristic(Eigen::VectorXd x1, Eigen::VectorXd x2, double &optimal_time);

    void getFlatState(Eigen::Vector4d state, Eigen::Vector2d control_input,
                      Eigen::MatrixXd &flat_state, int singul);
    std::vector<Eigen::Vector2d> getFlatInputs();

    bool is_shot_sucess(Eigen::Vector3d state1, Eigen::Vector3d state2, int &singual);
    bool is_shot_sucess_use_circle(Eigen::Vector3d state1, Eigen::Vector3d state2, int &singual);
    double computeShotTraj(Eigen::Vector3d &state1, Eigen::Vector3d &state2,
                           std::vector<Eigen::Vector3d> &path_list,
                           double &len);
    void stateTransit(Eigen::Vector3d &state0, Eigen::Vector3d &state1,
                      Eigen::Vector2d &ctrl_input);

    double evaluateLength(double curt, double locallength, double localtime, double startV = 0.0, double endV = 0.0);
    double evaluateDuration(double length, double startV = 0.0, double endV = 0.0);
    void resetEndState(Eigen::Vector4d &start_state, Eigen::Vector4d &end_state);

    double non_siguav;
    ompl::base::StateSpacePtr shotptr;

    inline int getSingularity(double vel)
    {
      int singul = 0;
      if (fabs(vel) > 1e-2)
      {
        if (vel >= 0.0)
        {
          singul = 1;
        }
        else
        {
          singul = -1;
        }
      }

      return singul;
    }
    // getHeu
    inline double getHeu(Eigen::VectorXd x1, Eigen::VectorXd x2)
    {
      double dx = abs(x1(0) - x2(0));
      double dy = abs(x1(1) - x2(1));
      return tie_breaker_ * sqrt(dx * dx + dy * dy);
    }

  public:
    KinoAstar(){};
    KinoAstar(const std::string &config_string, VoxelMap::Ptr &map_ptr, VehicleGeometry::Ptr &vehicle_geo_ptr);
    ~KinoAstar();
    typedef std::shared_ptr<KinoAstar> Ptr;

    VoxelMap::Ptr map_ptr_;
    VehicleGeometry::Ptr vehicle_geo_ptr_;
    std::vector<Eigen::Vector2d> car_vertex_;
    std::vector<std::vector<Eigen::Vector2d>> aux_vertices_list_;
    std::vector<Eigen::Vector2d> aux_vertices;

    enum
    {
      REACH_HORIZON = 1,
      REACH_END = 2,
      NO_PATH = 3,
      REACH_END_BUT_SHOT_FAILS = 4
    };

    /* main API */
    void reset();
    void findNearestNode(Eigen::Vector2d &start_pos, bool first_search);
    int search(Eigen::Vector4d start_state, Eigen::Vector2d init_ctrl,
               Eigen::Vector4d end_state, bool use3d = false);
    bool searchTime(plan_utils::KinoTrajData &flat_trajs);
    void getTruncatedposLists(double truncate_len = 30);
    void getFullposeLists();
    void getSingulNodes();
    Eigen::Vector3d CalculateInitPos(double &t, int &singul);

    void setMap(const VoxelMap::Ptr &ptr) { map_ptr_ = ptr; };
    void setVehicleGeometry(const VehicleGeometry::Ptr &ptr) { vehicle_geo_ptr_ = ptr; };

    // get kino traj for optimization
    void getKinoNode(plan_utils::KinoTrajData &flat_trajs);
    bool checkCollisionUsingPosAndYaw(const Eigen::Vector3d &state, bool &res);
    int searchUseCircle(Eigen::Vector4d start_state, Eigen::Vector2d init_ctrl,
                        Eigen::Vector4d end_state, bool use3d);
    /* circle point check collision */
    bool checkVirtualCollisionAtPoint(const Eigen::Vector2d &point) const;
    bool checkVirtualCollisionsAtPoints(const std::vector<Eigen::Vector2d> &points) const;
    bool circleCheckCollisionUsingPosAndYaw(const Eigen::Vector3d &state, bool &res);
    std::vector<Eigen::Vector2d> circlePointTransformation(const geometry_cover::CoverCircleResult &result, const Eigen::Vector3d &state);

    const geometry_cover::CoverCircleResult getCoverCircle() const;

    /*hzchzc*/
    Eigen::Vector3d evaluatePos(double t);
    std::vector<Eigen::Vector4d> SamplePosList(int N); // px py yaw t
    std::vector<Eigen::Vector3d> getSampleTraj();
    double totalTrajTime;
    double checkl = 0.2;

    std::vector<PathNodePtr> getVisitedNodes();

    std::vector<Eigen::Vector4d> state_list;
    std::vector<Eigen::Vector3d> acc_list;
    std::vector<Eigen::Vector3d> pose_list;
    std::vector<Eigen::Vector2d> control_input_list;
    geometry_cover::CoverCircleResult cover_circle_result_;
    std::vector<Eigen::Vector2d> flatInputs_; // To store the inputs at each pose point

  };

} // namespace path_searching

#endif // _KINODYNAMIC_ASTAR_HPP_
