#include <sstream>
#include "path_searching/kino_astar.hpp"

namespace path_searching
{

using namespace std;
using namespace map_utils;
  KinoAstar::~KinoAstar()
  {
    // for (int i = 0; i < allocate_num_; i++)
    // {
    //   delete path_node_pool_[i];
    // }
  }

  KinoAstar::KinoAstar(const std::string& config_string, VoxelMap::Ptr& map_ptr, VehicleGeometry::Ptr& vehicle_geo_ptr)
  {    
      // set parameters
      YAML::Node config = YAML::Load(config_string);
      auto search_config = config["search"];

      max_horizon_ = search_config["max_horizon"].as<double>(50.0);
      end_horizon_ = search_config["end_horizon"].as<double>(10.0);
      yaw_resolution_ = search_config["yaw_resolution"].as<double>(0.3);
      lambda_heu_ = search_config["lambda_heu"].as<double>(5.0);
      allocate_num_ = search_config["allocate_num"].as<int>(100000);
      check_num_ = search_config["check_num"].as<int>(5);
      max_seach_time = search_config["max_search_time"].as<double>(1000.1);
      traj_forward_penalty = search_config["traj_forward_penalty"].as<double>(1.0);
      traj_back_penalty = search_config["traj_back_penalty"].as<double>(2.5);
      traj_gear_switch_penalty = search_config["traj_gear_switch_penalty"].as<double>(15.0);
      traj_steer_penalty = search_config["traj_steer_penalty"].as<double>(0.5);
      traj_steer_change_penalty = search_config["traj_steer_change_penalty"].as<double>(0.0);
      step_arc_ = search_config["step_arc"].as<double>(1.0);
      step_arc_res_ = search_config["step_arc_res"].as<double>(0.5);
      step_steer_res_ = search_config["step_steer_res"].as<double>(0.5);
      checkl = search_config["checkl"].as<double>(0.2);
      max_vel_ = search_config["max_vel"].as<double>(10.0);
      max_acc_ = search_config["max_acc"].as<double>(5.0);
      max_cur_ = search_config["max_cur"].as<double>(0.3);
      max_steer_ = search_config["max_steer"].as<double>(50.0);
      time_resolution_ = search_config["time_resolution"].as<double>(0.1);
      distance_resolution_ = search_config["distance_resolution"].as<double>(0.5);
      velocity_resolution_ = search_config["velocity_resolution"].as<double>(0.5);
      non_siguav = search_config["non_siguav"].as<double>(0.05);
      collision_check_type_ = search_config["collision_check_type"].as<int>(0);

      // set map
      setMap(map_ptr);

      // set vehicle geometry
      car_wheelbase_ = vehicle_geo_ptr->getWheelbase();
      setVehicleGeometry(vehicle_geo_ptr);

      // set

      /* ---------- pre-allocated node ---------- */
      path_node_pool_.resize(allocate_num_);
      middle_node_pool_.resize(allocate_num_);
      for (int i = 0; i < allocate_num_; i++)
      {
        path_node_pool_[i] = std::make_shared<PathNode>();
        middle_node_pool_[i] = std::make_shared<MiddleNode>();
      }
      path_nodes_.clear();
      middle_nodes_.clear();
      last_path_pos_.clear();
      last_path_pos_temp_.clear();

      use_node_num_ = 0;
      use_time_node_num_ = 0;
      iter_num_ = 0;

      resolution_ = map_ptr_->getResolution();
      yaw_origin_ = -M_PI;
      
      max_steer_ = max_steer_ * M_PI / 180.0;

      min_vel_ = -max_vel_;
      min_acc_ = -max_acc_;

      inv_resolution_ = 1.0 / resolution_;
      inv_yaw_resolution_ = 1.0 / yaw_resolution_;

      shotptr =std::make_shared<ompl::base::ReedsSheppStateSpace>(1.0 / max_cur_);
      
      non_siguav = fabs(non_siguav);

      // stores the vertexs of the car
      car_vertex_.clear();
      for (const auto &vertex: vehicle_geo_ptr->getCarVertices())
      {
          car_vertex_.push_back(vertex);
      }
      
      if (car_vertex_.front() != car_vertex_.back()) {
          car_vertex_.push_back(car_vertex_.front());
      }
      // stores the AuxVertices of the car
      // aux_vertices.clear();
      for (int i = 0; i < vehicle_geo_ptr->getAuxCount(); i++) {
          std::vector<Eigen::Vector2d> aux_vertices_part = vehicle_geo_ptr->getAuxVertices(i);
          aux_vertices_list_.push_back(aux_vertices_part);
      }

      if (collision_check_type_ == 1)
      {
        cover_circle_result_=vehicle_geo_ptr->getCoverCircle();
        // if (!aux_vertices.empty())
        //   {
        //     ;
        //     // std::cout << "[Kino Astar]: Computing aux cover circles..." << std::endl;
        //     // addAuxVerticesCoverCircleInInitPose();
        //   }
      }
  }

  void KinoAstar::stateTransit(Eigen::Vector3d &state0,  Eigen::Vector3d &state1,
              Eigen::Vector2d &ctrl_input)
  {
      //helpful var
      double psi = ctrl_input[0]; 
      double s = ctrl_input[1]; 
      if(fabs(psi) > 1e-3){
        double k = car_wheelbase_ / tan(psi);
        state1[0] = state0[0] + k*(sin(state0[2]+s/k)-sin(state0[2]));
        state1[1] = state0[1] - k*(cos(state0[2]+s/k)-cos(state0[2]));
        state1[2] = state0[2] + s/k;
      }
      else{
        state1[0] = state0[0] + s * cos(state0[2]);
        state1[1] = state0[1] + s * sin(state0[2]);
        state1[2] = state0[2]; 
      }
  }

  bool KinoAstar::checkVirtualCollisionAtPoint(const Eigen::Vector2d &point) const 
  {
      // Check if point is within map boundaries
      if (!map_ptr_->isInMap(point))
      {
          return true; // Collision detected if point is outside map
      }

      // Check collision with obstacles in the map
      // !! Note: we consider virtual dilation here
      if (map_ptr_->getVoxelState(point) > map_utils::UNOCCUPIED)
      {
          return true; // Collision detected if point is inside an obstacle
      }
      else
      {
          return false; // No collision detected
      }
  }

  bool KinoAstar::checkVirtualCollisionsAtPoints(const std::vector<Eigen::Vector2d> &points) const 
  {
      for (const auto& point : points) {
          if (checkVirtualCollisionAtPoint(point))
            return true;
      }

      return false;
  }

  /*input state and cover cicles struct return states*/
  bool KinoAstar::circleCheckCollisionUsingPosAndYaw(const Eigen::Vector3d &state, bool &res)
  {  
      res = false;
      if (!map_ptr_->isInMap(Eigen::Vector2d(state[0], state[1]))) {
        res = true;
        return res;
      }

      std::vector<Eigen::Vector2d> points_to_check = circlePointTransformation(cover_circle_result_, state);
      
      // Check collisions at each point
      res = checkVirtualCollisionsAtPoints(points_to_check);
     
      return res;    
  }

  std::vector<Eigen::Vector2d> KinoAstar::circlePointTransformation(const geometry_cover::CoverCircleResult& result, const Eigen::Vector3d& state)
  {
    std::vector<Eigen::Vector2d> points_to_check;

    Eigen::Vector2d translation_vector(state[0], state[1]);
    Eigen::Matrix2d rotation_matrix;
    rotation_matrix << std::cos(state[2]), -std::sin(state[2]),
                       std::sin(state[2]),  std::cos(state[2]);

    // main body
    for (const auto& point : result.centers) {
        Eigen::Vector2d rotated_point = rotation_matrix * point;
        Eigen::Vector2d transformed_point = rotated_point + translation_vector;
        points_to_check.push_back(transformed_point);
    }
    // aux
    for (const auto& point : result.aux_centers) {
        Eigen::Vector2d rotated_point = rotation_matrix * point;
        Eigen::Vector2d transformed_point = rotated_point + translation_vector;
        points_to_check.push_back(transformed_point);
    }

    return points_to_check;
  }

  bool KinoAstar::checkCollisionUsingPosAndYaw(const Eigen::Vector3d &state, bool &res)
  {
      res = false;
      bool car_res=false;
      bool aux_res=false;
      if (!map_ptr_->isInMap(Eigen::Vector2d(state[0], state[1]))) {
        res = true;
        return res;
      }
      Eigen::Vector2d pos = state.head(2);
      double yaw = state[2];
      Eigen::Matrix2d Rotation_matrix;
      Rotation_matrix << cos(yaw),  -sin(yaw),
                         sin(yaw),  cos(yaw);
      
      for(int i = 0; i < car_vertex_.size()-1; i++)
      {
          Eigen::Vector2d start_point = pos + Rotation_matrix * car_vertex_[i];
          Eigen::Vector2d end_point = pos + Rotation_matrix * car_vertex_[i+1];
          if (!map_ptr_->isInMap(start_point) || !map_ptr_->isInMap(end_point))
          {
            car_res = true;
            break;
          }
          RayCaster car_raycaster;
          bool need_ray = car_raycaster.setInput(start_point / resolution_, end_point / resolution_);
          Eigen::Vector2d half(0.5, 0.5);
          Eigen::Vector2d ray_pt;
          while(car_raycaster.step(ray_pt))
          {
              Eigen::Vector2d tmp = (ray_pt + half) * resolution_;
              if (map_ptr_->checkPointCollision(tmp))
              {
                  car_res = true;
                  break;

              }
          }
      }
       if (!aux_vertices_list_.empty())
       {
          for(int i = 0; i < aux_vertices_list_.size(); i++)
          {
            // if(aux_res)
            // {
            //   break;
            // }
            for(int j = 0; j < aux_vertices_list_[i].size()-1;j++)
            {
                Eigen::Vector2d start_point = pos + Rotation_matrix * aux_vertices_list_[i][j];
                Eigen::Vector2d end_point = pos + Rotation_matrix * aux_vertices_list_[i][j+1];
                if (!map_ptr_->isInMap(start_point) || !map_ptr_->isInMap(end_point))
                {
                  aux_res = true;
                  break;
                }
                RayCaster aux_raycaster;
                bool need_ray = aux_raycaster.setInput(start_point / resolution_, end_point / resolution_);
                Eigen::Vector2d half(0.5, 0.5);
                Eigen::Vector2d ray_pt;
                while(aux_raycaster.step(ray_pt))
                {
                    Eigen::Vector2d tmp = (ray_pt + half) * resolution_;
                    if (map_ptr_->checkPointCollision(tmp))
                    {
                        aux_res = true;
                        break;
                    }
                }
            }
          }
       }
  
      if (aux_res==true || car_res==true)
      {
        res=true;
      }
      return res;
  }

  const geometry_cover::CoverCircleResult KinoAstar:: getCoverCircle() const
  {
    return cover_circle_result_;
  }

  int KinoAstar::searchUseCircle(Eigen::Vector4d start_state, Eigen::Vector2d init_ctrl,
                               Eigen::Vector4d end_state,bool use3d)
  {
    bool isocc = false;
    bool initsearch = false;
    auto t1 = chrono::high_resolution_clock::now();
    circleCheckCollisionUsingPosAndYaw(start_state.head(3), isocc);
    // checkCollisionUsingPosAndYaw(start_state.head(3), isocc);
    if(isocc){
      printf("[Kino Astar]: start is not free!\n");
      return NO_PATH;
    }
    // checkCollisionUsingPosAndYaw(end_state.head(3), isocc);
    circleCheckCollisionUsingPosAndYaw(end_state.head(3), isocc);
    if(isocc){
      printf("[Kino Astar]: end is not free!\n");
      return NO_PATH;
    }

    start_state_ = start_state;
    start_ctrl_ = init_ctrl;
    end_state_ = end_state;
    Eigen::Vector2i end_index;
    map_ptr_->posToIndex(end_state.head(2), end_index);
    /* ---------- initialize ---------- */
    
    PathNodePtr cur_node = path_node_pool_[0];
    cur_node->parent = NULL;
    cur_node->state = start_state.head(3);
    map_ptr_->posToIndex(start_state.head(2), cur_node->index);
    cur_node->yaw_idx = yawToIndex(start_state[2]);
    cur_node->g_score = 0.0;
    cur_node->input = Eigen::Vector2d(0.0,0.0);
    cur_node->singul = getSingularity(start_state(3));
    cur_node->number = 0;
    cur_node->f_score = lambda_heu_ * getHeu(cur_node->state, end_state);
    cur_node->node_state = IN_OPEN_SET;
    open_set_.push(cur_node);
    use_node_num_ += 1;
    if(!use3d)
      expanded_nodes_.insert(cur_node->index, cur_node);
    else
      expanded_nodes_.insert(cur_node->index, yawToIndex(start_state[2]),cur_node);
    PathNodePtr terminate_node = NULL;
    if(cur_node->singul == 0){ initsearch = true;}
    /* ---------- search loop ---------- */
    while (!open_set_.empty())
    {
      /* ---------- get lowest f_score node ---------- */
      cur_node = open_set_.top();

      /* ---------- determine termination ---------- */
      // to decide the near end
      // do not shot until near the target
      bool reach_end_horizon = (cur_node->state.head(2) - end_state_.head(2)).norm() <= end_horizon_;
      if(reach_end_horizon && initsearch) {
        is_shot_sucess_use_circle(cur_node->state, end_state_.head(3), cur_node->singul);
      }

      if (is_shot_succ_)
      {
        terminate_node = cur_node;
        retrievePath(terminate_node);
        has_path_ = true;
        /* one shot trajectory */
        // printf("[Kino Astar]: one shot! iter num: %d\n", iter_num_);
        return REACH_END;
      } 

      bool reach_horizon = (cur_node->state.head(2) - end_state_.head(2)).norm() > max_horizon_;
      if (reach_horizon)
      {
        cout << "[Kino Astar]: Reach horizon_" << endl;
        return REACH_HORIZON;
      }

      auto t2 = std::chrono::high_resolution_clock::now();
      double runTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() * 1000; // ms
      if(runTime > max_seach_time){
        terminate_node = cur_node;
        retrievePath(terminate_node);
        has_path_ = true;
        if (terminate_node->parent == NULL)
        {

          printf("[Kino Astar]: terminate_node->parent == NULL\n");
          printf("[Kino Astar]: NO_PATH\n");
          return NO_PATH;
        }
        else
        {
          printf("[Kino Astar]: Reach the max seach time\n");
          return REACH_END;
        }
      }
      /* ---------- pop node and add to close set ---------- */
      open_set_.pop();
      cur_node->node_state = IN_CLOSE_SET;
      iter_num_ += 1;
      /* ---------- init state propagation ---------- */
      Eigen::Vector3d cur_state = cur_node->state;
      Eigen::Vector3d pro_state;
      Eigen::Vector2d ctrl_input;
      std::vector<Eigen::Vector2d> inputs;

      if(!initsearch ){
        if(start_state_[3]>0){
          for(double arc = resolution_; arc <= 2*resolution_ + 1e-3; arc += resolution_){
            for(double steer = -max_steer_; steer <= max_steer_ + 1e-3; steer += step_steer_res_ * max_steer_ * 1.0){
              ctrl_input<< steer, arc;
              inputs.push_back(ctrl_input);
            }
          }
        }
        else{
          for(double arc = -resolution_; arc >= -2*resolution_ - 1e-3; arc-= resolution_){
            for(double steer = -max_steer_; steer <= max_steer_ + 1e-3; steer += step_steer_res_ * max_steer_ * 1.0){
              ctrl_input<< steer, arc;
              inputs.push_back(ctrl_input);
            }
          } 
        }
        initsearch = true;
      }
      else{
        for (double arc = -step_arc_; arc <= step_arc_ + 1e-3; arc += step_arc_res_*step_arc_){
          if(fabs(arc)<1.0e-2) continue;
          for (double steer = -max_steer_; steer <= max_steer_ + 1e-3; steer += step_steer_res_ * max_steer_ * 1.0)
          {
            ctrl_input << steer, arc;
            inputs.push_back(ctrl_input);
          }
        }
      }
      
      /* ---------- state propagation loop ---------- */
      for (auto& input:inputs){
        int singul = input[1]>0?1:-1;
        stateTransit(cur_state, pro_state, input);

        /* inside map range */
        Eigen::Vector2d pro_pos = pro_state.head(2);
        if( !map_ptr_->isInMap(pro_pos) )
        {
            // printf("[Kino Astar]: out of map range\n");
            continue;
        }
        /* not in close set */
        Eigen::Vector2i pro_id;
        map_ptr_->posToIndex(pro_state.head(2), pro_id);
        double pro_yaw_id = yawToIndex(pro_state[2]);
        PathNodePtr pro_node;
        if(use3d)
          pro_node = expanded_nodes_.find(pro_id, pro_yaw_id);
        else
          pro_node = expanded_nodes_.find(pro_id);
        if (pro_node != NULL && pro_node->node_state == IN_CLOSE_SET)
        {
          continue;
        }

        /* not in the same voxel */
        Eigen::Vector2i diff = pro_id - cur_node->index;
        int diff_yaw = pro_yaw_id - cur_node->yaw_idx;
        if (diff.norm() == 0 && ((!use3d) || diff_yaw == 0))
        {
          continue;
        }
        /* collision free */
        Eigen::Vector3d xt;
        bool is_occ = false;
        for (int k = 1; k <= check_num_; ++k)
        {
          double tmparc = input[1] * double(k) / double(check_num_);
          Eigen::Vector2d tmpctrl; tmpctrl << input[0],tmparc;
          stateTransit(cur_state, xt, tmpctrl);
          circleCheckCollisionUsingPosAndYaw(xt, is_occ);
          if (is_occ)
          {
            break;
          }

        }

        if (is_occ)  continue;

        /* ---------- compute cost ---------- */
        double tmp_g_score = 0.0;
        double tmp_f_score = 0.0;
        int lastDir = cur_node->singul;
        if(singul>0){
          tmp_g_score +=  std::fabs(input[1]) * traj_forward_penalty;
        }
        else{
          tmp_g_score += std::fabs(input[1]) * traj_back_penalty;
        }

        if(singul * lastDir < 0){
          tmp_g_score += traj_gear_switch_penalty;
        }
        tmp_g_score += traj_steer_penalty * std::fabs(input[0]) * std::fabs(input[1]);
        tmp_g_score += traj_steer_change_penalty * std::fabs(input[0]-cur_node->input[0]);
        tmp_g_score += cur_node->g_score;
        if(last_path_pos_.size() != 0)
        {
          if(cur_node->number < 20 && cur_node->number + nearest_idx_ < last_path_pos_.size())
          {
            Eigen::Vector2d last_nearest_pos = last_path_pos_[nearest_idx_ + cur_node->number];
            bool collision_between_two_nodes;
            collision_between_two_nodes = map_ptr_->checkLineCollision(pro_pos, last_nearest_pos);   
            if(collision_between_two_nodes && map_ptr_->getVoxelState(last_nearest_pos) == map_utils::UNOCCUPIED)   
            {
              // tmp_g_score += 100;
              continue;
            }      
          }
        }

        tmp_f_score = tmp_g_score + lambda_heu_ * getHeu(pro_state, end_state);
        // tmp_f_score = tmp_g_score + lambda_heu_ * getObsHeu(pro_state.head(2));

        /* ---------- compare expanded node in this loop ---------- */
        if (pro_node == NULL)
        {
          pro_node = path_node_pool_[use_node_num_];
          pro_node->index = pro_id;
          pro_node->state = pro_state;
          pro_node->yaw_idx = pro_yaw_id;
          pro_node->f_score = tmp_f_score;
          pro_node->g_score = tmp_g_score;
          pro_node->input = input;
          pro_node->parent = cur_node;
          pro_node->node_state = IN_OPEN_SET;
          pro_node->singul = singul;
          pro_node->number = cur_node->number + 1;
          open_set_.push(pro_node);
          if(use3d)
            expanded_nodes_.insert(pro_id, pro_yaw_id, pro_node);
          else
            expanded_nodes_.insert(pro_id,pro_node);
          use_node_num_ += 1;
          if (use_node_num_ == allocate_num_)
          {
            printf("[Kino Astar]: run out of memory.\n");
            return NO_PATH;
          }
        }
        else if (pro_node->node_state == IN_OPEN_SET)
        {

          if (tmp_g_score < pro_node->g_score)
          {
            pro_node->index = pro_id;
            pro_node->state = pro_state;
            pro_node->yaw_idx = pro_yaw_id;
            pro_node->f_score = tmp_f_score;
            pro_node->g_score = tmp_g_score;
            pro_node->input = input;
            pro_node->parent = cur_node;
            pro_node->singul = singul;
            pro_node->number = cur_node->number + 1;
          }
        }
        else
        {
          printf("[Kino Astar]: error type in searching: %d\n", pro_node->node_state);
        }
      }
    }

    /* ---------- open set empty, no path ---------- */
    printf("[Kino Astar]: open set empty, no path.\n");
    return NO_PATH;
  }


  int KinoAstar::search(Eigen::Vector4d start_state, Eigen::Vector2d init_ctrl,
                               Eigen::Vector4d end_state,bool use3d)
  {
    bool isocc = false;
    bool initsearch = false;
    auto t1 = chrono::high_resolution_clock::now();
    checkCollisionUsingPosAndYaw(start_state.head(3), isocc);
    if(isocc){
      printf("[Kino Astar]: start is not free!\n");
      return NO_PATH;
    }
    checkCollisionUsingPosAndYaw(end_state.head(3), isocc);
    if(isocc){
      printf("[Kino Astar]: end is not free!\n");
      return NO_PATH;
    }
    flatInputs_.clear();
    start_state_ = start_state;
    start_ctrl_ = init_ctrl;
    end_state_ = end_state;
    Eigen::Vector2i end_index;
    map_ptr_->posToIndex(end_state.head(2), end_index);
    /* ---------- initialize ---------- */
    
    PathNodePtr cur_node = path_node_pool_[0];
    cur_node->parent = NULL;
    cur_node->state = start_state.head(3);
    map_ptr_->posToIndex(start_state.head(2), cur_node->index);
    cur_node->yaw_idx = yawToIndex(start_state[2]);
    cur_node->g_score = 0.0;
    cur_node->input = Eigen::Vector2d(0.0,0.0);
    cur_node->singul = getSingularity(start_state(3));
    cur_node->number = 0;
    cur_node->f_score = lambda_heu_ * getHeu(cur_node->state, end_state);
    cur_node->node_state = IN_OPEN_SET;
    open_set_.push(cur_node);
    use_node_num_ += 1;
    if(!use3d)
      expanded_nodes_.insert(cur_node->index, cur_node);
    else
      expanded_nodes_.insert(cur_node->index, yawToIndex(start_state[2]),cur_node);
    PathNodePtr terminate_node = NULL;
    if(cur_node->singul == 0){ initsearch = true;}
    /* ---------- search loop ---------- */
    while (!open_set_.empty())
    {
      /* ---------- get lowest f_score node ---------- */
      cur_node = open_set_.top();

      /* ---------- determine termination ---------- */
      // to decide the near end
      // do not shot until near the target
      bool reach_end_horizon = (cur_node->state.head(2) - end_state_.head(2)).norm() <= end_horizon_;
      if(reach_end_horizon && initsearch) {
        is_shot_sucess(cur_node->state, end_state_.head(3), cur_node->singul);
      }

      if (is_shot_succ_)
      {
        terminate_node = cur_node;
        retrievePath(terminate_node);
        has_path_ = true;
        /* one shot trajectory */
        // printf("[Kino Astar]: one shot! iter num: %d\n", iter_num_);
        return REACH_END;
      } 

      bool reach_horizon = (cur_node->state.head(2) - end_state_.head(2)).norm() > max_horizon_;
      if (reach_horizon)
      {
        cout << "[Kino Astar]: Reach horizon_" << endl;
        return REACH_HORIZON;
      }

      auto t2 = std::chrono::high_resolution_clock::now();
      double runTime = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() * 1000; // ms
      if(runTime > max_seach_time){
        terminate_node = cur_node;
        retrievePath(terminate_node);
        has_path_ = true;
        if (terminate_node->parent == NULL)
        {

          printf("[Kino Astar]: terminate_node->parent == NULL\n");
          printf("[Kino Astar]: NO_PATH\n");
          return NO_PATH;
        }
        else
        {
          printf("[Kino Astar]: Reach the max seach time\n");
          return REACH_END;
        }
      }
      /* ---------- pop node and add to close set ---------- */
      open_set_.pop();
      cur_node->node_state = IN_CLOSE_SET;
      iter_num_ += 1;
      /* ---------- init state propagation ---------- */
      Eigen::Vector3d cur_state = cur_node->state;
      Eigen::Vector3d pro_state;
      Eigen::Vector2d ctrl_input;
      std::vector<Eigen::Vector2d> inputs;

      if(!initsearch ){
        if(start_state_[3]>0){
          for(double arc = resolution_; arc <= 2*resolution_ + 1e-3; arc += resolution_){
            for(double steer = -max_steer_; steer <= max_steer_ + 1e-3; steer += step_steer_res_ * max_steer_ * 1.0){
              ctrl_input<< steer, arc;
              inputs.push_back(ctrl_input);
            }
          }
        }
        else{
          for(double arc = -resolution_; arc >= -2*resolution_ - 1e-3; arc-= resolution_){
            for(double steer = -max_steer_; steer <= max_steer_ + 1e-3; steer += step_steer_res_ * max_steer_ * 1.0){
              ctrl_input<< steer, arc;
              inputs.push_back(ctrl_input);
            }
          } 
        }
        initsearch = true;
      }
      else{
        for (double arc = -step_arc_; arc <= step_arc_ + 1e-3; arc += step_arc_res_*step_arc_){
          if(fabs(arc)<1.0e-2) continue;
          for (double steer = -max_steer_; steer <= max_steer_ + 1e-3; steer += step_steer_res_ * max_steer_ * 1.0)
          {
            ctrl_input << steer, arc;
            inputs.push_back(ctrl_input);
          }
        }
      }
      
      /* ---------- state propagation loop ---------- */
      for (auto& input:inputs){
        int singul = input[1]>0?1:-1;
        stateTransit(cur_state, pro_state, input);

        /* inside map range */
        Eigen::Vector2d pro_pos = pro_state.head(2);
        if( !map_ptr_->isInMap(pro_pos) )
        {
            // printf("[Kino Astar]: out of map range\n");
            continue;
        }
        /* not in close set */
        Eigen::Vector2i pro_id;
        map_ptr_->posToIndex(pro_state.head(2), pro_id);
        double pro_yaw_id = yawToIndex(pro_state[2]);
        PathNodePtr pro_node;
        if(use3d)
          pro_node = expanded_nodes_.find(pro_id, pro_yaw_id);
        else
          pro_node = expanded_nodes_.find(pro_id);
        if (pro_node != NULL && pro_node->node_state == IN_CLOSE_SET)
        {
          continue;
        }

        /* not in the same voxel */
        Eigen::Vector2i diff = pro_id - cur_node->index;
        int diff_yaw = pro_yaw_id - cur_node->yaw_idx;
        if (diff.norm() == 0 && ((!use3d) || diff_yaw == 0))
        {
          continue;
        }
        /* collision free */
        Eigen::Vector3d xt;
        bool is_occ = false;
        for (int k = 1; k <= check_num_; ++k)
        {
          double tmparc = input[1] * double(k) / double(check_num_);
          Eigen::Vector2d tmpctrl; tmpctrl << input[0],tmparc;
          stateTransit(cur_state, xt, tmpctrl);
          checkCollisionUsingPosAndYaw(xt, is_occ);
          if (is_occ)
          {
            break;
          }

        }

        if (is_occ)  continue;

        /* ---------- compute cost ---------- */
        double tmp_g_score = 0.0;
        double tmp_f_score = 0.0;
        int lastDir = cur_node->singul;
        if(singul>0){
          tmp_g_score +=  std::fabs(input[1]) * traj_forward_penalty;
        }
        else{
          tmp_g_score += std::fabs(input[1]) * traj_back_penalty;
        }

        if(singul * lastDir < 0){
          tmp_g_score += traj_gear_switch_penalty;
        }
        tmp_g_score += traj_steer_penalty * std::fabs(input[0]) * std::fabs(input[1]);
        tmp_g_score += traj_steer_change_penalty * std::fabs(input[0]-cur_node->input[0]);
        tmp_g_score += cur_node->g_score;
        if(last_path_pos_.size() != 0)
        {
          if(cur_node->number < 20 && cur_node->number + nearest_idx_ < last_path_pos_.size())
          {
            Eigen::Vector2d last_nearest_pos = last_path_pos_[nearest_idx_ + cur_node->number];
            bool collision_between_two_nodes;
            collision_between_two_nodes = map_ptr_->checkLineCollision(pro_pos, last_nearest_pos);   
            if(collision_between_two_nodes && map_ptr_->getVoxelState(last_nearest_pos) == map_utils::UNOCCUPIED)   
            {
              // tmp_g_score += 100;
              continue;
            }      
          }
        }

        tmp_f_score = tmp_g_score + lambda_heu_ * getHeu(pro_state, end_state);
        // tmp_f_score = tmp_g_score + lambda_heu_ * getObsHeu(pro_state.head(2));

        /* ---------- compare expanded node in this loop ---------- */
        if (pro_node == NULL)
        {
          pro_node = path_node_pool_[use_node_num_];
          pro_node->index = pro_id;
          pro_node->state = pro_state;
          pro_node->yaw_idx = pro_yaw_id;
          pro_node->f_score = tmp_f_score;
          pro_node->g_score = tmp_g_score;
          pro_node->input = input;
          pro_node->parent = cur_node;
          pro_node->node_state = IN_OPEN_SET;
          pro_node->singul = singul;
          pro_node->number = cur_node->number + 1;
          open_set_.push(pro_node);
          if(use3d)
            expanded_nodes_.insert(pro_id, pro_yaw_id, pro_node);
          else
            expanded_nodes_.insert(pro_id,pro_node);
          use_node_num_ += 1;
          if (use_node_num_ == allocate_num_)
          {
            printf("[Kino Astar]: run out of memory.\n");
            return NO_PATH;
          }
        }
        else if (pro_node->node_state == IN_OPEN_SET)
        {

          if (tmp_g_score < pro_node->g_score)
          {
            pro_node->index = pro_id;
            pro_node->state = pro_state;
            pro_node->yaw_idx = pro_yaw_id;
            pro_node->f_score = tmp_f_score;
            pro_node->g_score = tmp_g_score;
            pro_node->input = input;
            pro_node->parent = cur_node;
            pro_node->singul = singul;
            pro_node->number = cur_node->number + 1;
          }
        }
        else
        {
          printf("[Kino Astar]: error type in searching: %d\n", pro_node->node_state);
        }
      }
    }

    /* ---------- open set empty, no path ---------- */
    printf("[Kino Astar]: open set empty, no path.\n");
    return NO_PATH;
  }

  bool checkDirectionChange(const std::vector<Eigen::Vector3d> &path_list, int &singul)
  {
    if (path_list.size() < 2)  return false;
      
    Eigen::Vector3d prev_pose, current_pose;
    Eigen::Vector2d motion_dir, vehicle_dir;

    for (int i = 0; i < path_list.size() - 1; i++) 
    {
      prev_pose = path_list[i];
      current_pose = path_list[i + 1];
      
      motion_dir = current_pose.head(2) - prev_pose.head(2);
      motion_dir.normalize();
      vehicle_dir << cos(prev_pose[2]), sin(prev_pose[2]);
      double dir = motion_dir.dot(vehicle_dir);

      if (dir * singul < 0) return true;
    }
    return false;
  }

  bool KinoAstar::is_shot_sucess(Eigen::Vector3d state1,Eigen::Vector3d state2, int &singul){
    
    std::vector<Eigen::Vector3d> path_list;
    double len,st;
    st = computeShotTraj(state1,state2,path_list,len);
    bool is_occ = false;
    for(unsigned int i = 0; i < path_list.size(); ++i){
        checkCollisionUsingPosAndYaw(path_list[i], is_occ);
        if(is_occ) return false;
    }

    // if (checkDirectionChange(path_list, singul)) return false;

    is_shot_succ_ = true;
    return true;
  }
  bool KinoAstar::is_shot_sucess_use_circle(Eigen::Vector3d state1,Eigen::Vector3d state2, int &singul){
  
  std::vector<Eigen::Vector3d> path_list;
  double len,st;
  st = computeShotTraj(state1,state2,path_list,len);
  bool is_occ = false;
  for(unsigned int i = 0; i < path_list.size(); ++i){
      circleCheckCollisionUsingPosAndYaw(path_list[i], is_occ);
      if(is_occ) return false;
  }
  
  // if (checkDirectionChange(path_list, singul)) return false;

  is_shot_succ_ = true;
  return true;
  }

  double KinoAstar::computeShotTraj(Eigen::Vector3d &state1, Eigen::Vector3d &state2,
                                    std::vector<Eigen::Vector3d> &path_list,
                                    double& len){
    namespace ob = ompl::base;
    namespace og = ompl::geometric;
    ob::ScopedState<> from(shotptr), to(shotptr), s(shotptr);
    from[0] = state1[0]; from[1] = state1[1]; from[2] = state1[2];
    to[0] = state2[0]; to[1] = state2[1]; to[2] = state2[2];
    std::vector<double> reals;
    len = shotptr->distance(from(), to());
    double sum_T = len/max_vel_;    

    for (double l = 0.0; l <=len; l += checkl)
    {
      shotptr->interpolate(from(), to(), l/len, s());
      reals = s.reals();
      path_list.push_back(Eigen::Vector3d(reals[0], reals[1], reals[2]));        
    }
  
    return sum_T;
  }



  // to retrieve the path to the correct order
  void KinoAstar::retrievePath(PathNodePtr end_node)
  {
    PathNodePtr cur_node = end_node;
    path_nodes_.push_back(cur_node);

    while (cur_node->parent != NULL)
    {
      flatInputs_.push_back(cur_node->input);
      cur_node = cur_node->parent;
      path_nodes_.push_back(cur_node);
    }

    reverse(path_nodes_.begin(), path_nodes_.end());
    reverse(flatInputs_.begin(), flatInputs_.end()); // Ensure inputs are in the correct order

    last_path_pos_temp_.clear();
    last_path_pos_temp_.resize(path_nodes_.size());
    for(int i = 0; i < last_path_pos_temp_.size(); i++)
    {
      last_path_pos_temp_[i] = path_nodes_[i]->state.head(2);
    }
  }

  void KinoAstar::reset()
  {
    expanded_nodes_.clear();
    path_nodes_.clear();

    std::priority_queue<PathNodePtr, std::vector<PathNodePtr>, NodeComparator> empty_queue;
    open_set_.swap(empty_queue);

    for (int i = 0; i < use_node_num_; i++)
    {
      PathNodePtr node = path_node_pool_[i];
      node->parent = NULL;
      node->node_state = NOT_EXPAND;
    }

    use_node_num_ = 0;
    iter_num_ = 0;
    is_shot_succ_ = false;
  }

  void KinoAstar::findNearestNode(Eigen::Vector2d& start_pos, bool first_search)
  {
    if(!first_search)
    {
      last_path_pos_.clear();
      last_path_pos_.resize(last_path_pos_temp_.size());
      for(int i = 0; i < last_path_pos_temp_.size(); i++)
      {
        last_path_pos_[i] = last_path_pos_temp_[i];
      }

      double min_distance = 10000.0;
      nearest_idx_ = 0;
      for(int i = 0; i < last_path_pos_.size(); i++)
      {
        double distance = (last_path_pos_[i] - start_pos).norm();
        if(distance < min_distance)
        {
          min_distance = distance;
          nearest_idx_ = i;
          nearest_init_distance_ = distance;
        }
      }      
    }
    else
    {
      last_path_pos_.clear();
    }
  }

  /*hzchzc px py yaw*/
  Eigen::Vector3d KinoAstar::evaluatePos(double t){
    t = std::min<double>(std::max<double>(0,t),totalTrajTime);
    double startvel = fabs(start_state_[3]);
    double endvel = fabs(end_state_[3]);
    int index = -1;
    double tmpT = 0;
    double CutTime = 0;
    //locate the local traj
    for(int i = 0;i<shot_timeList.size();i++){
      tmpT+=shot_timeList[i];
      if(tmpT>=t) {
        index = i; 
        CutTime = t-tmpT+shot_timeList[i];
        break;
        }
    }

    double initv = non_siguav,finv = non_siguav;
    if(index==0) {initv  = startvel;}
    if(index==shot_lengthList.size() - 1) finv = endvel;

    double localtime = shot_timeList[index];
    double locallength = shot_lengthList[index];
    int front = shotindex[index]; int back =  shotindex[index+1];
    std::vector<Eigen::Vector3d> localTraj;localTraj.assign(SampleTraj.begin()+front,SampleTraj.begin()+back+1);
    //find the nearest point
    double arclength = evaluateLength(CutTime,locallength,localtime,initv,finv);
    double tmparc = 0;
    for(int i = 0; i < localTraj.size()-1;i++){
      tmparc += (localTraj[i+1]-localTraj[i]).head(2).norm();
      if(tmparc>=arclength){
        double l1 = tmparc-arclength;
        double l = (localTraj[i+1]-localTraj[i]).head(2).norm();
        double l2 = l-l1;//l2
        Eigen::Vector3d state = l1/l*localTraj[i]+l2/l*localTraj[i+1];
        if(fabs(localTraj[i+1][2]-localTraj[i][2])>=M_PI){   
          double normalize_yaw = 0;
          if(localTraj[i+1][2]<=0){
            normalize_yaw = l1/l*localTraj[i][2]+l2/l*(localTraj[i+1][2]+2*M_PI);
          }
          else if(localTraj[i][2]<=0){
            normalize_yaw = l1/l*(localTraj[i][2]+2*M_PI)+l2/l*localTraj[i+1][2];
          }
          state[2] = normalize_yaw;
        }
        return state;
      }
    }
    return localTraj.back();
  }

  /*This function is used to get the full path list, which only includes the points along the path*/
  /*------"positionList_" is the main vector we are going to use in later code-------------------------*/
  void KinoAstar::getFullposeLists()
  {
      
      // start position and end position are included in this list
      positionList_.clear();
      pose_list.clear();
      control_input_list.clear();
      Eigen::Vector3d start_pos = start_state_.head(3);

      for(int i = 0; i < path_nodes_.size(); i++)
      {
          Eigen::Vector3d pose = path_nodes_[i]->state.head(3);
          // Eigen::Vector3d last_pos;
          // if(i == 0)
          // {
          //     last_pos = start_pos;
          // }
          // else
          // {
          //     last_pos = path_nodes_[i-1]->state.head(3);
          // }
          positionList_.push_back(pose.head(2));
          pose_list.push_back(pose);
          if (i > 0) {
            control_input_list.push_back(path_nodes_[i]->input);
          }
      }

      if(is_shot_succ_)
      {
          ompl::base::ScopedState<> from(shotptr), to(shotptr), s(shotptr);
          Eigen::Vector3d state1, state2;
          state1 = path_nodes_.back()->state.head(3);
          state2 = end_state_.head(3);
          from[0] = state1[0]; from[1] = state1[1]; from[2] = state1[2];
          to[0] = state2[0]; to[1] = state2[1]; to[2] = state2[2];
          double shotLength = shotptr->distance(from(), to());
          std::vector<double> reals;
          Eigen::Vector3d prev_pose = state1;

          for(double l = checkl; l < shotLength; l += checkl)
          {
              shotptr->interpolate(from(), to(), l/shotLength, s());
              reals = s.reals();
              positionList_.push_back(Eigen::Vector2d(reals[0], reals[1]));
              pose_list.push_back(Eigen::Vector3d(reals[0], reals[1], reals[2]));
              
              Eigen::Vector3d current_pose(reals[0], reals[1], reals[2]);
              Eigen::Vector2d ctrl_input;
              double s = checkl;
              Eigen::Vector2d motion_dir = current_pose.head(2) - prev_pose.head(2);
              motion_dir.normalize();
              Eigen::Vector2d vehicle_dir(cos(prev_pose[2]), sin(prev_pose[2]));
              if (motion_dir.dot(vehicle_dir) < 0) {
                  s = -s; 
              }
              
              if (fabs(current_pose[2] - prev_pose[2]) > 1e-3)
              {
                  double R = s/(current_pose[2] - prev_pose[2]);
                  double psi = atan(car_wheelbase_ / R);
                  ctrl_input[0] = psi;
                  ctrl_input[1] = s;
              }
              else
              {
                  ctrl_input[0] = 0;
                  ctrl_input[1] = s;
              }
              control_input_list.push_back(ctrl_input);
              prev_pose = current_pose;
          }
          if((positionList_.back() - end_state_.head(2)).norm() < checkl)
          {
              positionList_.pop_back();
              pose_list.pop_back();
          }
          positionList_.push_back(end_state_.head(2));
          pose_list.push_back(end_state_.head(3));
          control_input_list.push_back(Eigen::Vector2d(0, 0));
      }

  }
  
  /*This function is used to get the truncated path list, which only includes the points along the path*/
  /*------"positionList_" is the main vector we are going to use in later code-------------------------*/
  void KinoAstar::getTruncatedposLists(double truncate_len)
  {
      
      // double truncate_len = 30;
      bool exceed_len = false;
      // start position and end position are included in this list
      positionList_.clear();
      Eigen::Vector2d start_pos = start_state_.head(2);
      double tmp_len = 0.0;

      for(int i = 0; i < path_nodes_.size(); i++)
      {
          Eigen::Vector2d pos = path_nodes_[i]->state.head(2);
          Eigen::Vector2d last_pos;
          if(i == 0)
          {
              last_pos = start_pos;
          }
          else
          {
              last_pos = path_nodes_[i-1]->state.head(2);
          }
          tmp_len += (pos - last_pos).norm();
          if(tmp_len > truncate_len)
          {
              exceed_len = true;
              break;
          }
          else
          {
              positionList_.push_back(pos);
          }
      }
      if(exceed_len)
      {
          return;
      }

      if(is_shot_succ_)
      {
          ompl::base::ScopedState<> from(shotptr), to(shotptr), s(shotptr);
          Eigen::Vector3d state1, state2;
          state1 = path_nodes_.back()->state.head(3);
          state2 = end_state_.head(3);
          from[0] = state1[0]; from[1] = state1[1]; from[2] = state1[2];
          to[0] = state2[0]; to[1] = state2[1]; to[2] = state2[2];
          double shotLength = shotptr->distance(from(), to());
          std::vector<double> reals;
          for(double l = checkl; l < shotLength; l += checkl)
          {
              shotptr->interpolate(from(), to(), l/shotLength, s());
              reals = s.reals();
              positionList_.push_back(Eigen::Vector2d(reals[0], reals[1]));
              // tmp_len += checkl;
              // if(tmp_len > truncate_len)
              // {
              //     exceed_len = true;
              //     break;
              // }
          }
          if((positionList_.back() - end_state_.head(2)).norm() < checkl)
          {
              positionList_.pop_back();
          }
          positionList_.push_back(end_state_.head(2));
      }

  }

  // This function is used to find the sinularity point of the kinodynamic trajectory
  /*"direction_change_idx_" is used in later code, which stores the index of the singularity point*/
  // attention: direction_change_idx_ not only includes the index of the direction change point, 
  // but also includes the indexes of the first point and the final point
  void KinoAstar::getSingulNodes()
  {
      direction_change_idx_.clear();
      direction_change_idx_.push_back(0);
      if(positionList_.size() <= 2)
      {
          direction_change_idx_.push_back(positionList_.size() - 1);
          return;
      }
      for(int i = 1; i < positionList_.size() - 1; i++)
      {
          Eigen::Vector2d direction = positionList_[i] - positionList_[i-1];
          Eigen::Vector2d next_direction = positionList_[i+1] - positionList_[i];
          if(next_direction.dot(direction) < 0)
          {
              // The i-th point in path list is the sinuglarity point
              direction_change_idx_.push_back(i);
          }
      }
      direction_change_idx_.push_back(positionList_.size() - 1);

      //******************************************************************//
      /*new version comment all below*/
      std::vector<int> Vec_node_need_to_be_erased;
      for(int i = 1; i < direction_change_idx_.size(); i++)
      {
          if(direction_change_idx_[i] - direction_change_idx_[i-1] > 1)
          {
              int singul_idx = direction_change_idx_[i];
              int pre_singul_idx = direction_change_idx_[i] - 1;
              double distance = (positionList_[singul_idx] - positionList_[pre_singul_idx]).norm();
              if(distance < checkl)
              {
                  Vec_node_need_to_be_erased.push_back(pre_singul_idx);
              }
          }
      }
      for(int i = 0; i < Vec_node_need_to_be_erased.size(); i++)
      {
          positionList_.erase(positionList_.begin() + Vec_node_need_to_be_erased[i]);
          // Vec_node_need_to_be_erased[i]--;
          for(int j = i + 1; j < Vec_node_need_to_be_erased.size(); j++)
          {
              Vec_node_need_to_be_erased[j]--;
          }
      }


      direction_change_idx_.clear();
      direction_change_idx_.push_back(0);
      if(positionList_.size() <= 2)
      {
          direction_change_idx_.push_back(positionList_.size() - 1);
          return;
      }
      for(int i = 1; i < positionList_.size() - 1; i++)
      {
          Eigen::Vector2d direction = positionList_[i] - positionList_[i-1];
          Eigen::Vector2d next_direction = positionList_[i+1] - positionList_[i];
          if(next_direction.dot(direction) < 0)
          {
              // The i-th point in path list is the sinuglarity point
              direction_change_idx_.push_back(i);
          }
      }      

      direction_change_idx_.push_back(positionList_.size() - 1);
  }

  bool KinoAstar::searchTime(plan_utils::KinoTrajData &flat_trajs)
  {
      // TODO: generate flat_trajs for each segment
      /*-------Preliminaries: Judge the cars' moving direction------------*/
      int start_direction;
      double init_yaw = normalize_angle(start_state_(2));
      Eigen::Vector2d initAngleDir(cos(init_yaw), sin(init_yaw));
      Eigen::Vector2d initdir = positionList_[1] - positionList_[0];
      if(initAngleDir.dot(initdir) >= 0)
          start_direction = 1;
      else
          start_direction = -1;

      /*------------------------------------------------------------------*/
      // if(direction_change_idx_.size() > 2)
      //     return false;
      /*build s-t graph step1: discrete every point from the position list*/
      std::vector<double> distance2start;
      std::vector<Eigen::Vector2d> discrete_positionList;
      int end_index = direction_change_idx_[1];
       
      double distance = 0.0;
      distance2start.push_back(distance);
      for(int i = 1; i <= end_index; i++)
      {
          Eigen::Vector2d start_pos = positionList_[i-1];
          Eigen::Vector2d end_pos = positionList_[i];
          double node_distance = (end_pos - start_pos).norm();
          distance += node_distance;
          distance2start.push_back(distance);
      }
      total_s_ = distance;

      discrete_positionList.push_back(positionList_[0]);
      double car_c2r_dist = vehicle_geo_ptr_->getCenterToRearDist();
      for(double dis = distance_resolution_; dis < total_s_; dis += distance_resolution_)
      {
          int locate_idx = 0;
          for(int i = 1; i <= end_index; i++)
          {
              if(dis > distance2start[i-1] && dis <= distance2start[i])
              {
                  locate_idx = i;
              }
          }
          double d1 = distance2start[locate_idx-1];
          double deltad = dis - d1;
          Eigen::Vector2d pos1 = positionList_[locate_idx - 1];
          Eigen::Vector2d pos2 = positionList_[locate_idx];
          Eigen::Vector2d pos_direction = pos2 - pos1;
          double node_dis = (pos2 - pos1).norm();
          Eigen::Vector2d node_direction = pos2 - pos1;
          Eigen::Vector2d discrete_point = pos1 + deltad / node_dis * node_direction;
          double yaw_angle = atan2(start_direction * pos_direction(1), start_direction * pos_direction(0));
          Eigen::Matrix2d Rot_mat;
          Rot_mat << cos(yaw_angle), -sin(yaw_angle),
                     sin(yaw_angle), cos(yaw_angle);          
          discrete_point = discrete_point + Rot_mat * Eigen::Vector2d(car_c2r_dist, 0.0);
          discrete_positionList.push_back(discrete_point);
      }
      discrete_positionList.push_back(positionList_[end_index]);

      /*------------------------------------------------------------------*/
      
      /*build s-t graph step2: enumurate every point of other car's trajectory and check collision*/
      std::vector<Eigen::Vector3i> s_t_graph;
      int max_int_time = 0;
      
      // std::cout << "The number of s-t-collision parts are: " << s_t_graph.size() << std::endl;
      /*------------------------------------------------------------------------------------------*/
      /*Build s-t graph step3: allocate memeory for the s-t graph*/
      int max_int_s = discrete_positionList.size() - 2;
      int max_int_t = max_int_time;
      int s_voxels_num = max_int_s + 1;
      int t_voxels_num = max_int_t + 1;
      int max_voxels = s_voxels_num * t_voxels_num;
      std::vector<int> s_t_graph_vector;
      s_t_graph_vector.resize(max_voxels);
      fill(s_t_graph_vector.begin(), s_t_graph_vector.end(), 0);
      // it was not executed???
      for(int i = 0; i < s_t_graph.size(); i++)
      {
          int time_int = s_t_graph.at(i)(0);
          int distance_int = s_t_graph.at(i)(1);
          int time_distance_int = time_int + distance_int * t_voxels_num;
          s_t_graph_vector[time_distance_int] = 1;
      }
      /*----------------------------------------------------------*/

      /*-----------------------start to search!-------------------*/
      for(int i = 0; i < use_time_node_num_; i++)
      {
          MiddleNodePtr node = middle_node_pool_[i];
          node->parent = NULL;
      }
      expanded_middle_nodes_.clear();
      middle_nodes_.clear();
      use_time_node_num_ = 0;
      openSet_.clear();

      double plan_start_vel = abs(start_state_(3));
      MiddleNodePtr start_node;
      start_node = middle_node_pool_[0];
      start_node->start_s = 0.0;
      start_node->end_s = 0.0;
      start_node->start_vel = plan_start_vel;
      start_node->end_vel = plan_start_vel;
      start_node->start_t = 0.0;
      start_node->end_t = 0.0;
      start_node->s_t_v_idx = Eigen::Vector3i(t2idx(0.0), s2idx(0.0), v2idx(plan_start_vel));
      start_node->length = 0.0;
      start_node->distance2end = total_s_;
      start_node->acc = 0.0;
      start_node->g_score = 0.0;
      start_node->f_score = 1 * calculateMiniTime2singulPoint(plan_start_vel, total_s_);
      start_node->node_state = IN_OPEN_SET;
      start_node->parent = NULL;
      use_time_node_num_++;
      
      openSet_.insert(make_pair(start_node->f_score, start_node));
      expanded_middle_nodes_.insert(start_node->s_t_v_idx, start_node);

      MiddleNodePtr cur_node = NULL;
      double max_input = max_acc_;
      if(plan_start_vel * plan_start_vel > 2 * max_input * total_s_)
          max_input = plan_start_vel * plan_start_vel / (2 * total_s_) + 0.1;

      // for debugging //
      std::vector<MiddleNodePtr> debug_middle_nodes;
      debug_middle_nodes.push_back(start_node);
      //**************//

      while(!openSet_.empty())
      {
          cur_node = openSet_.begin()->second;
          openSet_.erase(openSet_.begin());
          cur_node->node_state = IN_CLOSE_SET;

          double pro_start_t = cur_node->end_t;
          double pro_start_s = cur_node->end_s;
          double pro_start_vel = cur_node->end_vel;       
          double pro_start2singuldistance = cur_node->distance2end;   

          /*---------one shot, which needs the end velcity and distance to singul point---------*/
          bool one_shot_success = true;
          double shot_duration = 0.0;
          std::vector<double> oneShot_trans_t;
          double d1 = (max_vel_ * max_vel_ - pro_start_vel * pro_start_vel) / (2 * max_input);
          double d2 = max_vel_ * max_vel_ / (2 * max_input);
          if(pro_start2singuldistance >= d1 + d2)
          {
              double t1 = (max_vel_ - pro_start_vel) / max_input;
              double t2 = (pro_start2singuldistance - d2 - d1) / max_vel_;
              double t3 = max_vel_ / max_input;
              shot_duration = t1 + t2 + t3;
              oneShot_trans_t.resize(3);
              oneShot_trans_t.at(0) = t1; oneShot_trans_t.at(1) = t1 + t2; oneShot_trans_t.at(2) = t1 + t2 + t3;
          }
          
          double d3 = pro_start_vel * pro_start_vel / (2 * max_input);
          if(pro_start2singuldistance >= d3 && pro_start2singuldistance < (d1 + d2))
          {
              double max_v = sqrt((2 * max_input * pro_start2singuldistance + pro_start_vel * pro_start_vel) / 2);
              double t1 = (max_v - pro_start_vel) / max_input;
              double t2 = max_v / max_input;
              shot_duration = t1 + t2;
              oneShot_trans_t.resize(2);
              oneShot_trans_t.at(0) = t1; oneShot_trans_t.at(1) = t1 + t2;
          }
          else if(pro_start2singuldistance < d3)
          {
              printf("[Kino Astar]: Too fast init speed, not possible to stop!\n");
              continue;
              // return false;
          }

          for(double shot_check_t = time_resolution_; shot_check_t <= shot_duration + 1e-3; shot_check_t += time_resolution_)
          {
              double cur_t, cur_s;
              cur_t = pro_start_t + shot_check_t;
              if(oneShot_trans_t.size() == 3)
              {
                  if(shot_check_t < oneShot_trans_t.at(0))
                  {
                      cur_s = pro_start_s + pro_start_vel * shot_check_t + 0.5 * max_input * shot_check_t * shot_check_t;
                  }
                  else if(shot_check_t >= oneShot_trans_t.at(0) && shot_check_t < oneShot_trans_t.at(1))
                  {
                      double s1 = (max_vel_ * max_vel_ - pro_start_vel * pro_start_vel) / (2 * max_input);
                      double s2 = max_vel_ * (shot_check_t - oneShot_trans_t.at(0));
                      cur_s = pro_start_s + s1 + s2;
                  }
                  else
                  {
                      double s1 = (max_vel_ * max_vel_ - pro_start_vel * pro_start_vel) / (2 * max_input);
                      double s2 = max_vel_ * (oneShot_trans_t.at(1) - oneShot_trans_t.at(0));
                      double s3 = max_vel_ * (shot_check_t - oneShot_trans_t.at(1)) - 0.5 * max_input * pow(shot_check_t - oneShot_trans_t.at(1), 2);
                      cur_s = pro_start_s + s1 + s2 + s3;
                  }

              }
              else if(oneShot_trans_t.size() == 2)
              {
                  if(shot_check_t < oneShot_trans_t.at(0))
                  {
                      cur_s = pro_start_s + pro_start_vel * shot_check_t + 0.5 * max_input * shot_check_t * shot_check_t;
                  }
                  else
                  {
                      double max_v = pro_start_vel + max_input * oneShot_trans_t.at(0);
                      double s1 = pro_start_vel * oneShot_trans_t.at(0) + 0.5 * max_input * pow(oneShot_trans_t.at(0), 2);
                      double s2 = max_v * (shot_check_t - oneShot_trans_t.at(0)) - 0.5 * max_input * pow(shot_check_t - oneShot_trans_t.at(0), 2);
                      cur_s = pro_start_s + s1 + s2;
                  }
              }
              int cur_t_int = t2idx(cur_t);
              if(cur_t_int >= max_int_t)
              {
                  one_shot_success = true;
                  break;
              }
              int cur_s_int = s2idx(cur_s);
              int cur_s_t_int = cur_t_int + cur_s_int * t_voxels_num;
              if(s_t_graph_vector.at(cur_s_t_int) > 0)
              {
                  one_shot_success = false;
                  break;
              }              
          }

          if(one_shot_success)
          {
              double shot_start_vel = cur_node->end_vel;
              double shot_start_s = cur_node->end_s;
              double shot_start_t = cur_node->end_t;              
              
              middle_nodes_.clear();
              middle_nodes_.push_back(cur_node);
              while(cur_node->parent != NULL)
              {
                  cur_node = cur_node->parent;
                  middle_nodes_.push_back(cur_node);
              }
              reverse(middle_nodes_.begin(), middle_nodes_.end()); // middle_nodes_ include the first no distance node
              
              flat_trajs.clear();
              
              double segment_duration = middle_nodes_.back()->end_t + shot_duration;
              double sampletime = 0.1;
              if(sampletime > segment_duration)
                  sampletime = segment_duration / 2.0;
              
              std::vector<Eigen::Vector3d> traj_pts;
              std::vector<double> thetas;
              double samplet;
              for(samplet = sampletime; samplet < segment_duration; samplet += sampletime)
              {
                  Eigen::Vector3d sample_pos = CalculateInitPos(samplet, start_direction);
                  // traj_pts.push_back(sample_pos.head(2));
                  Eigen::Vector2d traj_pt_pos = sample_pos.head(2);
                  Eigen::Vector3d traj_pt; traj_pt << traj_pt_pos, sampletime;
                  traj_pts.push_back(traj_pt);
                  thetas.push_back(sample_pos(2));
              }
              samplet -= sampletime; // to make sure segment_duration > samplet

              Eigen::Vector2d last_last_pos = positionList_.at(end_index - 1);
              Eigen::Vector2d last_pos = positionList_.at(end_index);
              Eigen::Vector2d last_direction = last_pos - last_last_pos;
              double segment_end_yaw = atan2(start_direction * last_direction(1), start_direction * last_direction(0));
              double segment_end_vel = start_direction > 0 ? non_siguav : -non_siguav;
              traj_pts.push_back(Eigen::Vector3d(last_pos(0), last_pos(1), segment_duration - samplet));
              thetas.push_back(normalize_angle(segment_end_yaw));

              Eigen::Vector4d segment_start_state, segment_end_state;

              segment_start_state = start_state_;
              segment_end_state << last_pos, segment_end_yaw, segment_end_vel;
              plan_utils::FlatTrajData flat_traj;
              Eigen::MatrixXd startS, endS;
              Eigen::Vector2d Initctrlinput(0, 0), Finctrlinput(0, 0);
              Initctrlinput = start_ctrl_;
              getFlatState(segment_start_state, Initctrlinput, startS, start_direction);
              getFlatState(segment_end_state, Finctrlinput, endS, start_direction);   
              flat_traj.traj_pts = traj_pts;
              flat_traj.thetas = thetas;
              flat_traj.start_state = startS;
              flat_traj.final_state = endS;
              flat_traj.singul = start_direction;
              flat_trajs.push_back(flat_traj);      
      
              /*-----------------------------------------------*/
              return true;
          }
          /*-----------------------------------------------------------------*/

          std::vector<double> ctrl_inputs, durations;
          double res = 1.0;
          for(double input = -max_input; input <= max_input + 1e-3; input += res * max_input)
          {
              ctrl_inputs.push_back(input);
          }
          durations.push_back(2 * time_resolution_);
          durations.push_back(4 * time_resolution_);
          durations.push_back(8 * time_resolution_);
          // for(double tau = 2 * time_resolution_; tau <= 3 * 2 * time_resolution_ + 1e-3; tau += 2 * time_resolution_)
          // {
          //     durations.push_back(tau);
          // }

          std::vector<MiddleNodePtr> expanded_nodes_from_one_node;
          expanded_nodes_from_one_node.clear();
          for(int i = 0; i < ctrl_inputs.size(); i++)
          {
              for(int j = 0; j < durations.size(); j++)
              {
                  double input = ctrl_inputs[i];
                  double duration = durations[j];

                  double pro_end_t = pro_start_t + duration;
                  double pro_end_s = pro_start_s + pro_start_vel * duration + 0.5 * input * duration * duration;
                  double pro_end_vel = pro_start_vel + input * duration;
                  double pro_end2singuldistance = total_s_ - pro_end_s;
                  Eigen::Vector3i pro_end_stv_idx(t2idx(pro_end_t), s2idx(pro_end_s), v2idx(pro_end_vel));
                  
                  if(pro_end_vel > max_vel_ || pro_end_vel < non_siguav)
                      continue;
                  if(pro_end_vel * pro_end_vel / (2 * max_input) > (pro_end2singuldistance - 1e-3))
                      continue;
                  if(pro_end_s > total_s_ || pro_end_s < pro_start_s)
                      continue;

                  MiddleNodePtr pro_node = expanded_middle_nodes_.find(pro_end_stv_idx);
                  if(pro_node != NULL && pro_node->node_state == IN_CLOSE_SET)
                      continue;

                  // check collision and dynamics
                  bool collision = false;
                  for(double check_t = pro_end_t - pro_start_t; check_t >= -1e-3; check_t -= time_resolution_)
                  {
                      double check_s = pro_start_s + pro_start_vel * check_t + 0.5 * input * check_t * check_t;
                      if(check_s > total_s_ || check_s < pro_start_s)
                      {
                          collision = true;
                          break;
                      }
                      // no need to check velocity, only check end velocity is enough. But it's necessary to check s.
                      // if t is beyond the max t, then no need to check any more, because there will not be collision
                      int check_int_t = t2idx(check_t + pro_start_t);
                      if(check_int_t >= max_int_t)
                          continue;
              
                      int check_int_s = s2idx(check_s);
                      int check_s_t_idx = check_int_t + t_voxels_num * check_int_s;

                      if(s_t_graph_vector.at(check_s_t_idx) > 0)
                      {
                          collision = true;
                          break;
                      }
                          
                  }
                  if(collision)
                      break;

                  double tmp_g_score = pro_end_t;
                  double tmp_h_score = 1 * calculateMiniTime2singulPoint(pro_end_vel, pro_end2singuldistance);
                  double tmp_f_score = tmp_g_score + tmp_h_score;
                  bool prune = false;
                  for(int k = 0; k < expanded_nodes_from_one_node.size(); k++)
                  {
                      MiddleNodePtr node_from_cur_node = expanded_nodes_from_one_node.at(k);
                      if(pro_end_stv_idx == node_from_cur_node->s_t_v_idx)
                      {
                          prune = true;
                          if(node_from_cur_node->f_score > tmp_f_score)
                          {
                              node_from_cur_node->end_s = pro_end_s;
                              node_from_cur_node->end_t = pro_end_t;
                              node_from_cur_node->end_vel = pro_end_vel;
                              node_from_cur_node->length = pro_start_vel * duration + 0.5 * input * duration * duration;
                              node_from_cur_node->distance2end = pro_end2singuldistance;
                              node_from_cur_node->acc = input;
                              node_from_cur_node->g_score = tmp_g_score;
                              node_from_cur_node->f_score = tmp_f_score;
                          }  
                          break;
                      }
                  }

                  if(!prune)
                  {
                      if(pro_node == NULL)
                      {
                          pro_node = middle_node_pool_[use_time_node_num_];
                          use_time_node_num_++;
                          pro_node->start_s = pro_start_s;
                          pro_node->start_t = pro_start_t;
                          pro_node->start_vel = pro_start_vel;
                          pro_node->end_s = pro_end_s;
                          pro_node->end_t = pro_end_t;
                          pro_node->end_vel = pro_end_vel;
                          pro_node->s_t_v_idx = pro_end_stv_idx;
                          pro_node->length = pro_start_vel * duration + 0.5 * input * duration * duration;
                          pro_node->distance2end = total_s_ - pro_end_s;
                          pro_node->acc = input;
                          pro_node->g_score = tmp_g_score;
                          pro_node->f_score = tmp_f_score;
                          pro_node->node_state = IN_OPEN_SET;
                          pro_node->parent = cur_node;

                          if(use_time_node_num_ > 20000)
                          {
                              printf("Too many nodes! Middle node search fails!\n");
                              return false;
                          }

                          debug_middle_nodes.push_back(pro_node);

                          openSet_.insert(make_pair(pro_node->f_score, pro_node));
                          expanded_middle_nodes_.insert(pro_node->s_t_v_idx, pro_node);
                          expanded_nodes_from_one_node.push_back(pro_node);
                      }
                      else if(pro_node->node_state == IN_OPEN_SET && pro_node->f_score > tmp_f_score)
                      {
                          pro_node->start_s = pro_start_s;
                          pro_node->start_t = pro_start_t;
                          pro_node->start_vel = pro_start_vel;
                          pro_node->end_s = pro_end_s;
                          pro_node->end_t = pro_end_t;
                          pro_node->end_vel = pro_end_vel;
                          pro_node->length = pro_start_vel * duration + 0.5 * input * duration * duration;
                          pro_node->distance2end = total_s_ - pro_end_s;
                          pro_node->acc = input;
                          pro_node->g_score = tmp_g_score;
                          pro_node->f_score = tmp_f_score;
                          pro_node->parent = cur_node;
                      }
                  }
              }
          }

      }
      printf("[Kino Astar]: MiddleNode open set empty! Search fails!\n");
      printf("[Kino Astar]: Number of expanded nodes: %zu\n", debug_middle_nodes.size());
      /*----------------------------------------------------------*/
      return false;
  }

  double KinoAstar::calculateMiniTime2singulPoint(double initspeed, double distance2end)
  {
      double d1 = (max_vel_ * max_vel_ - initspeed * initspeed) / (2 * max_acc_) + max_vel_ * max_vel_ / (2 * max_acc_);
      double d2 = distance2end - d1;
      if(d2 > 0)
      {
          double t1 = (max_vel_ - initspeed) / max_acc_;
          double t2 = d2 / max_vel_;
          double t3 = max_vel_ / max_acc_;
          return t1 + t2 + t3;
      }

      double d3 = initspeed * initspeed / (2 * max_acc_);
      if(distance2end >= d3)
      {
          double vmax = sqrt((2 * max_acc_ * distance2end + initspeed * initspeed) / 2.0);
          double t1 = (vmax - initspeed) / max_acc_;
          double t2 = vmax / max_acc_;
          return t1 + t2;
      }
      else
      {
          return 2 * distance2end / initspeed;
      }
  }

  int KinoAstar::s2idx(double s)
  {
      return floor((s + 1e-3) / distance_resolution_);
  }

  int KinoAstar::t2idx(double t)
  {
      return floor((t + 1e-3) / time_resolution_);
  }

  int KinoAstar::v2idx(double v)
  {
      return floor(v / velocity_resolution_);
  }

  Eigen::Vector3d KinoAstar::CalculateInitPos(double& t, int& singul)
  {
      double t_bar = t;
      double shot_start_t = middle_nodes_.back()->end_t;
      double shot_start_s = middle_nodes_.back()->end_s;
      double shot_start_v = middle_nodes_.back()->end_vel;
      double shot_length = total_s_ - shot_start_s;
      double cur_s;
      if(t_bar < shot_start_t - 1e-3)
      {
          for(MiddleNodePtr middle_node : middle_nodes_)
          {
              if(t_bar >= (middle_node->start_t - 1e-5) && t_bar < middle_node->end_t)
              {
                  double node_t_bar = t_bar - middle_node->start_t;
                  double node_start_v = middle_node->start_vel;
                  double node_acc = middle_node->acc;
                  cur_s = middle_node->start_s + node_start_v * node_t_bar + 0.5 * node_acc * pow(node_t_bar, 2);
                  break;
              }
          }
      }
      else
      {
          double shot_t_bar = t_bar - shot_start_t;
          double d1 = (max_vel_ * max_vel_ - shot_start_v * shot_start_v) / (2 * max_acc_);
          double d2 = max_vel_ * max_vel_ / (2 * max_acc_);
          double d3 = shot_start_v * shot_start_v / (2 * max_acc_);
          if(shot_length > (d1 + d2))
          {
              double t1 = (max_vel_ - shot_start_v) / max_acc_;
              double t2 = (shot_length - d2 - d1) / max_vel_;
              double t3 = max_vel_ / max_acc_;
              if(shot_t_bar < t1)
              {
                  cur_s = shot_start_s + shot_start_v * shot_t_bar + 0.5 * max_acc_ * pow(shot_t_bar, 2);
              }
              else if(shot_t_bar >= t1 && shot_t_bar < (t1 + t2))
              {
                  cur_s = shot_start_s + d1 + max_vel_ * (shot_t_bar - t1);
              }
              else
              {
                  cur_s = shot_start_s + d1 + max_vel_ * t2 + max_vel_ * (shot_t_bar - t1 - t2) - 0.5 * max_acc_ * pow(shot_t_bar - t1 - t2, 2);
              }
          }
          else if(shot_length <= d1 + d2 && shot_length > (d3 - 1e-3))
          {
              double max_v = sqrt((2 * max_acc_ * shot_length + shot_start_v * shot_start_v) / 2);
              double t1 = (max_v - shot_start_v) / max_acc_;
              double t2 = max_v / max_acc_;
              if(shot_t_bar < t1)
              {
                  cur_s = shot_start_s + shot_start_v * shot_t_bar + 0.5 * max_acc_ * shot_t_bar * shot_t_bar;
              }
              else
              {
                  cur_s = shot_start_s + (max_v * max_v - shot_start_v * shot_start_v) / (2 * max_acc_) + max_v * (shot_t_bar - t1) - 0.5 * max_acc_ * pow(shot_t_bar - t1, 2);
              }
          }
          else
          {
              // printf("[Kino Astar]: Wrong! Cannot calculate the init pos!\n");
              double max_input = shot_start_v * shot_start_v / (2 * shot_length);
              cur_s = shot_start_s + shot_start_v * shot_t_bar - 0.5 * max_input * pow(shot_t_bar, 2);
          }
      }

      double path_distance = 0.0;
      Eigen::Vector2d path_node_start_pos, path_node_end_pos;
      for(int i = 0; i < positionList_.size() - 1; i++)
      {
          double cur_path_distance = (positionList_.at(i+1) - positionList_.at(i)).norm();
          if(cur_s >= path_distance - 1e-3 && cur_s < path_distance + cur_path_distance)
          {
              path_node_start_pos = positionList_.at(i);
              path_node_end_pos = positionList_.at(i+1);
              break;
          }
          path_distance += cur_path_distance;
      }
      double s_bar = cur_s - path_distance;
      Eigen::Vector2d path_node_direction = path_node_end_pos - path_node_start_pos;
      double path_node_dis = path_node_direction.norm();
      Eigen::Vector2d sample_pos = path_node_start_pos + s_bar / path_node_dis * path_node_direction;
      double sample_yaw = atan2(singul * path_node_direction(1), singul * path_node_direction(0));
      
      return Eigen::Vector3d(sample_pos(0), sample_pos(1), normalize_angle(sample_yaw));
  } 

  // This function returns the distance between current node's "start" to the nearest singul point
  double KinoAstar::calculateCurNode2nearestSingul(int idx)
  {
      for(int i = 0; i < direction_change_idx_.size() - 1; i++)
      {
          if(idx > direction_change_idx_[i] && idx < direction_change_idx_[i+1])
          {
              double length_to_nearest_singul_point = 0.0;
              for(int j = idx - 1; j < direction_change_idx_[i+1]; j++)
              {
                  length_to_nearest_singul_point += (positionList_[j+1] - positionList_[j]).norm();
              }
              return length_to_nearest_singul_point;
          }
      }
  }

  // The time heuristic is the minimum time from this point to the final point
  double KinoAstar::gettimeHeu(int& idx, double& vel)
  {
      std::vector<double> remaining_length_vector;
      for(int i = 0; i < direction_change_idx_.size() - 1; i++)
      {
          if(idx > direction_change_idx_[i] && idx <= direction_change_idx_[i+1])
          {
              double length_to_nearest_singul_point = 0.0;
              for(int j = idx; j < direction_change_idx_[i+1]; j++)
              {
                  length_to_nearest_singul_point += (positionList_[j+1] - positionList_[j]).norm();
              }
              remaining_length_vector.push_back(length_to_nearest_singul_point);
          }
          else if(direction_change_idx_[i] > idx)
          {
              int start_idx = direction_change_idx_[i];
              int end_idx = direction_change_idx_[i+1];
              double segment_length = 0.0;
              for(int j = start_idx; j < end_idx; j++)
              {
                  segment_length += (positionList_[j+1] - positionList_[j]).norm();
              }
              remaining_length_vector.push_back(segment_length);
          }
      }
      double min_time_to_end = 0.0;
      for(int i = 0; i < remaining_length_vector.size(); i++)
      {
          double length = remaining_length_vector[i];
          if(i == 0)
          {
              double velocity = abs(vel);
              double l1 = (max_vel_ * max_vel_ - velocity * velocity) / (2 * max_acc_);
              double l2 = (max_vel_ * max_vel_) / (2 * max_acc_);
              double l3 = velocity * velocity / (2 * max_acc_);
              if(length > l1 + l2)
              {
                  double t1 = (max_vel_ - velocity) / max_acc_;
                  double t2 = (length - l1 - l2) / max_vel_;
                  double t3 = max_vel_ / max_acc_;
                  min_time_to_end += (t1 + t2 + t3);
              }
              else if(length > l3)
              {
                  double max_speed = sqrt(max_acc_ * length + 0.5 * velocity * velocity);
                  double t1 = (max_speed - velocity) / max_acc_;
                  double t2 = max_speed / max_acc_;
                  min_time_to_end += (t1 + t2);
              }
              else
              {
                  if(length < 0.05)
                  {
                      min_time_to_end += 0.0;
                  }
                  else
                  {
                      double a_max = velocity * velocity / (2 * length);
                      double t_dec = velocity / a_max;
                      min_time_to_end += t_dec;
                  }
              }
          }
          else
          {
              double l = max_vel_ * max_vel_ / max_acc_;
              if(length > l)
              {
                  double t1 = max_vel_ / max_acc_;
                  double t2 = (length - l) / max_vel_;
                  double t3 = t1;
                  min_time_to_end += (t1 + t2 + t3);
              }
              else
              {
                  double v_max = sqrt(length * max_acc_);
                  min_time_to_end += 2 * (v_max / max_acc_);
              }
          }
      }
      return 5 * min_time_to_end;
  }
  
  void KinoAstar::getKinoNode(plan_utils::KinoTrajData &flat_trajs)
  {
    // double truncate_len = 30;
    // bool exceed_len = false;

    flat_trajs.clear();
    std::vector<Eigen::Vector3d> roughSampleList;
    double startvel = fabs(start_state_[3]);
    double endvel = fabs(end_state_[3]);
    PathNodePtr node = path_nodes_.back();
    std::vector<Eigen::Vector3d> traj_pts;  // 3, N
    std::vector<double> thetas;
    Eigen::Vector4d x0, xt;
    Eigen::Vector2d ut, u0;
    while(node->parent != NULL){
      for (int k = check_num_; k >0; k--)
      {
        Eigen::Vector3d state;
        double tmparc = node->input[1] * double(k) / double(check_num_);
        Eigen::Vector2d tmpctrl; tmpctrl << node->input[0],tmparc;
        stateTransit(node->parent->state, state, tmpctrl);
        state[2] = normalize_angle(state[2]);
        roughSampleList.push_back(state);

      }
      node = node->parent;
    } 
    start_state_[2] = normalize_angle(start_state_[2]);
    roughSampleList.push_back(start_state_.head(3));
    reverse(roughSampleList.begin(),roughSampleList.end());

    if(is_shot_succ_){
      ompl::base::ScopedState<> from(shotptr), to(shotptr), s(shotptr);
      Eigen::Vector3d state1,state2;
      state1 = roughSampleList.back();
      state2 = end_state_.head(3);
      from[0] = state1[0]; from[1] = state1[1]; from[2] = state1[2];
      to[0] = state2[0]; to[1] = state2[1]; to[2] = state2[2];
      double shotLength = shotptr->distance(from(), to());
      std::vector<double> reals;
      for(double l = checkl; l < shotLength; l += checkl){
        shotptr->interpolate(from(), to(), l/shotLength, s());
        reals = s.reals();
        roughSampleList.push_back(Eigen::Vector3d(reals[0], reals[1], normalize_angle(reals[2])));        
      }
      end_state_[2] = normalize_angle(end_state_[2]);
      roughSampleList.push_back(end_state_.head(3));
    }
    // //truncate the init trajectory
    // double tmp_len = 0;
    // int truncate_idx = 0;
    // for(truncate_idx = 0;truncate_idx <roughSampleList.size()-1;truncate_idx++){  
    //   tmp_len += (roughSampleList[truncate_idx+1]-roughSampleList[truncate_idx]).norm();
    //   if(tmp_len>truncate_len){
    //     break;
    //   }
    // }
    // roughSampleList.assign(roughSampleList.begin(),roughSampleList.begin()+truncate_idx+1);
    SampleTraj = roughSampleList;
    /*divide the whole shot traj into different segments*/   
    shot_lengthList.clear();
    shot_timeList.clear();
    shotindex.clear();
    shot_SList.clear(); 
    double tmpl = 0;
    bool ifnewtraj = false;
    int lastS = (SampleTraj[1]-SampleTraj[0]).head(2).dot(Eigen::Vector2d(cos(SampleTraj[0][2]),sin(SampleTraj[0][2])))>=0?1:-1;
    //hzchzc attention please!!! max_v must = min_v  max_acc must = min_acc
    shotindex.push_back(0);
    for(int i = 0; i<SampleTraj.size()-1; i++){
      Eigen::Vector3d state1 = SampleTraj[i];
      Eigen::Vector3d state2 = SampleTraj[i+1];
      int curS = (state2-state1).head(2).dot(Eigen::Vector2d(cos(state1[2]),sin(state1[2]))) >=0 ? 1:-1;
      if(curS*lastS >= 0){
        tmpl += (state2-state1).head(2).norm();
      }
      else{  
        shotindex.push_back(i);
        shot_SList.push_back(lastS);
        shot_lengthList.push_back(tmpl);
        shot_timeList.push_back(evaluateDuration(tmpl,non_siguav,non_siguav));

        tmpl = (state2-state1).head(2).norm();
      }       
      lastS = curS;
    }
    shot_SList.push_back(lastS);
    shot_lengthList.push_back(tmpl);
    shot_timeList.push_back(evaluateDuration(tmpl,non_siguav,non_siguav));
    shotindex.push_back(SampleTraj.size()-1);
    if(shot_timeList.size()>=2){
      shot_timeList[0] = evaluateDuration(shot_lengthList[0],startvel,non_siguav);
      shot_timeList[shot_timeList.size()-1] = evaluateDuration(shot_lengthList.back(),non_siguav,endvel);
    }
    else{
      shot_timeList[0] = evaluateDuration(shot_lengthList[0],startvel,endvel);
    }
    /*extract flat traj  
    the flat traj include the end point but not the first point*/
    for(int i=0;i<shot_lengthList.size();i++){
      double initv = non_siguav,finv = non_siguav;
      Eigen::Vector2d Initctrlinput,Finctrlinput;
      Initctrlinput<<0,0;Finctrlinput<<0,0;
      if(i==0) {initv  = startvel; Initctrlinput = start_ctrl_;}
      if(i==shot_lengthList.size() - 1 && fabs(endvel) > non_siguav) finv = endvel;

      double locallength = shot_lengthList[i];
      int sig = shot_SList[i];
      std::vector<Eigen::Vector3d> localTraj;localTraj.assign(SampleTraj.begin()+shotindex[i],SampleTraj.begin()+shotindex[i+1]+1);
      traj_pts.clear();
      thetas.clear();        
      double samplet;
      double tmparc = 0;
      int index = 0;
      double sampletime = 0.1;
      if(shot_timeList[i]<=sampletime){
        sampletime = shot_timeList[i] / 2.0;
      }
      for(samplet = sampletime; samplet<shot_timeList[i]; samplet+=sampletime){
        double arc = evaluateLength(samplet,locallength,shot_timeList[i], initv, finv);
        //find the nearest point
        for(int k = index; k<localTraj.size()-1 ;k++)
        {
          tmparc+= (localTraj[k+1]-localTraj[k]).head(2).norm();
          if(tmparc>=arc){
            index = k;
            double l1 = tmparc-arc;
            double l = (localTraj[k+1]-localTraj[k]).head(2).norm();
            double l2 = l-l1;//l2
            double px = (l1/l*localTraj[k]+l2/l*localTraj[k+1])[0];
            double py = (l1/l*localTraj[k]+l2/l*localTraj[k+1])[1];
            double yaw= (l1/l*localTraj[k]+l2/l*localTraj[k+1])[2];
            if(fabs(localTraj[k+1][2]-localTraj[k][2])>=M_PI){   
              double normalize_yaw = 0;
              if(localTraj[k+1][2]<=0){
                normalize_yaw = l1/l*localTraj[k][2]+l2/l*(localTraj[k+1][2]+2*M_PI);
              }
              else if(localTraj[k][2]<=0){
                normalize_yaw = l1/l*(localTraj[k][2]+2*M_PI)+l2/l*localTraj[k+1][2];
              }
              yaw = normalize_yaw;
            }
            traj_pts.push_back(Eigen::Vector3d(px,py,sampletime));
            thetas.push_back(yaw);
            tmparc -=(localTraj[k+1]-localTraj[k]).head(2).norm();
            break;
          }
        }      
      }
      traj_pts.push_back(Eigen::Vector3d(localTraj.back()[0],localTraj.back()[1],shot_timeList[i]-(samplet-sampletime)));
      thetas.push_back(localTraj.back()[2]);
      plan_utils::FlatTrajData flat_traj;
      Eigen::MatrixXd startS;
      Eigen::MatrixXd endS;
      getFlatState(Eigen::Vector4d(localTraj.front()[0],localTraj.front()[1],localTraj.front()[2],sig*initv),Initctrlinput,startS,sig);
      getFlatState(Eigen::Vector4d(localTraj.back()[0],localTraj.back()[1],localTraj.back()[2],sig*finv),Finctrlinput,endS,sig);     
      flat_traj.traj_pts = traj_pts;
      flat_traj.thetas = thetas;
      flat_traj.start_state = startS;
      flat_traj.final_state = endS;
      flat_traj.singul = sig;
      flat_trajs.push_back(flat_traj);
    }
    totalTrajTime = 0.0;
    for(const auto dt : shot_timeList){
       totalTrajTime += dt; 
    }
  }

  double KinoAstar::evaluateDuration(double length, double startV, double endV){
   double critical_len; //the critical length of two-order optimal control, acc is the input
    if(startV>max_vel_||endV>max_vel_){
      printf("[Kino Astar]: evaluateDuration:start or end vel is larger that the limit!\n");
    }
    double startv2 = pow(startV,2);
    double endv2 = pow(endV,2);
    double maxv2 = pow(max_vel_,2);
    critical_len = (maxv2-startv2)/(2*max_acc_)+(maxv2-endv2)/(2*max_acc_);
    if(length>=critical_len){
      return (max_vel_-startV)/max_acc_+(max_vel_-endV)/max_acc_+(length-critical_len)/max_vel_;
    }
    else{
      double tmpv = sqrt(0.5*(startv2+endv2+2*max_acc_*length));
      return (tmpv-startV)/max_acc_ + (tmpv-endV)/max_acc_;
    }

  }

  double KinoAstar::evaluateLength(double curt,double locallength,double localtime,double startV, double endV){
    double critical_len; //the critical length of two-order optimal control, acc is the input
    if(startV>max_vel_||endV>max_vel_){
      printf("[Kino Astar]: evaluateLength:start or end vel is larger that the limit!\n");
    }
    double startv2 = pow(startV,2);
    double endv2 = pow(endV,2);
    double maxv2 = pow(max_vel_,2);
    critical_len = (maxv2-startv2)/(2*max_acc_)+(maxv2-endv2)/(2*max_acc_);
    if(locallength>=critical_len){
      double t1 = (max_vel_-startV)/max_acc_;
      double t2 = t1+(locallength-critical_len)/max_vel_;
      if(curt<=t1){
        return startV*curt + 0.5*max_acc_*pow(curt,2);
      }
      else if(curt<=t2){
        return startV*t1 + 0.5*max_acc_*pow(t1,2)+(curt-t1)*max_vel_;
      }
      else{
        return startV*t1 + 0.5*max_acc_*pow(t1,2) + (t2-t1)*max_vel_ + max_vel_*(curt-t2)-0.5*max_acc_*pow(curt-t2,2);
      }
    }
    else{
      double tmpv = sqrt(0.5*(startv2+endv2+2*max_acc_*locallength));
      double tmpt = (tmpv-startV)/max_acc_;
      if(curt<=tmpt){
        return startV*curt+0.5*max_acc_*pow(curt,2);
      }
      else{
        return startV*tmpt+0.5*max_acc_*pow(tmpt,2) + tmpv*(curt-tmpt)-0.5*max_acc_*pow(curt-tmpt,2);
      }
    }
  }

  const double kPi = acos(-1.0);
  double KinoAstar::normalize_angle(const double& theta)
  {
    double theta_tmp = theta;
    theta_tmp -= (theta >= kPi) * 2 * kPi;
    theta_tmp += (theta < -kPi) * 2 * kPi;
    return theta_tmp;
  }

  int KinoAstar::yawToIndex(double& yaw)
  {
    yaw = normalize_angle(yaw);
    int idx = floor((yaw - yaw_origin_) * inv_yaw_resolution_);
    return idx;
  }

  void KinoAstar::getFlatState(Eigen::Vector4d state, Eigen::Vector2d control_input,
                                  Eigen::MatrixXd &flat_state, int singul)
  {

    flat_state.resize(2, 3);

    double angle = state(2);
    double vel   = state(3); // vel > 0 

    Eigen::Matrix2d init_R;
    init_R << cos(angle),  -sin(angle),
              sin(angle),   cos(angle);

    // if (abs(vel) <= non_siguav){
    //   vel = singul * non_siguav;
    // }
    // else{
    //   vel = singul * vel;
    // }
    
    flat_state << state.head(2), init_R*Eigen::Vector2d(vel, 0.0), 
                  init_R*Eigen::Vector2d(control_input(1), std::tan(control_input(0)) / car_wheelbase_ * std::pow(vel, 2));

  }

} // namespace path_searching