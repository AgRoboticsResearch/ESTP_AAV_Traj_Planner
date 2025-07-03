#ifndef _TRAJ_CONTAINER_H_
#define _TRAJ_CONTAINER_H_

#include <Eigen/Eigen>
#include <vector>

#include "plan_utils/poly_traj_utils.hpp"

namespace plan_utils
{

    struct FlatTrajData
    {
        int singul;                            // direction
        std::vector<Eigen::Vector3d> traj_pts; // x, y, duration (3, N)
        std::vector<double> thetas;            // angle (N)
        Eigen::MatrixXd start_state;           // start flat state (2, 3)
        Eigen::MatrixXd final_state;           // end flat state (2, 3)
    };

    struct LocalTrajData
    {
        Trajectory traj;
        int traj_id;
        double duration;
        double start_time; // world time
        double end_time;   // world time
        Eigen::Vector2d start_pos;
        Eigen::Vector2d end_pos;
        double init_angle;
    };

    typedef std::vector<LocalTrajData> SingulTrajData;
    typedef std::vector<FlatTrajData> KinoTrajData;

    class TrajContainer
    {
    public:
        TrajContainer() {}
        ~TrajContainer() {}

        int traj_id = 0;
        SingulTrajData singul_traj;

        void addSingulTraj(const Trajectory &trajectory, const double &world_time)
        {
            LocalTrajData local_traj;
            local_traj.traj_id = ++traj_id;
            local_traj.duration = trajectory.getTotalDuration();
            local_traj.start_pos = trajectory.getJuncPos(0);
            local_traj.end_pos = trajectory.getJuncPos(trajectory.getPieceNum());
            local_traj.start_time = world_time;
            local_traj.end_time = world_time + local_traj.duration;
            local_traj.traj = trajectory;

            singul_traj.push_back(local_traj);
        }

        int locateSingulId(const double &t)
        {
            int number_of_singul_trajs = singul_traj.size();
            if (t < singul_traj[0].start_time)
            {
                return 0;
            }
            else if (t >= singul_traj[number_of_singul_trajs - 1].end_time)
            {
                return number_of_singul_trajs - 1;
            }
            for (int i = 0; i < number_of_singul_trajs; i++)
            {
                if (t >= singul_traj[i].start_time && t < singul_traj[i].end_time)
                    return i;
            }
            return -1;
        }

        void clearSingul()
        {
            singul_traj.clear();
            traj_id = 0;
        }
        
    };

} // namespace plan_utils

#endif // _TRAJ_CONTAINER_H_