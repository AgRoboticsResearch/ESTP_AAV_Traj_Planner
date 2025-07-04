# Only autonomous agricultural vehicle (AAV)
import os, sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import yaml
import pickle

from MincoCar.utils.voxel_map_utils import (
    FeatureMap2D,
    create_voxel_map_from_grid_map,
    draw_voxel_map_in_2d,
)
from MincoCar.utils.vehicle_utils import VehicleModel
from MincoCar.utils.polygon_utils import *
from MincoCar.utils.minco_utils import *

from pymincocar import VehicleGeometry, TrajPlanner

COLLISION_CHECK_TYPE = 1  # 0 use raycaster,1 use circle
USE_FIRI = False


def plot_circle_result(ax, circle_result):
    for center in circle_result.centers:
        circle = patches.Circle(
            center, circle_result.radius, edgecolor="b", facecolor="none"
        )
        ax.add_patch(circle)
    ax.set_aspect("equal")


def plot_kino_poses(ax, poses):
    for pose in poses:
        x, y, yaw = pose[:3]
        rect = patches.Rectangle(
            (x, y), 1, 0.5, angle=np.rad2deg(yaw), edgecolor="r", facecolor="none"
        )
        ax.add_patch(rect)
    ax.set_aspect("equal")


if __name__ == "__main__":

    """Generate a voxel map"""
    pkl_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "MincoCar", "test", "headland.pkl")
    with open(pkl_path, "rb") as handle:
        obstacles = pickle.load(handle)

    polys = [Polygon(obstacle) for obstacle in obstacles]

    # Create a feature map
    feature_map = FeatureMap2D(x_min=-15, x_max=15, y_min=-5, y_max=15)
    feature_map.add_polygons(polys)

    # Plot feature map
    fig, ax = plt.subplots(1, 1)
    feature_map.draw(ax)

    # Generate grid map
    resolution = 0.1  # m
    grid_map = feature_map.to_gridmap(resolution, use_raster=True)

    # Generate voxel map
    voxel_map = create_voxel_map_from_grid_map(grid_map)
    if USE_FIRI:
        voxel_map.dilate(1, is_virtual=False)

    # Vehicle geometry model
    axle_to_front = 2.85
    axle_to_back = 0.5
    car_width = 1.48
    car_wheelbase = 1.9
    car_max_steer = 31.5  # deg
    vehicle_model = VehicleModel(
        max_steer=np.deg2rad(car_max_steer),
        wheel_base=car_wheelbase,
        width=car_width,
        axle_to_front=axle_to_front,
        axle_to_back=axle_to_back,
    )
    if COLLISION_CHECK_TYPE == 0:
        convex_vertices = vehicle_model.get_conservative_convex_poly_vertices_in_body()
    else:
        convex_vertices = vehicle_model.get_car_convex_poly_vertices_in_body() #1
    vehicle_geo = VehicleGeometry(convex_vertices, car_wheelbase)

    """ Init Trajectory Planner """
    # config parameters
    config_dict = {
        "search": {
            "yaw_resolution": 0.3,
            "lambda_heu": 5.0,
            "allocate_num": 100000,
            "check_num": 5,
            "max_search_time": 5000.0,
            "traj_forward_penalty": 1.0,
            "traj_back_penalty": 2.5,
            "traj_gear_switch_penalty": 15.0,
            "traj_steer_penalty": 0.5,
            "traj_steer_change_penalty": 0.0,
            "step_arc": 1.0,
            "step_steer_res": 0.5,
            "checkl": 0.5,
            "max_vel": 1.5,
            "max_acc": 1.0,
            "max_cur": 0.3,
            "max_steer": 31.5,  # deg
            "time_resolution": 0.1,
            "distance_resolution": 0.5,
            "velocity_resolution": 0.5,
            "non_siguav": 0.05,
            "collision_check_type": COLLISION_CHECK_TYPE,
        },
        "optimization": {
            "use_firi": USE_FIRI,
            "traj_piece_duration": 1.0,
            "traj_resolution": 10,
            "dense_traj_resolution": 20,
            "wei_sta_obs": 2000.0,
            "wei_feas": 2000.0,
            "wei_time": 100.0,
            "max_vel": 1.5,
            "max_acc": 1.0,
            "max_cur": 0.350877,
            "max_omg": 1.0,
            "lbfgs_delta": 1.0e-10,
            "gear_opt": True,
            "log_cost": True,
            "verbose": True,
        },
    }

    # Circle check
    if COLLISION_CHECK_TYPE == 1:
        safety_distance = 0.35
        row_width = 2.5
        max_sagitta = row_width / 2 - car_width / 2 - safety_distance
        vehicle_geo.setMaxSagitta(max_sagitta)
        circle_result = vehicle_geo.getCoverCircle()

        voxel_map.dilate(int(circle_result.radius / resolution) + 1, is_virtual=True)

    # Draw voxel map
    draw_voxel_map_in_2d(ax, voxel_map)

    # Init TrajPlanner
    config_yaml = yaml.dump(config_dict)
    traj_planner = TrajPlanner(config_yaml, voxel_map, vehicle_geo)

    # Assign start and end state
    # state = [x, y, yaw (rad), vel]
    start_state = np.array([1.54286108, 8.75, -np.pi, 0])
    end_state = np.array([0.66122618, 3.75, 0, 0])

    """ Generate Kino Astar path """
    # Front end search
    front_end_success = traj_planner.generateKinoPath(start_state, end_state)
    if not front_end_success:
        exit("Exit!")

    kino_search_poses = np.array(traj_planner.getKinoAstarPathFinder().pose_list)
    ax.plot(
        kino_search_poses[:, 0],
        kino_search_poses[:, 1],
        "b-o",
        markersize=3,
        alpha=0.5,
        label="Kino search",
    )

    # Get flat trajectories
    flat_trajs = traj_planner.getFlatTraj()
    (
        kino_searchTime_poses,
        kino_searchTime_time,
        kino_searchTime_dir,
    ) = extract_from_flat_trajs(flat_trajs)

    # Plot flat trajectories
    ax.plot(
        kino_searchTime_poses[:, 0],
        kino_searchTime_poses[:, 1],
        "r-s",
        markersize=1,
        alpha=0.5,
        label="Kino searchTime",
    )

    fig2, ax2 = plt.subplots(8, 1, sharex=True, figsize=(6, 6))
    ax2[0].set_ylabel("x")
    ax2[1].set_ylabel("y")
    ax2[2].set_ylabel("theta")
    ax2[3].set_ylabel("dir")
    ax2[4].set_ylabel("vel")
    ax2[5].set_ylabel("acc")
    ax2[6].set_ylabel("curv")
    ax2[7].set_ylabel("omg")
    ax2[-1].set_xlabel("time")

    ax2[4].axhline(
        config_dict["optimization"]["max_vel"], color="k", linestyle="--", linewidth=0.5
    )
    ax2[4].axhline(
        -config_dict["optimization"]["max_vel"],
        color="k",
        linestyle="--",
        linewidth=0.5,
    )
    ax2[5].axhline(
        config_dict["optimization"]["max_acc"], color="k", linestyle="--", linewidth=0.5
    )
    ax2[5].axhline(
        -config_dict["optimization"]["max_acc"],
        color="k",
        linestyle="--",
        linewidth=0.5,
    )
    ax2[6].axhline(
        config_dict["optimization"]["max_cur"], color="k", linestyle="--", linewidth=0.5
    )
    ax2[6].axhline(
        -config_dict["optimization"]["max_cur"],
        color="k",
        linestyle="--",
        linewidth=0.5,
    )
    ax2[7].axhline(
        config_dict["optimization"]["max_omg"], color="k", linestyle="--", linewidth=0.5
    )
    ax2[7].axhline(
        -config_dict["optimization"]["max_omg"],
        color="k",
        linestyle="--",
        linewidth=0.5,
    )

    for i in range(3):
        if i == 2:
            ax2[i].plot(
                kino_searchTime_time,
                normalize_angle(kino_searchTime_poses[:, i]),
                "ro",
                markersize=1,
            )
        else:
            ax2[i].plot(
                kino_searchTime_time, kino_searchTime_poses[:, i], "ro", markersize=1
            )

    ax2[3].plot(
        kino_searchTime_time,
        kino_searchTime_dir,
        "ro",
        markersize=1,
        label="Kino searchTime",
    )

    """ Run MINCO """
    # Back end optimization
    back_end_success = traj_planner.RunMINCO()
    if back_end_success:
        # Visualize safe convex corridors
        plot_corridors(ax, traj_planner.getCorridor(), fill=False)

        # Extract optimal trajectories
        opt_trajs_list = traj_planner.getOptTraj()
        sample_ts, sample_poses, sample_v_a_cur_omg, sample_dirs = (
            extract_from_opt_trajs(opt_trajs_list, 0.1)
        )
        
        # Plot optimal trajectories
        for pose in sample_poses:
            vehicle_model.draw_car(
                ax, pose[0], pose[1], pose[2], car_color="C1", aux_color="C0"
            )
            if traj_planner.checkCollisionUsingPosAndYaw(pose):
                print(f"Collision detected at: {pose}")

        ax.plot(
            sample_poses[:, 0],
            sample_poses[:, 1],
            "g-*",
            markersize=1,
            alpha=0.5,
            label="OptTraj",
        )

        for i in range(3):
            ax2[i].plot(sample_ts, sample_poses[:, i], "go", markersize=1)

        for i in range(4):
            ax2[i + 4].plot(sample_ts, sample_v_a_cur_omg[:, i], "go", markersize=1)

        ax2[3].plot(
            sample_ts,
            sample_dirs,
            "go",
            markersize=1,
            label="OptTraj",
        )

        if config_dict["optimization"]["log_cost"]:
            costs_history = traj_planner.getCostsHistory()
            costs_history = np.array(costs_history)

            fig3, ax3 = plt.subplots(1, 1)
            ax3.grid(True)
            ax3.set_xlabel("iteration")
            ax3.set_ylabel("log(cost)")
            ax3.set_title("Costs History")
            epis = 1.0e-6
            ax3.plot(
                costs_history[:, 0],
                np.log10(costs_history[:, 1] + epis),
                label="smooth",
            )
            ax3.plot(
                costs_history[:, 0], np.log10(costs_history[:, 2] + epis), label="time"
            )
            ax3.plot(
                costs_history[:, 0],
                np.log10(costs_history[:, 3] + epis),
                label="obstacle",
            )
            ax3.plot(
                costs_history[:, 0],
                np.log10(costs_history[:, 4] + epis),
                label="feasibility",
            )
            ax3.legend(fontsize=9)

    ax.legend(fontsize=9)
    ax2[3].legend(fontsize=9)

    fig.tight_layout()
    fig2.tight_layout()
    plt.show()
