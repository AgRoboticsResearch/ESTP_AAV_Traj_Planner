from typing import List, Union
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.spatial import ConvexHull


def normalize_angle(angles):
    # output angle between -pi to pi
    wrap_angle = (angles + math.pi) % (2 * math.pi) - math.pi
    return wrap_angle


def extract_from_flat_trajs(flat_traj_list: List):
    basetime = 0.0
    poses = []
    times = []
    dirs = []
    for idx, flat_traj in enumerate(flat_traj_list):
        singul = flat_traj["singul"]
        start_state = flat_traj["start_state"]
        final_state = flat_traj["final_state"]
        init_length = len(poses)
        # pose
        if idx == 0:
            poses.append(
                [
                    start_state[0][0],
                    start_state[1][0],
                    np.arctan2(singul * start_state[1][1], singul * start_state[0][1]),
                ]
            )
            times.append(basetime)

        poses.extend(
            [
                [pts[0], pts[1], theta]
                for pts, theta in zip(flat_traj["traj_pts"], flat_traj["thetas"])
            ]
        )
        # time
        times.extend(
            (np.cumsum(np.array(flat_traj["traj_pts"])[:, 2]) + basetime).tolist()
        )
        basetime = times[-1]
        # dir
        added_length = len(poses) - init_length
        dirs.extend(singul * np.ones(added_length))

    poses = np.array(poses)
    times = np.array(times).flatten()
    dirs = np.array(dirs).flatten()

    return poses, times, dirs

def extract_from_opt_trajs(opt_trajs_list: List, dt=0.1, wheelbase=1.9):
    basetime = 0.0
    sample_ts = []
    sample_poses = []
    sample_v_a_cur_omg_steer = []
    sample_dirs = []
    for idx, opt_traj in enumerate(opt_trajs_list):
        duration = opt_traj["duration"]
        trajectory = opt_traj["traj"]
        dir = trajectory.getDirection()
        sample_t = np.arange(0.0, duration, dt)
        if sample_t[-1] != duration:
            sample_t = np.append(sample_t, duration)
        # pose and v_a_cur
        for t in sample_t:
            t = min(t, duration)
            pos = trajectory.getPos(t)
            angle = trajectory.getAngle(t)
            angle = normalize_angle(angle)
            sample_poses.append([pos[0], pos[1], angle])
            sample_v_a_cur_omg_steer.append(
                [
                    trajectory.getVel(t),
                    trajectory.getAcc(t),
                    trajectory.getCurv(t),
                    trajectory.getOmega(t),
                    trajectory.getSteerRate(t, wheelbase),
                ]
            )
        # time
        sample_ts.extend(sample_t + basetime)
        basetime += duration
        # dir
        sample_dirs.extend(dir * np.ones_like(sample_t))

    sample_poses = np.array(sample_poses)
    sample_v_a_cur_omg_steer = np.array(sample_v_a_cur_omg_steer)
    sample_ts = np.array(sample_ts).flatten()
    sample_dirs = np.array(sample_dirs).flatten()

    return sample_ts, sample_poses, sample_v_a_cur_omg_steer, sample_dirs

def plot_corridors(ax: plt.Axes, hpolys_list: Union[np.ndarray, List], color="C0", fill=True, alpha=0.1):
    if isinstance(hpolys_list, np.ndarray):
        plot_corridor(ax, hpolys_list, color, fill, alpha)
        return

    if isinstance(hpolys_list, list):
        for item in hpolys_list:
            plot_corridors(ax, item, color, fill, alpha)


def plot_corridor(ax: plt.Axes, hpoly: np.ndarray, color="C0", fill=True, alpha=0.1):
    if not hpoly.any():
        print("corridor is empty!")
        return
    
    vertices = hpoly[2:, :].T
    hull = ConvexHull(vertices)
    hull_points = vertices[hull.vertices]
    polygon = patches.Polygon(
        hull_points,
        closed=True,
        facecolor=color if fill else "none",
        edgecolor="none" if fill else color,
        linewidth=1.0,
        alpha=alpha,
    )
    ax.add_patch(polygon)
