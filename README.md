# Efficient And Safe Trajectory (EAST) Planner for Autonomous Agricultural Vehicle with Implements

[![IEEE Xplore](https://img.shields.io/badge/IEEE%20Xplore-10854653-blue)](https://ieeexplore.ieee.org/document/10854653)
[![arXiv VAMP](https://img.shields.io/badge/arXiv-2501.10636-b31b1b.svg)](https://arxiv.org/abs/2501.10636)

[EAST-planner project website](https://agroboticsresearch.github.io/east_planner/)

This is the project page of the paper "**Efficient and Safe Trajectory Planning for Autonomous Agricultural Vehicle Headland Turning in Cluttered Orchard Environments**". This repo introduces a new method of efficient trajectory planning for agricultural vehicle on the headland. The environment, as well as the agricultural vehicle, are modeled with convex polygons. The planning algorithm is implemented with C++ and the main functions were visualized with Python leveraging Pybind modules.

If you found this work useful for your own research, please consider citing the following:
```bibtex
@article{Wei_2025,
   title={Efficient and Safe Trajectory Planning for Autonomous Agricultural Vehicle Headland Turning in Cluttered Orchard Environments}, 
   author={Wei, Peng and Peng, Chen and Lu, Wenwu and Zhu, Yuankai and Vougioukas, Stavros and Fei, Zhenghao and Ge, Zhikang},
   journal={IEEE Robotics and Automation Letters}, 
   year={2025},
   volume={10},
   number={3},
   pages={2574-2581},
   doi={10.1109/LRA.2025.3534056}
}
```

If you have any question, please feel free to contact us via `penwei@ucdavis.edu (Peng Wei)` and `chen.peng@zju.edu.cn (Chen Peng)`.

The main contribution of the proposed methodology is as follows:

1. An efficient way of modeling the vehicle's configuration space and implemented a fast way to do the collision checking.
2. Use combinatoral safe corridors to express the manuvering space of the vehicle, as well as its added on implements. 
3. A differential based method was applied to solve the path planning problem online.

### Front end searching

1. The vehicle (w/o implement) is modeled with rectangles covered by circles. 
<p align="center">
   <img src="images/circle_cover_and_max.png" width="400"/>
</p>

2. Expanding the polygon geometry
<p align="center">
   <img src="images/config_space.png" width="500"/>
</p>

4. Searching result for a test case
<p align="center">
   <img src="images/front_end_search.png" width="500"/>
</p>

### Back end smoothing

For the back end smoothing, we applied combined corridors for represent the feasibility constraints. In this way, the path can be smoothed in the constrained area.

<p align="center">
   <img src="images/combined_corridors.png" width="400"/>
   <img src="images/combined_corridor2.png" width="280"/>
   <img src="images/kms_sprayer.png" width="480"/>
</p>
