# Efficient and Safe Trajectory (EAST) Planner for Autonomous Agricultural Vehicle Headland Turning in Cluttered Orchard Environments 

[![IEEE Xplore](https://img.shields.io/badge/IEEE%20Xplore-10854653-blue)](https://ieeexplore.ieee.org/document/10854653)
[![arXiv VAMP](https://img.shields.io/badge/arXiv-2501.10636-b31b1b.svg)](https://arxiv.org/abs/2501.10636)
[![Project Website](https://img.shields.io/badge/Project_Page-green)](https://agroboticsresearch.github.io/east_planner/)


## Updates
- **[July 3, 2025]** We have released our codebase! 

- **[January 27, 2025]** Our paper "[Efficient and Safe Trajectory Planning for Autonomous Agricultural Vehicle Headland Turning in Cluttered Orchard Environments](https://ieeexplore.ieee.org/document/10854653)" has been published in *IEEE Robotics and Automation Letters (RA-L)*. 


## Overview
This repository introduces a novel method for efficient and safe trajectory planning for agricultural vehicles performing headland turning in cluttered orchard environments. Both the environment and the vehicle (including its implement) are modeled using convex polygons. The planning algorithm is implemented in C++, and the Python bindings of the core functions are provided using pybind11.

### Main Contributions

1. An efficient way of modeling the vehicle's configuration space and implemented a fast way to do the collision checking.
2. Use combinatoral safe corridors to express the manuvering space of the vehicle, as well as its added on implements. 
3. A differential based method was applied to solve the path planning problem online.

### Front-End Searching

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

### Back-End Smoothing

For the back end smoothing, we applied combined corridors for represent the feasibility constraints. In this way, the path can be smoothed in the constrained area.

<p align="center">
   <img src="images/combined_corridors.png" width="400"/>
   <img src="images/combined_corridor2.png" width="280"/>
   <img src="images/kms_sprayer.png" width="480"/>
</p>

If you have any questions, feel free to contact us: [Peng Wei](mailto:penwei@ucdavis.edu) and [Chen Peng](mailto:chen.peng@zju.edu.cn).

## Installation Guide

Here's how to install the project, following the steps you provided:

### 1\. Clone the repository

```bash
git clone https://github.com/AgRoboticsResearch/east_planner.git
```

### 2\. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 3\. Initialize Pybind11 submodule

```bash
git submodule update --init
```

### 4\. Install system dependencies

```bash
# Update package lists
sudo apt update
# Install Eigen3 and OMPL development libraries
sudo apt install libeigen3-dev libompl-dev
# Install YAML-CPP development library
sudo apt install libyaml-cpp-dev
# Install essential build tools
sudo apt install build-essential cmake libssl-dev python3-pip
```

### 5\. Install scikit-build

```bash
pip3 install scikit-build
```

### 6\. Install MincoCar Python package

```bash
# Navigate back to the MincoCar directory
cd MincoCar
pip3 install .
```

### 7\. Verify installation
After installation, you can verify that everything works correctly by running an example:
```bash
# Navigate to the example directory
cd example/
# Run the provided example script
python3 east_planner_example.py
```
This will run a basic demonstration of the planner. For more comprehensive examples and tests, please explore the MincoCar/test/ directory.

## BibTex
If you find this work useful for your own research, please cite the following:

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

## Acknowledgment
We gratefully acknowledge the use of the following open-source tools:
- [pybind11](https://github.com/pybind/pybind11) —  for Python bindings of C++ modules. 
- [qpOASES](https://github.com/coin-or/qpOASES) — for solving quadratic programming (QP) problems
- [Dftpav](https://github.com/ZJU-FAST-Lab/Dftpav/tree/main) — for building safety corridor

## License
This project is released under the Apache License 2.0 License.