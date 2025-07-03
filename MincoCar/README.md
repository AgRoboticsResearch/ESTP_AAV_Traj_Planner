# MincoCar

Python binding of the [Dftpav](https://github.com/ZJU-FAST-Lab/Dftpav) and [Car-like-Robotic-swarm](https://github.com/ZJU-FAST-Lab/Car-like-Robotic-swarm) packages.

## Dependency
```bash
# pybind11
cd <root folder>
git submodule update --init

# other
sudo apt update
sudo apt install libeigen3-dev libompl-dev
sudo apt install libyaml-cpp-dev
sudo apt install build-essential cmake libssl-dev python3-pip

# install scikit-build
pip3 install scikit-build
```

## To install **pymincocar**
```bash
cd MincoCar
pip3 install .
```

## Examples
Try the python examples in the [test](test/) folder.