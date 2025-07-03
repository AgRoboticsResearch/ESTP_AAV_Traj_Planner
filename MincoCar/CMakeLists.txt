cmake_minimum_required(VERSION 3.12)
project(pymincocar)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Eigen3 REQUIRED)
find_package(ompl REQUIRED)
find_package(yaml-cpp REQUIRED)

# external pybind11 library
set(PYBIND11_PATH ${PROJECT_SOURCE_DIR}/../extern/pybind11)
add_subdirectory(${PYBIND11_PATH} ${PROJECT_BINARY_DIR}/pybind11)

# external qpOASES library
set(qpOASES_PATH ${PROJECT_SOURCE_DIR}/../extern/qpOASES)
add_subdirectory(${qpOASES_PATH} ${PROJECT_BINARY_DIR}/qpOASES)

include_directories(
    include
    ${EIGEN3_INCLUDE_DIR}
    ${OMPL_INCLUDE_DIRS}
    ${YAML_CPP_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/../extern/qpOASES/include
)

# Solve the issue to compile on macOS
set(DYLD_LIBRARY_PATH $ENV{DYLD_LIBRARY_PATH})
if (DYLD_LIBRARY_PATH)
    string(REPLACE ":" ";" DYLD_LIBRARY_PATH_LIST ${DYLD_LIBRARY_PATH})
    foreach(LIB_PATH IN LISTS DYLD_LIBRARY_PATH_LIST)
        link_directories(${LIB_PATH})
    endforeach()
endif()

# map_utils
add_library(map_utils
    src/map_utils/voxel_map.cpp
    src/map_utils/raycast.cpp
)

# corridor_utils
add_library(corridor_utils
    src/corridor_utils/rectangle_corridor.cpp
    src/corridor_utils/firi_corridor.cpp
)
target_link_libraries(corridor_utils PRIVATE qpOASES)

# path_searching
add_library(path_searching
    src/path_searching/kino_astar.cpp
)
target_link_libraries(path_searching PRIVATE
    map_utils
    circle_cover
    ${OMPL_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
)

# circle_cover
add_library(circle_cover
    src/circle_cover/circle_cover.cpp
)
target_link_libraries(circle_cover PRIVATE
    map_utils
    corridor_utils
    path_searching
    ${OMPL_LIBRARIES}
    ${YAML_CPP_LIBRARIES}
)

# plan_manage
add_library(plan_manage
    src/plan_manage/traj_planner.cpp
    src/plan_manage/traj_optimizer.cpp
)
target_link_libraries(plan_manage PRIVATE
    map_utils
    corridor_utils
    path_searching
    circle_cover
    ${YAML_CPP_LIBRARIES}
)

# pybind
pybind11_add_module(pymincocar src/pybind/mincocar_pybind.cpp)
target_link_libraries(pymincocar PRIVATE
    plan_manage
)

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION .)