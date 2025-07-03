#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include "circle_cover/circle_cover.hpp"

#include "map_utils/voxel_map.hpp"
#include "vehicle_utils/vehicle_geometry.hpp"
#include "path_searching/kino_astar.hpp"
#include "plan_manage/traj_planner.hpp"
#include "plan_utils/traj_container.hpp"

namespace py = pybind11;

// Debug print out
class PyStdoutStream : public std::streambuf
{
protected:
    virtual std::streamsize xsputn(const char *s, std::streamsize n) override
    {
        py::gil_scoped_acquire acquire;
        py::print(py::str(s, n), py::arg("end") = "");
        return n;
    }

    virtual int overflow(int c) override
    {
        if (c != EOF)
        {
            char z = c;
            xsputn(&z, 1);
        }
        return c;
    }
};

// VoxelMap Bindings
void createVoxelMapBindings(py::module &m)
{
    py::class_<map_utils::VoxelMap, std::shared_ptr<map_utils::VoxelMap>>(m, "VoxelMap")
        .def(py::init<const Eigen::Vector3d &, const Eigen::Vector3d &, double>(),
             py::arg("map_size"), py::arg("origin"), py::arg("resolution"))
        .def(py::init<const Eigen::Vector3i &, const Eigen::Vector3d &, double>(),
             py::arg("grid_size"), py::arg("origin"), py::arg("resolution"))
        .def("dilate", &map_utils::VoxelMap::dilate,
             py::arg("num"), py::arg("is_virtual") = false)
        .def("setOccupancy", (bool(map_utils::VoxelMap::*)(const Eigen::Vector3d &, uint8_t)) & map_utils::VoxelMap::setOccupancy)
        .def("setOccupancy", (bool(map_utils::VoxelMap::*)(const Eigen::Vector3i &, uint8_t)) & map_utils::VoxelMap::setOccupancy)
        .def("setOccupancy", (bool(map_utils::VoxelMap::*)(const Eigen::Vector2d &, uint8_t)) & map_utils::VoxelMap::setOccupancy)
        .def("setOccupancy", (bool(map_utils::VoxelMap::*)(const Eigen::Vector2i &, uint8_t)) & map_utils::VoxelMap::setOccupancy)
        .def("getResolution", &map_utils::VoxelMap::getResolution)
        .def("getOrigin", &map_utils::VoxelMap::getOrigin)
        .def("getMapSize", &map_utils::VoxelMap::getMapSize)
        .def("getGridSize", &map_utils::VoxelMap::getGridSize)
        .def("getVoxelState", (uint8_t(map_utils::VoxelMap::*)(const Eigen::Vector3d &)) & map_utils::VoxelMap::getVoxelState)
        .def("getVoxelState", (uint8_t(map_utils::VoxelMap::*)(const Eigen::Vector3i &)) & map_utils::VoxelMap::getVoxelState)
        .def("getVoxelState", (uint8_t(map_utils::VoxelMap::*)(const Eigen::Vector2d &)) & map_utils::VoxelMap::getVoxelState)
        .def("getVoxelState", (uint8_t(map_utils::VoxelMap::*)(const Eigen::Vector2i &)) & map_utils::VoxelMap::getVoxelState)
        .def("getVoxels", &map_utils::VoxelMap::getVoxels)
        .def("getSurfIdx", &map_utils::VoxelMap::getSurfIdx)
        .def("getSurf", [](map_utils::VoxelMap &self)
             {
            std::vector<Eigen::Vector3d> points;
            self.getSurf(points);
            return points; })
        .def("checkPointCollision", &map_utils::VoxelMap::checkPointCollision)
        .def("isValidState", &map_utils::VoxelMap::isValidState);
}

// circle cover Bindings
void circleCoverBindings(py::module &m)
{
    py::class_<geometry_cover::CoverCircleResult>(m, "CoverCircleResult")
        .def(py::init<>())
        .def_readwrite("radius", &geometry_cover::CoverCircleResult::radius)
        .def_readwrite("centers", &geometry_cover::CoverCircleResult::centers)
        .def_readwrite("aux_radius", &geometry_cover::CoverCircleResult::aux_radius)
        .def_readwrite("aux_centers", &geometry_cover::CoverCircleResult::aux_centers);
}

// VehicleGeometry Bindings
void createVehicleGeometryBindings(py::module &m)
{
    py::class_<vehicle_utils::VehicleGeometry, std::shared_ptr<vehicle_utils::VehicleGeometry>>(m, "VehicleGeometry")
        .def(py::init<const std::vector<Eigen::Vector2d> &, double>(),
             py::arg("vertices"), py::arg("car_wheelbase"))
        .def(py::init<const std::vector<Eigen::Vector2d> &, double, double>(),
             py::arg("vertices"), py::arg("car_wheelbase"), py::arg("max_sagitta"))
        .def("setVehicleVertices", &vehicle_utils::VehicleGeometry::setVehicleVertices)
        .def("setAuxVertices", [](vehicle_utils::VehicleGeometry &self, py::list vertices_python_list)
             {
            std::vector<std::vector<Eigen::Vector2d>> aux_vertices_list;
            for (auto vertices : vertices_python_list) {
                std::vector<Eigen::Vector2d> aux_vertices;
                for (auto v : vertices) {
                    aux_vertices.push_back(v.cast<Eigen::Vector2d>());
                }
                aux_vertices_list.push_back(aux_vertices);
            }
            self.setAuxVertices(aux_vertices_list); })
        .def("setMaxSagitta", &vehicle_utils::VehicleGeometry::setMaxSagitta)
        .def("getCarVertices", (const std::vector<Eigen::Vector2d> (vehicle_utils::VehicleGeometry::*)() const) & vehicle_utils::VehicleGeometry::getCarVertices)
        .def("getCarVertices", (std::vector<Eigen::Vector2d>(vehicle_utils::VehicleGeometry::*)(const Eigen::Vector3d &)) & vehicle_utils::VehicleGeometry::getCarVertices)
        .def("getAuxCount", &vehicle_utils::VehicleGeometry::getAuxCount)
        .def("getAllAuxVertices", (const std::vector<std::vector<Eigen::Vector2d>> (vehicle_utils::VehicleGeometry::*)() const) & vehicle_utils::VehicleGeometry::getAllAuxVertices)
        .def("getAllAuxVertices", (std::vector<std::vector<Eigen::Vector2d>>(vehicle_utils::VehicleGeometry::*)(const Eigen::Vector3d &)) & vehicle_utils::VehicleGeometry::getAllAuxVertices)
        .def("getAuxVertices", (const std::vector<Eigen::Vector2d> (vehicle_utils::VehicleGeometry::*)(int) const) & vehicle_utils::VehicleGeometry::getAuxVertices)
        .def("getAuxVertices", (std::vector<Eigen::Vector2d>(vehicle_utils::VehicleGeometry::*)(int, const Eigen::Vector3d &)) & vehicle_utils::VehicleGeometry::getAuxVertices)
        .def("getCoverCircle", &vehicle_utils::VehicleGeometry::getCoverCircle);
}

// KinoAstar Bindings
void createKinoAstarBindings(py::module &m)
{
    // KinoAstar class
    py::class_<path_searching::KinoAstar, std::shared_ptr<path_searching::KinoAstar>>(m, "KinoAstar")
        .def(py::init<const std::string &, map_utils::VoxelMap::Ptr &, vehicle_utils::VehicleGeometry::Ptr &>(),
             py::arg("config_string"), py::arg("map_ptr"), py::arg("vehicle_geo_ptr"))
        .def_readwrite("pose_list", &path_searching::KinoAstar::pose_list)
        .def_readwrite("control_input_list", &path_searching::KinoAstar::control_input_list)
        .def("checkVirtualCollisionAtPoint", &path_searching::KinoAstar::checkVirtualCollisionAtPoint)
        .def("checkCollisionUsingPosAndYaw", [](path_searching::KinoAstar &self, Eigen::Vector3d &state)
             {
                bool res;
                self.path_searching::KinoAstar::checkCollisionUsingPosAndYaw(state, res);
                return res; })
        .def("circleCheckCollisionUsingPosAndYaw", [](path_searching::KinoAstar &self, Eigen::Vector3d &state)
             {
            bool res;
            self.path_searching::KinoAstar::circleCheckCollisionUsingPosAndYaw(state, res);
            return res; })
        .def("circlePointTransformation", &path_searching::KinoAstar::circlePointTransformation)
        .def("getFullposeLists", &path_searching::KinoAstar::getFullposeLists)
        .def("getCoverCircle", &path_searching::KinoAstar::getCoverCircle)
        // .def("getSingulNodes", &path_searching::KinoAstar::getSingulNodes)
        .def("search", &path_searching::KinoAstar::search,
             py::arg("start_state"), py::arg("init_ctrl"),
             py::arg("end_state"), py::arg("use3d") = false)
        .def("searchTime", [](path_searching::KinoAstar &self, py::list flat_trajs_list)
             {
                plan_utils::KinoTrajData flat_trajs;
                bool result = self.searchTime(flat_trajs);
                flat_trajs_list.attr("clear")();
                for (const auto& item : flat_trajs)
                {
                    py::dict flat_traj_data;
                    flat_traj_data["singul"] = item.singul;
                    flat_traj_data["traj_pts"] = item.traj_pts;
                    flat_traj_data["thetas"] = item.thetas;
                    flat_traj_data["start_state"] = item.start_state;
                    flat_traj_data["final_state"] = item.final_state;
                    flat_trajs_list.append(flat_traj_data);
                }
                return result; })
        .def("getKinoNode", [](path_searching::KinoAstar &self, py::list flat_trajs_list)
             {
                plan_utils::KinoTrajData flat_trajs;
                self.getKinoNode(flat_trajs);
                flat_trajs_list.attr("clear")();
                for (const auto& item : flat_trajs)
                {
                    py::dict flat_traj_data;
                    flat_traj_data["singul"] = item.singul;
                    flat_traj_data["traj_pts"] = item.traj_pts;
                    flat_traj_data["thetas"] = item.thetas;
                    flat_traj_data["start_state"] = item.start_state;
                    flat_traj_data["final_state"] = item.final_state;
                    flat_trajs_list.append(flat_traj_data);
                } });
}

// Corridor Bindings
void CreateCorridorBindings(py::module &m)
{
    // Rectangle Corridor
    m.def("calcRectangleCorridor", &corridor_utils::calcRectangleCorridor,
          py::arg("hPolys"), py::arg("state_list"), py::arg("map_ptr"), py::arg("vehicle_geo_ptr"),
          py::arg("limitBound") = 10.0);
    m.def(
        "calcRectangleCorridor", [](const std::vector<Eigen::Vector3d> &state_list, const map_utils::VoxelMap &map, const vehicle_utils::VehicleGeometry &vehicle_geo, double limitBound) -> std::vector<Eigen::MatrixXd>
        {
        std::vector<Eigen::MatrixXd> hPolys;
        corridor_utils::calcRectangleCorridor(
            hPolys, 
            state_list, 
            std::make_shared<map_utils::VoxelMap>(map), 
            std::make_shared<vehicle_utils::VehicleGeometry>(vehicle_geo),
            limitBound);
        return hPolys; },
        py::arg("state_list"), py::arg("map"), py::arg("vehicle_geo"), py::arg("limitBound") = 10.0);
    m.def("calcAuxRectangleCorridor", &corridor_utils::calcAuxRectangleCorridor,
          py::arg("hPolys_list"), py::arg("state_list"), py::arg("map_ptr"), py::arg("vehicle_geo_ptr"),
          py::arg("limitBound") = 10.0);
    m.def(
        "calcAuxRectangleCorridor", [](const std::vector<Eigen::Vector3d> &state_list, const map_utils::VoxelMap &map, const vehicle_utils::VehicleGeometry &vehicle_geo, double limitBound) -> std::vector<std::vector<Eigen::MatrixXd>>
        {
        std::vector<std::vector<Eigen::MatrixXd>> hPolys_list;
        corridor_utils::calcAuxRectangleCorridor(
            hPolys_list, 
            state_list, 
            std::make_shared<map_utils::VoxelMap>(map), 
            std::make_shared<vehicle_utils::VehicleGeometry>(vehicle_geo),
            limitBound);
        return hPolys_list; },
        py::arg("state_list"), py::arg("map"), py::arg("vehicle_geo"), py::arg("limitBound") = 10.0);

    // FIRI Corridor
    m.def("calcFIRICorridor", &corridor_utils::calcFIRICorridor,
          py::arg("hPolys"), py::arg("state_list"), py::arg("map_ptr"), py::arg("vehicle_geo_ptr"), py::arg("range") = 10.0);
    m.def(
        "calcFIRICorridor", [](const std::vector<Eigen::Vector3d> &state_list, const map_utils::VoxelMap &map, const vehicle_utils::VehicleGeometry &vehicle_geo, double range) -> std::vector<Eigen::MatrixXd>
        {
        std::vector<Eigen::MatrixXd> hPolys;
        corridor_utils::calcFIRICorridor(
            hPolys,
            state_list, 
            std::make_shared<map_utils::VoxelMap>(map), 
            std::make_shared<vehicle_utils::VehicleGeometry>(vehicle_geo),
            range);

        return hPolys; },
        py::arg("state_list"), py::arg("map"), py::arg("vehicle_geo"), py::arg("range") = 10.0);
    m.def("calcAuxFIRICorridor", &corridor_utils::calcAuxFIRICorridor,
          py::arg("hPolys_list"), py::arg("state_list"), py::arg("map_ptr"), py::arg("vehicle_geo_ptr"), py::arg("range") = 10.0);
    m.def(
        "calcAuxFIRICorridor", [](const std::vector<Eigen::Vector3d> &state_list, const map_utils::VoxelMap &map, const vehicle_utils::VehicleGeometry &vehicle_geo, double range) -> std::vector<std::vector<Eigen::MatrixXd>>
        {
        std::vector<std::vector<Eigen::MatrixXd>> hPolys_list;
        corridor_utils::calcAuxFIRICorridor(
            hPolys_list,
            state_list, 
            std::make_shared<map_utils::VoxelMap>(map), 
            std::make_shared<vehicle_utils::VehicleGeometry>(vehicle_geo),
            range);

        return hPolys_list; },
        py::arg("state_list"), py::arg("map"), py::arg("vehicle_geo"), py::arg("range") = 10.0);
}

// TrajPlanner Bindings
void createTrajPlannerBindings(py::module &m)
{
    // Trajectory class
    py::class_<plan_utils::Trajectory>(m, "Trajectory")
        .def(py::init<>())
        .def("getSingul", &plan_utils::Trajectory::getSingul)
        .def("getDirection", &plan_utils::Trajectory::getDirection)
        .def("getPieceNum", &plan_utils::Trajectory::getPieceNum)
        .def("getDurations", &plan_utils::Trajectory::getDurations)
        .def("getTotalDuration", &plan_utils::Trajectory::getTotalDuration)
        .def("getTotalLength", &plan_utils::Trajectory::getTotalLength)
        .def("getPositions", &plan_utils::Trajectory::getPositions)
        .def("getPos", &plan_utils::Trajectory::getPos)
        .def("getR", &plan_utils::Trajectory::getR)
        .def("getRdot", &plan_utils::Trajectory::getRdot)
        .def("getdSigma", &plan_utils::Trajectory::getdSigma)
        .def("getddSigma", &plan_utils::Trajectory::getddSigma)
        .def("getVel", &plan_utils::Trajectory::getVel)
        .def("getAcc", &plan_utils::Trajectory::getAcc)
        .def("getLatAcc", &plan_utils::Trajectory::getLatAcc)
        .def("getAngle", &plan_utils::Trajectory::getAngle)
        .def("getCurv", &plan_utils::Trajectory::getCurv)
        .def("getOmega", &plan_utils::Trajectory::getOmega)
        .def("getSteer", &plan_utils::Trajectory::getSteer)
        .def("getSteerRate", &plan_utils::Trajectory::getSteerRate)
        .def("getJuncPos", &plan_utils::Trajectory::getJuncPos)
        .def("getJuncdSigma", &plan_utils::Trajectory::getJuncdSigma)
        .def("getJuncddSigma", &plan_utils::Trajectory::getJuncddSigma);

    // TrajPlanner class
    py::class_<plan_manage::TrajPlanner>(m, "TrajPlanner")
        .def(py::init<const std::string &, map_utils::VoxelMap::Ptr &, vehicle_utils::VehicleGeometry::Ptr &>(),
             py::arg("config_string"), py::arg("map_ptr"), py::arg("vehicle_geo_ptr"))
        .def("generateKinoPath", &plan_manage::TrajPlanner::generateKinoPath,
             py::arg("start"), py::arg("end"))
        .def("RunMINCO", &plan_manage::TrajPlanner::RunMINCO)
        .def("checkCollisionUsingPosAndYaw", &plan_manage::TrajPlanner::checkCollisionUsingPosAndYaw)
        .def("getKinoAstarPathFinder", &plan_manage::TrajPlanner::getKinoAstarPathFinder)
        .def("getFlatTraj", [](plan_manage::TrajPlanner &self)
             {
            py::list flat_trajs_list;
            for (const auto& item : self.getFlatTraj())
            {
                py::dict flat_traj_data;
                flat_traj_data["singul"] = item.singul;
                flat_traj_data["traj_pts"] = item.traj_pts;
                flat_traj_data["thetas"] = item.thetas;
                flat_traj_data["start_state"] = item.start_state;
                flat_traj_data["final_state"] = item.final_state;
                flat_trajs_list.append(flat_traj_data);
            }

            return flat_trajs_list; })
        .def("getCorridorState", &plan_manage::TrajPlanner::getCorridorState)
        .def("getCorridor", &plan_manage::TrajPlanner::getCorridor)
        .def("getAuxCorridor", &plan_manage::TrajPlanner::getAuxCorridor)
        .def("getOptTraj", [](plan_manage::TrajPlanner &self)
             {
            py::list opt_trajs_list;
            for (const auto& item : self.getOptTraj())
            {
                py::dict opt_traj_data;
                opt_traj_data["duration"] = item.duration;
                opt_traj_data["start_time"] = item.start_time;
                opt_traj_data["end_time"] = item.end_time;
                opt_traj_data["start_pos"] = item.start_pos;
                opt_traj_data["end_pos"] = item.end_pos;
                opt_traj_data["traj"] = item.traj;
                opt_trajs_list.append(opt_traj_data);
            }

            return opt_trajs_list; })
        .def("getConstraintPts", &plan_manage::TrajPlanner::getConstraintPts)
        .def("getCostsHistory", &plan_manage::TrajPlanner::getCostsHistory)
        // debug
        .def("getIniStates", &plan_manage::TrajPlanner::getIniStates)
        .def("getFinStates", &plan_manage::TrajPlanner::getFinStates)
        .def("getInitInnerPts", &plan_manage::TrajPlanner::getInitInnerPts)
        .def("getInitTs", &plan_manage::TrajPlanner::getInitTs)
        .def("getSinguls", &plan_manage::TrajPlanner::getSinguls);
}

PYBIND11_MODULE(pymincocar, m)
{
    m.doc() = " ";

    static PyStdoutStream pyStdoutStream;
    std::cout.rdbuf(&pyStdoutStream);

    createVoxelMapBindings(m);
    createVehicleGeometryBindings(m);
    circleCoverBindings(m);
    createKinoAstarBindings(m);
    CreateCorridorBindings(m);
    createTrajPlannerBindings(m);
}