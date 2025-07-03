#include "corridor_utils/safe_corridor.hpp"

namespace corridor_utils
{

    void calcRectangleCorridor(std::vector<Eigen::MatrixXd> &hPolys,
                               const std::vector<Eigen::Vector3d> &state_list,
                               const map_utils::VoxelMap::Ptr &map_ptr,
                               const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                               double limitBound)
    {
        hPolys.clear();
        hPolys.reserve(state_list.size());
        double resolution = map_ptr->getResolution();
        double step = resolution * 1.0;

        double minX = vehicle_geo_ptr->getMinX();
        double minY = vehicle_geo_ptr->getMinY();
        double maxX = vehicle_geo_ptr->getMaxX();
        double maxY = vehicle_geo_ptr->getMaxY();

        for (const auto &state : state_list)
        {
            Eigen::Vector2d pos = state.head(2);
            double yaw = state(2);
            Eigen::Matrix2d ego_R;
            ego_R << cos(yaw), -sin(yaw),
                sin(yaw), cos(yaw);

            Eigen::Vector2d corners[4] = {
                {minX, maxY}, // top_left
                {maxX, maxY}, // top_right
                {maxX, minY}, // bottom_right
                {minX, minY}  // bottom_left
            };

            Eigen::Vector4d is_expanding(1.0, 1.0, 1.0, 1.0);
            Eigen::Vector2d point1, point2, new_point1, new_point2;

            while (is_expanding.any())
            {
                for (int i = 0; i < 4; ++i)
                {
                    if (!is_expanding[i])
                        continue;

                    bool isCollision = false;
                    switch (i)
                    {
                    case 0: // top
                        point1 = pos + ego_R * corners[0];
                        point2 = pos + ego_R * corners[1];
                        new_point1 = pos + ego_R * (corners[0] + Eigen::Vector2d(0, step));
                        new_point2 = pos + ego_R * (corners[1] + Eigen::Vector2d(0, step));
                        break;
                    case 1: // right
                        point1 = pos + ego_R * corners[1];
                        point2 = pos + ego_R * corners[2];
                        new_point1 = pos + ego_R * (corners[1] + Eigen::Vector2d(step, 0));
                        new_point2 = pos + ego_R * (corners[2] + Eigen::Vector2d(step, 0));
                        break;
                    case 2: // bottom
                        point1 = pos + ego_R * corners[2];
                        point2 = pos + ego_R * corners[3];
                        new_point1 = pos + ego_R * (corners[2] + Eigen::Vector2d(0, -step));
                        new_point2 = pos + ego_R * (corners[3] + Eigen::Vector2d(0, -step));
                        break;
                    case 3: // left
                        point1 = pos + ego_R * corners[3];
                        point2 = pos + ego_R * corners[0];
                        new_point1 = pos + ego_R * (corners[3] + Eigen::Vector2d(-step, 0));
                        new_point2 = pos + ego_R * (corners[0] + Eigen::Vector2d(-step, 0));
                        break;
                    }

                    isCollision = map_ptr->checkLineCollision(point1, new_point1) ||
                                  map_ptr->checkLineCollision(new_point1, new_point2) ||
                                  map_ptr->checkLineCollision(new_point2, point2);

                    if (isCollision)
                    {
                        is_expanding[i] = 0.0;
                    }
                    else
                    {
                        switch (i)
                        {
                        case 0: // top
                            corners[0][1] += step;
                            corners[1][1] += step;
                            if (abs(corners[0][1] - maxY) > limitBound)
                                is_expanding[i] = 0.0;
                            break;
                        case 1: // right
                            corners[1][0] += step;
                            corners[2][0] += step;
                            if (abs(corners[1][0] - maxX) > limitBound)
                                is_expanding[i] = 0.0;
                            break;
                        case 2: // bottom
                            corners[2][1] -= step;
                            corners[3][1] -= step;
                            if (abs(corners[2][1] - minY) > limitBound)
                                is_expanding[i] = 0.0;
                            break;
                        case 3: // left
                            corners[3][0] -= step;
                            corners[0][0] -= step;
                            if (abs(corners[3][0] - minX) > limitBound)
                                is_expanding[i] = 0.0;
                            break;
                        }
                    }
                }
            }

            Eigen::Vector2d new_point, norm;
            Eigen::MatrixXd hPoly(4, 4);

            new_point = pos + ego_R * corners[0];
            norm << -cos(yaw), -sin(yaw);
            hPoly.col(0).head<2>() = norm;
            hPoly.col(0).tail<2>() = new_point;

            new_point = pos + ego_R * corners[1];
            norm << -sin(yaw), cos(yaw);
            hPoly.col(1).head<2>() = norm;
            hPoly.col(1).tail<2>() = new_point;

            new_point = pos + ego_R * corners[2];
            norm << cos(yaw), sin(yaw);
            hPoly.col(2).head<2>() = norm;
            hPoly.col(2).tail<2>() = new_point;

            new_point = pos + ego_R * corners[3];
            norm << sin(yaw), -cos(yaw);
            hPoly.col(3).head<2>() = norm;
            hPoly.col(3).tail<2>() = new_point;

            hPolys.push_back(hPoly);
        }
    }

    void calcAuxRectangleCorridor(std::vector<std::vector<Eigen::MatrixXd>> &hPolys_list,
                                  const std::vector<Eigen::Vector3d> &state_list,
                                  const map_utils::VoxelMap::Ptr &map_ptr,
                                  const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                                  double limitBound)
    {
        hPolys_list.clear();
        hPolys_list.resize(vehicle_geo_ptr->getAuxCount());
        if (vehicle_geo_ptr->getAuxCount() == 0)
            return;

        double resolution = map_ptr->getResolution();
        double step = resolution * 1.0;

        double minX, maxX, minY, maxY;
        for (size_t k = 0; k < vehicle_geo_ptr->getAuxCount(); k++)
        {
            std::tuple<double, double, double, double> bbox =
                vehicle_geo_ptr->calculateBoundingBox(vehicle_geo_ptr->getAuxVertices(k));
            minX = std::get<0>(bbox);
            maxX = std::get<1>(bbox);
            minY = std::get<2>(bbox);
            maxY = std::get<3>(bbox);

            for (const auto &state : state_list)
            {
                Eigen::Vector2d pos = state.head(2);
                double yaw = state(2);
                Eigen::Matrix2d ego_R;
                ego_R << cos(yaw), -sin(yaw),
                    sin(yaw), cos(yaw);

                Eigen::Vector2d corners[4] = {
                    {minX, maxY}, // top_left
                    {maxX, maxY}, // top_right
                    {maxX, minY}, // bottom_right
                    {minX, minY}  // bottom_left
                };

                Eigen::Vector4d is_expanding(1.0, 1.0, 1.0, 1.0);
                Eigen::Vector2d point1, point2, new_point1, new_point2;

                while (is_expanding.any())
                {
                    for (int i = 0; i < 4; ++i)
                    {
                        if (!is_expanding[i])
                            continue;

                        bool isCollision = false;
                        switch (i)
                        {
                        case 0: // top
                            point1 = pos + ego_R * corners[0];
                            point2 = pos + ego_R * corners[1];
                            new_point1 = pos + ego_R * (corners[0] + Eigen::Vector2d(0, step));
                            new_point2 = pos + ego_R * (corners[1] + Eigen::Vector2d(0, step));
                            break;
                        case 1: // right
                            point1 = pos + ego_R * corners[1];
                            point2 = pos + ego_R * corners[2];
                            new_point1 = pos + ego_R * (corners[1] + Eigen::Vector2d(step, 0));
                            new_point2 = pos + ego_R * (corners[2] + Eigen::Vector2d(step, 0));
                            break;
                        case 2: // bottom
                            point1 = pos + ego_R * corners[2];
                            point2 = pos + ego_R * corners[3];
                            new_point1 = pos + ego_R * (corners[2] + Eigen::Vector2d(0, -step));
                            new_point2 = pos + ego_R * (corners[3] + Eigen::Vector2d(0, -step));
                            break;
                        case 3: // left
                            point1 = pos + ego_R * corners[3];
                            point2 = pos + ego_R * corners[0];
                            new_point1 = pos + ego_R * (corners[3] + Eigen::Vector2d(-step, 0));
                            new_point2 = pos + ego_R * (corners[0] + Eigen::Vector2d(-step, 0));
                            break;
                        }

                        isCollision = map_ptr->checkLineCollision(point1, new_point1) ||
                                      map_ptr->checkLineCollision(new_point1, new_point2) ||
                                      map_ptr->checkLineCollision(new_point2, point2);

                        if (isCollision)
                        {
                            is_expanding[i] = 0.0;
                        }
                        else
                        {
                            switch (i)
                            {
                            case 0: // top
                                corners[0][1] += step;
                                corners[1][1] += step;
                                if (abs(corners[0][1] - maxY) > limitBound)
                                    is_expanding[i] = 0.0;
                                break;
                            case 1: // right
                                corners[1][0] += step;
                                corners[2][0] += step;
                                if (abs(corners[1][0] - maxX) > limitBound)
                                    is_expanding[i] = 0.0;
                                break;
                            case 2: // bottom
                                corners[2][1] -= step;
                                corners[3][1] -= step;
                                if (abs(corners[2][1] - minY) > limitBound)
                                    is_expanding[i] = 0.0;
                                break;
                            case 3: // left
                                corners[3][0] -= step;
                                corners[0][0] -= step;
                                if (abs(corners[3][0] - minX) > limitBound)
                                    is_expanding[i] = 0.0;
                                break;
                            }
                        }
                    }
                }

                Eigen::Vector2d new_point, norm;
                Eigen::MatrixXd hPoly(4, 4);

                new_point = pos + ego_R * corners[0];
                norm << -cos(yaw), -sin(yaw);
                hPoly.col(0).head<2>() = norm;
                hPoly.col(0).tail<2>() = new_point;

                new_point = pos + ego_R * corners[1];
                norm << -sin(yaw), cos(yaw);
                hPoly.col(1).head<2>() = norm;
                hPoly.col(1).tail<2>() = new_point;

                new_point = pos + ego_R * corners[2];
                norm << cos(yaw), sin(yaw);
                hPoly.col(2).head<2>() = norm;
                hPoly.col(2).tail<2>() = new_point;

                new_point = pos + ego_R * corners[3];
                norm << sin(yaw), -cos(yaw);
                hPoly.col(3).head<2>() = norm;
                hPoly.col(3).tail<2>() = new_point;

                hPolys_list[k].push_back(hPoly);
            }
        }
    }

} // namespace corridor_utils