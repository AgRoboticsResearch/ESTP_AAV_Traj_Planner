#include "corridor_utils/safe_corridor.hpp"
#include "corridor_utils/firi/sfc_gen.hpp"

namespace corridor_utils
{

    // Calculate the centroid of the polygon
    Eigen::Vector2d calculateCentroid(const std::vector<Eigen::Vector3d> &vertices)
    {
        Eigen::Vector2d centroid(0, 0);
        int n = vertices.size();
        for (const Eigen::Vector3d &p : vertices)
        {
            centroid += Eigen::Vector2d(p.x(), p.y());
        }
        centroid /= n;
        return centroid;
    }

    // Calculate the angle between two points with respect to the origin
    double calculateAngle(const Eigen::Vector2d &origin, const Eigen::Vector2d &p)
    {
        return atan2(p.y() - origin.y(), p.x() - origin.x());
    }

    void calcFIRICorridor(std::vector<Eigen::MatrixXd> &hPolys,
                          const std::vector<Eigen::Vector3d> &statelist,
                          const map_utils::VoxelMap::Ptr &map_ptr,
                          const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                          double range)
    {
        hPolys.clear();
        double resolution = map_ptr->getResolution();

        std::vector<Eigen::Vector3d> startGoal;
        Eigen::Vector3d start = statelist.front();
        Eigen::Vector3d goal = statelist.back();
        start[2] = resolution / 2.0;
        goal[2] = resolution / 2.0;
        startGoal.push_back(start);
        startGoal.push_back(goal);

        std::vector<Eigen::Vector3d> route = statelist;
        for (auto &point : route)
        {
            point[2] = resolution / 2.0;
        }

        // get locations of vehicle vertices
        std::vector<std::vector<Eigen::Vector3d>> vertices_list;
        vertices_list.reserve(statelist.size());

        for (const auto &state : statelist)
        {
            std::vector<Eigen::Vector2d> vertices_2d = vehicle_geo_ptr->getCarVertices(state);
            std::vector<Eigen::Vector3d> vertices;
            vertices.reserve(vertices_2d.size());
            for (const auto &vertex : vertices_2d)
            {
                vertices.push_back(Eigen::Vector3d(vertex[0], vertex[1], resolution / 2.0));
            }

            vertices_list.push_back(vertices);
        }

        // Each row of hPolys_firi is defined by h0, h1, h2, h3 as
        // h0*x + h1*y + h2*z + h3 <= 0
        std::vector<Eigen::MatrixX4d> hPolys_firi;
        sfc_gen::generateSFC(hPolys_firi, map_ptr, startGoal, route, vertices_list, range);

        std::vector<std::vector<Eigen::Vector3d>> vPolys;
        vPolys = sfc_gen::convertHPolytoVertex(hPolys_firi);

        for (auto &vPoly : vPolys)
        {
            vPoly.erase(std::remove_if(
                            vPoly.begin(), vPoly.end(), [](const Eigen::Vector3d &vertex)
                            { return vertex.z() >= 1.0e-6; }),
                        vPoly.end());

            // reorder vertices clockwisely
            Eigen::Vector2d centroid = calculateCentroid(vPoly);
            std::sort(vPoly.begin(), vPoly.end(),
                      [&centroid](const Eigen::Vector3d &a, const Eigen::Vector3d &b)
                      {
                          Eigen::Vector2d pointA(a.x(), a.y());
                          Eigen::Vector2d pointB(b.x(), b.y());
                          double angleA = calculateAngle(centroid, pointA);
                          double angleB = calculateAngle(centroid, pointB);
                          return angleA > angleB;
                      });

            // for (auto &vertex : vPoly)
            // {
            //     vertex.z() = 0;
            // }
        }

        // Each col of hPoly denotes a facet (outter_normal^T,point^T)^T
        // The outter_normal is assumed to be NORMALIZED
        for (const auto &vPoly : vPolys)
        {
            Eigen::MatrixXd hPoly(4, vPoly.size());

            for (int i = 0; i < vPoly.size(); ++i)
            {
                // Calculate normal
                Eigen::Vector2d current = vPoly[i].head<2>();
                Eigen::Vector2d previous = vPoly[(i + vPoly.size() - 1) % vPoly.size()].head<2>();
                Eigen::Vector2d direction = current - previous;
                Eigen::Vector2d normal(-direction[1], direction[0]);
                normal.normalize();

                hPoly(0, i) = normal[0];
                hPoly(1, i) = normal[1];
                hPoly(2, i) = current[0];
                hPoly(3, i) = current[1];
            }

            hPolys.push_back(hPoly);
        }
    }

    void calcAuxFIRICorridor(std::vector<std::vector<Eigen::MatrixXd>> &hPolys_list,
                            const std::vector<Eigen::Vector3d> &state_list,
                            const map_utils::VoxelMap::Ptr &map_ptr,
                            const vehicle_utils::VehicleGeometry::Ptr &vehicle_geo_ptr,
                            double range)
    {
        hPolys_list.clear();
        hPolys_list.resize(vehicle_geo_ptr->getAuxCount());
        if (vehicle_geo_ptr->getAuxCount() == 0)
            return;
        
        double resolution = map_ptr->getResolution();

        std::vector<Eigen::Vector3d> startGoal;
        Eigen::Vector3d start = state_list.front();
        Eigen::Vector3d goal = state_list.back();
        start[2] = resolution / 2.0;
        goal[2] = resolution / 2.0;
        startGoal.push_back(start);
        startGoal.push_back(goal);

        std::vector<Eigen::Vector3d> route = state_list;
        for (auto &point : route)
        {
            point[2] = resolution / 2.0;
        }

        for (size_t k = 0; k < vehicle_geo_ptr->getAuxCount(); k++)
        {
            // get locations of vehicle vertices
            std::vector<std::vector<Eigen::Vector3d>> vertices_list;
            vertices_list.reserve(state_list.size());

            for (const auto &state : state_list)
            {
                std::vector<Eigen::Vector2d> vertices_2d = vehicle_geo_ptr->getAuxVertices(k, state);
                std::vector<Eigen::Vector3d> vertices;
                vertices.reserve(vertices_2d.size());
                for (const auto &vertex : vertices_2d)
                {
                    vertices.push_back(Eigen::Vector3d(vertex[0], vertex[1], resolution / 2.0));
                }

                vertices_list.push_back(vertices);
            }

            // Each row of hPolys_firi is defined by h0, h1, h2, h3 as
            // h0*x + h1*y + h2*z + h3 <= 0
            std::vector<Eigen::MatrixX4d> hPolys_firi;
            sfc_gen::generateSFC(hPolys_firi, map_ptr, startGoal, route, vertices_list, range);

            std::vector<std::vector<Eigen::Vector3d>> vPolys;
            vPolys = sfc_gen::convertHPolytoVertex(hPolys_firi);

            for (auto &vPoly : vPolys)
            {
                vPoly.erase(std::remove_if(
                                vPoly.begin(), vPoly.end(), [](const Eigen::Vector3d &vertex)
                                { return vertex.z() >= 1.0e-6; }),
                            vPoly.end());

                // reorder vertices clockwisely
                Eigen::Vector2d centroid = calculateCentroid(vPoly);
                std::sort(vPoly.begin(), vPoly.end(),
                          [&centroid](const Eigen::Vector3d &a, const Eigen::Vector3d &b)
                          {
                              Eigen::Vector2d pointA(a.x(), a.y());
                              Eigen::Vector2d pointB(b.x(), b.y());
                              double angleA = calculateAngle(centroid, pointA);
                              double angleB = calculateAngle(centroid, pointB);
                              return angleA > angleB;
                          });

                // for (auto &vertex : vPoly)
                // {
                //     vertex.z() = 0;
                // }
            }

            // Each col of hPoly denotes a facet (outter_normal^T,point^T)^T
            // The outter_normal is assumed to be NORMALIZED
            for (const auto &vPoly : vPolys)
            {
                Eigen::MatrixXd hPoly(4, vPoly.size());

                for (int i = 0; i < vPoly.size(); ++i)
                {
                    // Calculate normal
                    Eigen::Vector2d current = vPoly[i].head<2>();
                    Eigen::Vector2d previous = vPoly[(i + vPoly.size() - 1) % vPoly.size()].head<2>();
                    Eigen::Vector2d direction = current - previous;
                    Eigen::Vector2d normal(-direction[1], direction[0]);
                    normal.normalize();

                    hPoly(0, i) = normal[0];
                    hPoly(1, i) = normal[1];
                    hPoly(2, i) = current[0];
                    hPoly(3, i) = current[1];
                }

                hPolys_list[k].push_back(hPoly);
            }
        }
    }

} // namespace corridor_utils