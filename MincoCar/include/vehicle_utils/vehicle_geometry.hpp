#ifndef VEHICLE_GEOMETRY_HPP
#define VEHICLE_GEOMETRY_HPP

#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Eigen>
#include "circle_cover/circle_cover.hpp"

namespace vehicle_utils
{
    using namespace geometry_cover;
    class VehicleGeometry
    {
    public:
        VehicleGeometry(){};
        ~VehicleGeometry(){};
        VehicleGeometry(const std::vector<Eigen::Vector2d> &vertices, double car_wheelbase)
        {
            setVehicleVertices(vertices);
            car_vertex_.clear();
            for (const auto &vertex : getCarVertices())
            {
                car_vertex_.push_back(vertex);
            }
            if (car_vertex_.front() != car_vertex_.back())
            {
                car_vertex_.push_back(car_vertex_.front());
            }
            car_wheelbase_ = car_wheelbase;
            // kino_astar_path_finder_ = std::make_shared<KinoAstar>(config_string, map_ptr, vehicle_geo_ptr);
            circleCover = std::make_shared<CircleCover>();
        };

        VehicleGeometry(const std::vector<Eigen::Vector2d> &vertices, double car_wheelbase, double max_sagitta)
            : VehicleGeometry(vertices, car_wheelbase)
        {
            setMaxSagitta(max_sagitta);
        };

        typedef std::shared_ptr<VehicleGeometry> Ptr;
        std::vector<Eigen::Vector2d> car_vertex_;

        void setVehicleVertices(const std::vector<Eigen::Vector2d> &vertices)
        {
            car_vertices_.clear();
            car_vertices_ = vertices;
            std::tuple<double, double, double, double> bbox = calculateBoundingBox(vertices);
            minX_ = std::get<0>(bbox);
            maxX_ = std::get<1>(bbox);
            minY_ = std::get<2>(bbox);
            maxY_ = std::get<3>(bbox);

            car_length_ = maxX_ - minX_;
            car_width_ = maxY_ - minY_;
            car_c2r_dist_ = (maxX_ + minX_) / 2.0;
        }

        void setAuxVertices(const std::vector<std::vector<Eigen::Vector2d>> &aux_vertices_list)
        {
            aux_vertices_list_.clear();
            aux_vertices_list_ = aux_vertices_list;
        }

        void setMaxSagitta(double max_sagitta)
        {
            max_sagitta_ = max_sagitta;

            circle_result_ = circleCover->computeCoverCircleInInitPose(car_vertices_,max_sagitta_);
            aux_vertices.clear();
        // Vector to store all rectangle groups
            std::vector<std::vector<Eigen::Vector2d>> rectangle_groups;

            // Loop to get all auxiliary vertices groups
            for (int i = 0; i < getAuxCount(); ++i) {
                // Get the i-th group of auxiliary vertices
                std::vector<Eigen::Vector2d> aux_vertices_part = getAuxVertices(i);

                // Ensure we only add groups of 5 vertices (a rectangle with a repeated vertex)
                if (aux_vertices_part.size() >= 5) {
                    // Extract groups of 5 vertices
                    for (size_t j = 0; j + 4 < aux_vertices_part.size(); j += 5) {
                        // Create a new vector for the current rectangle (excluding the repeated vertex)
                        std::vector<Eigen::Vector2d> rectangle(aux_vertices_part.begin() + j, aux_vertices_part.begin() + j + 4);
                        rectangle_groups.push_back(rectangle);
                    }
                }
            }

            // Output each group of rectangle vertices for debugging purposes
            std::cout << "Rectangle groups:" << std::endl;
            int group_index = 1;
            for (const auto &group : rectangle_groups) {
                std::cout << "Group " << group_index++ << ":" << std::endl;
                for (const auto &vertex : group) {
                    std::cout << "(" << vertex.x() << ", " << vertex.y() << ")" << std::endl;
                }
            }

            // Add each group of rectangle vertices to the cover circles
            for (const auto &rectangle : rectangle_groups) {
                if (!rectangle.empty()) {
                    circle_result_ = circleCover->addAuxVerticesCoverCircleInInitPose(circle_result_, rectangle, circle_result_.radius);
                }
            }
           
        }

        /* Helper functions */
        const double getLength() const { return car_length_; }
        const double getWidth() const { return car_width_; }
        const double getCenterToRearDist() const { return car_c2r_dist_; }
        const double getWheelbase() const { return car_wheelbase_; }

        // car
        const std::vector<Eigen::Vector2d> getCarVertices() const { return car_vertices_; }
        std::vector<Eigen::Vector2d> getCarVertices(const Eigen::Vector3d &state)
        {
            Eigen::Vector2d pos = state.head(2);
            double yaw = state(2);
            Eigen::Matrix2d ego_R;
            ego_R << cos(yaw), -sin(yaw),
                sin(yaw), cos(yaw);

            std::vector<Eigen::Vector2d> vertices;
            vertices.reserve(car_vertices_.size());
            for (const auto &vertex : car_vertices_)
            {
                Eigen::Vector2d point;
                point = ego_R * vertex + pos;
                vertices.push_back(Eigen::Vector2d(point[0], point[1]));
            }
            return vertices;
        }

        // aux
        const int getAuxCount() const { return aux_vertices_list_.size(); }

        const std::vector<std::vector<Eigen::Vector2d>> getAllAuxVertices() const { return aux_vertices_list_; }

        std::vector<std::vector<Eigen::Vector2d>> getAllAuxVertices(const Eigen::Vector3d &state)
        {
            Eigen::Vector2d pos = state.head(2);
            double yaw = state(2);
            Eigen::Matrix2d ego_R;
            ego_R << cos(yaw), -sin(yaw),
                sin(yaw), cos(yaw);

            std::vector<std::vector<Eigen::Vector2d>> vertices_list;
            for (size_t i = 0; i < aux_vertices_list_.size(); i++)
            {
                std::vector<Eigen::Vector2d> vertices;
                vertices.reserve(aux_vertices_list_[i].size());
                for (const auto &vertex : aux_vertices_list_[i])
                {
                    Eigen::Vector2d point;
                    point = ego_R * vertex + pos;
                    vertices.push_back(Eigen::Vector2d(point[0], point[1]));
                }
                vertices_list.push_back(vertices);
            }
            return vertices_list;
        }

        const std::vector<Eigen::Vector2d> getAuxVertices(int i) const
        {
            if (i < 0 || i > aux_vertices_list_.size())
            {
                std::cerr << "Invalid Auxiliary index" << std::endl;
                return std::vector<Eigen::Vector2d>();
            }

            return aux_vertices_list_[i];
        }

        std::vector<Eigen::Vector2d> getAuxVertices(int i, const Eigen::Vector3d &state)
        {
            if (i < 0 || i > aux_vertices_list_.size())
            {
                std::cerr << "Invalid auxiliary index" << std::endl;
                return std::vector<Eigen::Vector2d>();
            }

            Eigen::Vector2d pos = state.head(2);
            double yaw = state(2);
            Eigen::Matrix2d ego_R;
            ego_R << cos(yaw), -sin(yaw),
                sin(yaw), cos(yaw);

            std::vector<Eigen::Vector2d> vertices;
            vertices.reserve(aux_vertices_list_[i].size());
            for (const auto &vertex : aux_vertices_list_[i])
            {
                Eigen::Vector2d point;
                point = ego_R * vertex + pos;
                vertices.push_back(Eigen::Vector2d(point[0], point[1]));
            }
            return vertices;
        }

        // cover circle
        const CoverCircleResult getCoverCircle() const { return circle_result_; }

        // bounding box
        const double getMinX() const { return minX_; }
        const double getMaxX() const { return maxX_; }
        const double getMinY() const { return minY_; }
        const double getMaxY() const { return maxY_; }

        std::tuple<double, double, double, double> calculateBoundingBox(const std::vector<Eigen::Vector2d> &vertices)
        {
            double minX = vertices[0][0];
            double maxX = vertices[0][0];
            double minY = vertices[0][1];
            double maxY = vertices[0][1];

            for (int i = 1; i < vertices.size(); i++)
            {
                if (vertices[i][0] < minX)
                    minX = vertices[i][0];
                if (vertices[i][0] > maxX)
                    maxX = vertices[i][0];
                if (vertices[i][1] < minY)
                    minY = vertices[i][1];
                if (vertices[i][1] > maxY)
                    maxY = vertices[i][1];
            }
            return std::make_tuple(minX, maxX, minY, maxY);
        }

    private:
        double car_length_;
        double car_width_;
        double car_c2r_dist_;
        double car_wheelbase_;
        double max_sagitta_;
        std::vector<Eigen::Vector2d> car_vertices_;
        std::vector<std::vector<Eigen::Vector2d>> aux_vertices_list_;
        std::shared_ptr<CircleCover> circleCover;
        // bounding box
        double minX_, maxX_, minY_, maxY_;
        std::vector<Eigen::Vector2d> aux_vertices;
        // cover circle
        CoverCircleResult circle_result_;
    };
} // namespace vehicle_utils

#endif // VEHICLE_GEOMETRY_HPP