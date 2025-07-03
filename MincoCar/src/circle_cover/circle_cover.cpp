#include "circle_cover/circle_cover.hpp"
#include <iostream>
namespace geometry_cover
{
    CircleCover::CircleCover() {
        // Constructor implementation
    }
    CircleCover::~CircleCover() {
        // Destructor implementation
    }

CoverCircleResult CircleCover::computeCoverCircleInInitPose(std::vector<Eigen::Vector2d> car_vertex, double max_sagitta) {
    CoverCircleResult result;
    
    // Initialize min and max values for x and y coordinates
    double min_x = car_vertex[0].x();
    double max_x = car_vertex[0].x();
    double min_y = car_vertex[0].y();
    double max_y = car_vertex[0].y();

    // Calculate min and max values for x and y coordinates
    for (int i = 1; i < 4; ++i) {
        double x = car_vertex[i].x();
        double y = car_vertex[i].y();
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }

    double width = max_x - min_x;
    double height = max_y - min_y;

    // Ensure short_side is the smaller dimension
    double short_side = std::min(width, height);
    double long_side = std::max(width, height);

    printf("[Vehicle Geometry]: circle max_sagitta: %f \n", max_sagitta);

    double final_radius = 0.0;
    int i = 1;
    while (true) {
        int pieces_short = std::pow(2, std::max(0, i - 1));
        int pieces_long = std::pow(2, i);
        double divided_short = short_side / pieces_short;
        double divided_long = long_side / pieces_long;
        final_radius = std::sqrt((divided_short / 2.0) * (divided_short / 2.0) + (divided_long / 2.0) * (divided_long / 2.0));

        if (max_sagitta < 0) {
            printf("\033[31m[Vehicle Geometry]: circle compute error, max_sagitta<0!!!!!!!!!!!!!\n");
            max_sagitta = 0.5;
            printf("\033[31m[Vehicle Geometry]: circle compute error, max_sagitta=0.5!!!!!!!!!!!!!\n");
        }
        printf("\033[0m");

        if (max_sagitta > final_radius - divided_short / 2) {
            final_radius = std::sqrt((divided_short / 2.0) * (divided_short / 2.0) + (divided_long / 2.0) * (divided_long / 2.0));

            for (int x = 0; x < pieces_long; x++) {
                if (width == long_side) {
                    double center_x = min_x + x * divided_long + divided_long / 2.0;
                    result.centers.push_back(Eigen::Vector2d(center_x, min_y + divided_short / 2.0));
                    result.centers.push_back(Eigen::Vector2d(center_x, max_y - divided_short / 2.0));
                } else {
                    double center_y = min_y + x * divided_long + divided_long / 2.0;
                    result.centers.push_back(Eigen::Vector2d(min_x + divided_short / 2.0, center_y));
                    result.centers.push_back(Eigen::Vector2d(max_x - divided_short / 2.0, center_y));
                }
            }

            for (int y = 1; y < pieces_short - 1; y++) {
                if (width == long_side) {
                    double center_y = min_y + y * divided_short + divided_short / 2.0;
                    result.centers.push_back(Eigen::Vector2d(min_x + divided_long / 2.0, center_y));
                    result.centers.push_back(Eigen::Vector2d(max_x - divided_long / 2.0, center_y));
                } else {
                    double center_x = min_x + y * divided_short + divided_short / 2.0;
                    result.centers.push_back(Eigen::Vector2d(center_x, min_y + divided_long / 2.0));
                    result.centers.push_back(Eigen::Vector2d(center_x, max_y - divided_long / 2.0));
                }
            }
            printf("[Vehicle Geometry]: circle final level: %d \n", i);
            break;
        }
        ++i;
        if (i > 5) {
            for (int x = 0; x < pieces_long; x++) {
                if (width == long_side) {
                    double center_x = min_x + x * divided_long + divided_long / 2.0;
                    result.centers.push_back(Eigen::Vector2d(center_x, min_y + divided_short / 2.0));
                    result.centers.push_back(Eigen::Vector2d(center_x, max_y - divided_short / 2.0));
                } else {
                    double center_y = min_y + x * divided_long + divided_long / 2.0;
                    result.centers.push_back(Eigen::Vector2d(min_x + divided_short / 2.0, center_y));
                    result.centers.push_back(Eigen::Vector2d(max_x - divided_short / 2.0, center_y));
                }
            }

            for (int y = 1; y < pieces_short - 1; y++) {
                if (width == long_side) {
                    double center_y = min_y + y * divided_short + divided_short / 2.0;
                    result.centers.push_back(Eigen::Vector2d(min_x + divided_long / 2.0, center_y));
                    result.centers.push_back(Eigen::Vector2d(max_x - divided_long / 2.0, center_y));
                } else {
                    double center_x = min_x + y * divided_short + divided_short / 2.0;
                    result.centers.push_back(Eigen::Vector2d(center_x, min_y + divided_long / 2.0));
                    result.centers.push_back(Eigen::Vector2d(center_x, max_y - divided_long / 2.0));
                }
            }
            printf("[Vehicle Geometry]: circle final level: %d \n", i);
            break;
        }
    }
    result.radius = final_radius;
    printf("[Vehicle Geometry]: circle final radius: %f m\n", final_radius);
    return result;
}



//    CoverCircleResult CircleCover::addAuxVerticesUseSixSideCoverCircleInInitPose(CircleCover::CoverCircleResult cover_circle_result_,std::vector<Eigen::Vector2d> aux_vertices,double max_radius)
//   {
//        if (aux_vertices.size() < 4) {
//         std::cerr << "Error: Exactly 4 vertices are required." << std::endl;
//         return cover_circle_result_;
//         // Handle the error appropriately, maybe throw an exception or return an error state
//     }

//     // Initialize min and max values
//     double min_x = aux_vertices[0].x();
//     double max_x = aux_vertices[0].x();
//     double min_y = aux_vertices[0].y();
//     double max_y = aux_vertices[0].y();

//     // Find the bounding box of the vertices
//     for (const auto& vertex : aux_vertices) {
//         double x = vertex.x();
//         double y = vertex.y();
//         if (x < min_x) min_x = x;
//         if (x > max_x) max_x = x;
//         if (y < min_y) min_y = y;
//         if (y > max_y) max_y = y;
//     }

//     // Calculate dimensions
//     double width = max_x - min_x;
//     double height = max_y - min_y;

//     int numCirclesX = static_cast<int>(std::ceil(width / (cover_circle_result_.radius)));
//     int numCirclesY = static_cast<int>(std::ceil(height / (std::sqrt(3) * cover_circle_result_.radius)));

//     // Vector to hold all circle centers
//     std::vector<Eigen::Matrix<double, 2, 1>> all_centers_eigen;

//     // Generate circle centers
//     for (int row = 0; row < numCirclesY; ++row) {
//         for (int col = 0; col < numCirclesX; ++col) {
//             // Calculate the offset for hexagonal packing
//             double offsetX = cover_circle_result_.radius * 0.5;
//             double centerX = min_x + col * 1 * cover_circle_result_.radius + offsetX;
//             double offsetY = std::sqrt(3) * cover_circle_result_.radius * 0.5;
//             double centerY = min_y + row * std::sqrt(3) * cover_circle_result_.radius + offsetY;

//             // Add all centers to the vector
//             Eigen::Matrix<double, 2, 1> center;
//             center << centerX, centerY;
//             all_centers_eigen.push_back(center);
//         }
//     }
//   }
CoverCircleResult CircleCover::addAuxVerticesCoverCircleInInitPose(CoverCircleResult geo_cover_circle_result_, std::vector<Eigen::Vector2d> aux_vertices, double max_radius) {
  
    // Initialize min and max values for x and y coordinates
    double min_x = aux_vertices[0].x();
    double max_x = aux_vertices[0].x();
    double min_y = aux_vertices[0].y();
    double max_y = aux_vertices[0].y();

    // Calculate min and max values for x and y coordinates
    for (int i = 1; i < aux_vertices.size(); ++i) {
        double x = aux_vertices[i].x();
        double y = aux_vertices[i].y();
        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }

    double width = max_x - min_x;
    double height = max_y - min_y;

    // Ensure short_side is the smaller dimension
    double short_side = std::min(width, height);
    double long_side = std::max(width, height);

    double final_radius = std::sqrt((short_side / 2.0) * (short_side / 2.0) + (long_side / 2.0) * (long_side / 2.0));
    
    // initialize the aux circle
    // geo_cover_circle_result_.aux_centers.clear();
    geo_cover_circle_result_.aux_radius = final_radius;

    if (max_radius >= final_radius) {
        geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d((min_x + max_x) / 2.0, (min_y + max_y) / 2.0));
        printf("[Vehicle Geometry]: aux circle final level: %d\n", 0);
        printf("[Vehicle Geometry]: aux circle final radius: %f m\n", final_radius);
        return geo_cover_circle_result_;
    }

    int i = 1;
    while (true) {
        int pieces_short = std::pow(2, std::max(0, i - 1));
        int pieces_long = std::pow(2, i);
        double divided_short = short_side / pieces_short;
        double divided_long = long_side / pieces_long;
        final_radius = std::sqrt((divided_short / 2.0) * (divided_short / 2.0) + (divided_long / 2.0) * (divided_long / 2.0));

        if (max_radius >= final_radius) {
            // Add centers along the long side
            for (int x = 0; x < pieces_long; ++x) {
                if (width == long_side) {
                    double center_x = min_x + x * divided_long + divided_long / 2.0;
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(center_x, min_y + divided_short / 2.0));
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(center_x, max_y - divided_short / 2.0));
                } else {
                    double center_y = min_y + x * divided_long + divided_long / 2.0;
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(min_x + divided_short / 2.0, center_y));
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(max_x - divided_short / 2.0, center_y));
                }
            }

            // Add centers along the short side, excluding the corners already covered
            for (int y = 1; y < pieces_short - 1; ++y) {
                if (width == long_side) {
                    double center_y = min_y + y * divided_short + divided_short / 2.0;
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(min_x + divided_long / 2.0, center_y));
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(max_x - divided_long / 2.0, center_y));
                } else {
                    double center_x = min_x + y * divided_short + divided_short / 2.0;
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(center_x, min_y + divided_long / 2.0));
                    geo_cover_circle_result_.aux_centers.push_back(Eigen::Vector2d(center_x, max_y - divided_long / 2.0));
                }
            }

            printf("[Vehicle Geometry]: aux circle final level: %d\n", i);
            printf("[Vehicle Geometry]: aux circle final radius: %f m\n", final_radius);

            break;
        }
        ++i;
    }
    geo_cover_circle_result_.aux_radius = final_radius;
    return geo_cover_circle_result_;
}

} // namespace circle_cover
