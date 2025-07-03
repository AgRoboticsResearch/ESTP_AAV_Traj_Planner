#ifndef CIRCLE_COVER_HPP
#define CIRCLE_COVER_HPP

#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigen>



namespace geometry_cover
{

    struct CoverCircleResult
    {
        double radius;
        double aux_radius;
        std::vector<Eigen::Vector2d> centers;
        std::vector<Eigen::Vector2d> aux_centers;
    };
    class CircleCover {
    public:
        CircleCover();
        ~CircleCover();

    typedef std::shared_ptr<CircleCover> Ptr;

    CoverCircleResult computeCoverCircleInInitPose(std::vector<Eigen::Vector2d> car_vertex,double max_sagitta);
    CoverCircleResult addAuxVerticesCoverCircleInInitPose(CoverCircleResult cover_circle_result_,std::vector<Eigen::Vector2d> aux_vertices,double max_radius);
    private:
    
    };

} // namespace circle_cover

#endif // CIRCLE_COVER_HPP
