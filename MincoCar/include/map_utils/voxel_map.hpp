#ifndef _MAPPING_HPP_
#define _MAPPING_HPP_

#include <iostream>
#include <math.h>
#include <Eigen/Eigen>
#include <queue>
#include <memory>

#include "raycast.hpp"
#include "voxel_dilater.hpp"

namespace map_utils
{
    enum VoxelState : uint8_t {
        UNOCCUPIED = 0,
        OCCUPIED = 1,
        DILATED = 2,
        VIRTUAL_DILATED = 3
    };

    class VoxelMap
    {
    public:
        VoxelMap(){};
        VoxelMap(const Eigen::Vector3d &map_size,
                 const Eigen::Vector3d &origin,
                 double resolution);
        VoxelMap(const Eigen::Vector3i &grid_size,
                 const Eigen::Vector3d &origin,
                 double resolution);
        ~VoxelMap(){};
        typedef std::shared_ptr<VoxelMap> Ptr;

        void reset_buffer();
        void dilate(int num, bool is_virtual);

        bool setOccupancy(const Eigen::Vector3d &pos, uint8_t occ);
        bool setOccupancy(const Eigen::Vector3i &id, uint8_t occ);
        bool setOccupancy(const Eigen::Vector2d &pos, uint8_t occ);
        bool setOccupancy(const Eigen::Vector2i &id, uint8_t occ);

        void posToIndex(const Eigen::Vector3d &pos, Eigen::Vector3i &id) const;
        void posToIndex(const Eigen::Vector2d &pos, Eigen::Vector2i &id) const;
        void indexToPos(const Eigen::Vector3i &id, Eigen::Vector3d &pos) const;
        void indexToPos(const Eigen::Vector2i &id, Eigen::Vector2d &pos) const;

        // check
        bool isValidState(const uint8_t &state);

        bool isInMap(const Eigen::Vector3d &pos);
        bool isInMap(const Eigen::Vector3i &id);
        bool isInMap(const Eigen::Vector2i &pos);
        bool isInMap(const Eigen::Vector2d &id);

        bool checkPointCollision(const Eigen::Vector2d &point); 
        bool checkLineCollision(const Eigen::Vector2d &start_pt, const Eigen::Vector2d &end_pt);

        // helper functions
        double getResolution() const { return resolution_; }
        const Eigen::Vector3d getOrigin() const { return origin_; }
        const Eigen::Vector3d getMapSize() const { return map_size_; }
        const Eigen::Vector3i getGridSize() const { return grid_size_; }
        const Eigen::Vector3d getCorner() const { return map_size_ + origin_; }

        uint8_t getVoxelState(const Eigen::Vector3d &pos);
        uint8_t getVoxelState(const Eigen::Vector3i &id);
        uint8_t getVoxelState(const Eigen::Vector2d &pos);
        uint8_t getVoxelState(const Eigen::Vector2i &id);
        std::vector<uint8_t> getVoxels() const { return occupancy_buffer_; }

        template <typename VectorType>
        uint8_t query(const VectorType &vector) { return getVoxelState(vector); };

        const std::vector<Eigen::Vector3i> getSurfIdx() const { return surf_; };
        void getSurf(std::vector<Eigen::Vector3d> &points) const;

    private:
        Eigen::Vector3d origin_, map_size_;
        Eigen::Vector3i grid_size_;
        double resolution_;
        int grid_size_x_multiply_y_;
        std::vector<Eigen::Vector3i> surf_;
        bool has_dilated_;

        std::vector<uint8_t> occupancy_buffer_; // This buffer stores the states of each grid in the "global" map
    };

} // namespace map_utils

#endif // _MAPPING_HPP_