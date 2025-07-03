#include "map_utils/voxel_map.hpp"

namespace map_utils
{
    VoxelMap::VoxelMap(const Eigen::Vector3d &map_size,
                       const Eigen::Vector3d &origin,
                       double resolution)
    {
        origin_ = origin;
        resolution_ = resolution;
        for (int i = 0; i < 3; i++)
        {
            double size = (map_size(i) <= 0) ? resolution_ : map_size(i);
            grid_size_(i) = ceil(size / resolution_);
        }

        map_size_ = grid_size_.cast<double>() * resolution_;
        reset_buffer();
    }

    VoxelMap::VoxelMap(const Eigen::Vector3i &grid_size,
                       const Eigen::Vector3d &origin,
                       double resolution)
    {
        grid_size_ = grid_size;
        origin_ = origin;
        resolution_ = resolution;
        for (int i = 0; i < 3; i++)
        {
            map_size_(i) = grid_size_(i) * resolution_;
        }

        map_size_ = grid_size_.cast<double>() * resolution_;
        reset_buffer();
    }

    void VoxelMap::reset_buffer()
    {
        // initialize size of buffer
        grid_size_x_multiply_y_ = grid_size_(0) * grid_size_(1);
        int buffer_size_ = grid_size_(2) * grid_size_x_multiply_y_;
        occupancy_buffer_.resize(buffer_size_);
        fill(occupancy_buffer_.begin(), occupancy_buffer_.end(), UNOCCUPIED);

        has_dilated_ = false;
    }

    void VoxelMap::dilate(int num, bool is_virtual = false)
    {
        if (num <= 0)
        {
            return;
        }
        else
        {
            int check_value;
            if (has_dilated_ && is_virtual)
            {
                check_value = DILATED;
            }
            else
            {
                check_value = OCCUPIED;
            }
            int dilate_value = is_virtual ? VIRTUAL_DILATED : DILATED;

            std::vector<Eigen::Vector3i> lvec, cvec;
            lvec.reserve(grid_size_.prod());
            cvec.reserve(grid_size_.prod());
            int i, j, k, idx;
            bool check;
            Eigen::Vector3d bounds;
            bounds << grid_size_(0) - 1,
                (grid_size_(1) - 1) * grid_size_(0),
                (grid_size_(2) - 1) * grid_size_x_multiply_y_;
            for (int x = 0; x <= bounds(0); x++)
            {
                for (int y = 0; y <= bounds(1); y += grid_size_(0))
                {
                    for (int z = 0; z <= bounds(2); z += grid_size_x_multiply_y_)
                    {

                        if (occupancy_buffer_[x + y + z] == check_value)
                        {
                            VOXEL_DILATER(i, j, k,
                                          x, y, z,
                                          grid_size_(0), grid_size_x_multiply_y_,
                                          bounds(0), bounds(1), bounds(2),
                                          check, occupancy_buffer_, idx, dilate_value, cvec)
                        }
                    }
                }
            }

            for (int loop = 1; loop < num; loop++)
            {
                std::swap(cvec, lvec);
                for (const Eigen::Vector3i &id : lvec)
                {
                    VOXEL_DILATER(i, j, k,
                                  id(0), id(1), id(2),
                                  grid_size_(0), grid_size_x_multiply_y_,
                                  bounds(0), bounds(1), bounds(2),
                                  check, occupancy_buffer_, idx, dilate_value, cvec)
                }
                lvec.clear();
            }

            if (!is_virtual)
            {
                has_dilated_ = true;
                surf_ = cvec;
            }            
        }
    }

    bool VoxelMap::setOccupancy(const Eigen::Vector3d &pos, uint8_t occ = OCCUPIED)
    {
        if (!isValidState(occ))
        {
            return false;
        }

        Eigen::Vector3i id;
        posToIndex(pos, id);

        if (!isInMap(id))
        {
            return false;
        }

        int idx_ctns = id(0) + id(1) * grid_size_(0) + id(2) * grid_size_x_multiply_y_;
        occupancy_buffer_[idx_ctns] = occ;

        return true;
    }

    bool VoxelMap::setOccupancy(const Eigen::Vector3i &id, uint8_t occ = OCCUPIED)
    {
        if (!isValidState(occ))
        {
            return false;
        }

        if (!isInMap(id))
        {
            return false;
        }

        int idx_ctns = id(0) + id(1) * grid_size_(0) + id(2) * grid_size_x_multiply_y_;
        occupancy_buffer_[idx_ctns] = occ;

        return true;
    }

    bool VoxelMap::setOccupancy(const Eigen::Vector2d &pos, uint8_t occ = OCCUPIED)
    {
        if (!isValidState(occ))
        {
            return false;
        }

        Eigen::Vector2i id;
        posToIndex(pos, id);

        if (!isInMap(id))
        {
            return false;
        }

        for (int z_id = 0; z_id < grid_size_(2); z_id++)
        {
            int idx_ctns = id(0) + id(1) * grid_size_(0) + z_id * grid_size_x_multiply_y_;
            occupancy_buffer_[idx_ctns] = occ;
        }

        return true;
    }

    bool VoxelMap::setOccupancy(const Eigen::Vector2i &id, uint8_t occ = OCCUPIED)
    {
        if (!isValidState(occ))
        {
            return false;
        }

        if (!isInMap(id))
        {
            return false;
        }

        for (int z_id = 0; z_id < grid_size_(2); z_id++)
        {
            int idx_ctns = id(0) + id(1) * grid_size_(0) + z_id * grid_size_x_multiply_y_;
            occupancy_buffer_[idx_ctns] = occ;
        }

        return true;
    }

    void VoxelMap::posToIndex(const Eigen::Vector3d &pos, Eigen::Vector3i &id) const
    {
        for (int i = 0; i < 3; i++)
            id(i) = floor((pos(i) - origin_(i)) / resolution_);
    }

    void VoxelMap::posToIndex(const Eigen::Vector2d &pos, Eigen::Vector2i &id) const
    {
        for (int i = 0; i < 2; i++)
            id(i) = floor((pos(i) - origin_(i)) / resolution_);
    }

    void VoxelMap::indexToPos(const Eigen::Vector3i &id, Eigen::Vector3d &pos) const
    {
        for (int i = 0; i < 3; i++)
            pos(i) = (id(i) + 0.5) * resolution_ + origin_(i);
    }

    void VoxelMap::indexToPos(const Eigen::Vector2i &id, Eigen::Vector2d &pos) const
    {
        for (int i = 0; i < 2; i++)
            pos(i) = (id(i) + 0.5) * resolution_ + origin_(i);
    }

    bool VoxelMap::isValidState(const uint8_t &state)
    {
        if (state != UNOCCUPIED && state != OCCUPIED && state != DILATED && state != VIRTUAL_DILATED)
        {
            return false;
        }
        return true;
    }

    bool VoxelMap::isInMap(const Eigen::Vector3d &pos)
    {
        Eigen::Vector3i idx;
        posToIndex(pos, idx);
        return isInMap(idx);
    }

    bool VoxelMap::isInMap(const Eigen::Vector3i &id)
    {
        if (id(0) < 0 || id(0) > grid_size_(0) - 1 || id(1) < 0 || id(1) > grid_size_(1) - 1 || id(2) < 0 || id(2) > grid_size_(2) - 1)
        {
            return false;
        }
        else
        {
            return true;
        }
    };

    bool VoxelMap::isInMap(const Eigen::Vector2d &pos)
    {
        Eigen::Vector2i idx;
        posToIndex(pos, idx);
        return isInMap(idx);
    }

    bool VoxelMap::isInMap(const Eigen::Vector2i &id)
    {
        if (id(0) < 0 || id(0) > grid_size_(0) - 1 || id(1) < 0 || id(1) > grid_size_(1) - 1)
        {
            return false;
        }
        else
            return true;
    };

    bool VoxelMap::checkPointCollision(const Eigen::Vector2d &point)
    {
        // Check if point is within map boundaries
        if (!isInMap(point))
        {
            return true; // Collision detected if point is outside map
        }

        // Check collision with obstacles in the map
        if (getVoxelState(point) > UNOCCUPIED && getVoxelState(point) != VIRTUAL_DILATED)
        {
            return true; // Collision detected if point is inside an obstacle
        }
        else
        {
            return false; // No collision detected
        }
    }

    bool VoxelMap::checkLineCollision(const Eigen::Vector2d &start_pt, const Eigen::Vector2d &end_pt)
    {
        if (checkPointCollision(start_pt) || checkPointCollision(end_pt))
        {
            return true;
        }

        RayCaster raycaster;
        bool need_ray = raycaster.setInput(start_pt / resolution_, end_pt / resolution_);
        // if(!need_ray)
        // {
        //     return false;
        // }
        Eigen::Vector2d half = Eigen::Vector2d(0.5, 0.5);
        Eigen::Vector2d ray_pt;
        // if(!raycaster.step(ray_pt))
        // {
        //     return false;
        // }
        while (raycaster.step(ray_pt))
        {
            Eigen::Vector2d tmp = (ray_pt + half) * resolution_;
            if (checkPointCollision(tmp))
            {
                return true;
            }
        }
        return false;
    }

    uint8_t VoxelMap::getVoxelState(const Eigen::Vector3d &pos)
    {
        Eigen::Vector3i id;
        posToIndex(pos, id);
        if (!isInMap(id))
            return -1;

        return occupancy_buffer_[id(0) + id(1) * grid_size_(0) + id(2) * grid_size_x_multiply_y_];
    }

    uint8_t VoxelMap::getVoxelState(const Eigen::Vector3i &id)
    {
        if (!isInMap(id))
            return -1;

        return occupancy_buffer_[id(0) + id(1) * grid_size_(0) + id(2) * grid_size_x_multiply_y_];
    }

    uint8_t VoxelMap::getVoxelState(const Eigen::Vector2d &pos)
    {
        Eigen::Vector2i id;
        posToIndex(pos, id);
        if (!isInMap(id))
            return -1;

        return occupancy_buffer_[id(0) + id(1) * grid_size_(0)];
    }

    uint8_t VoxelMap::getVoxelState(const Eigen::Vector2i &id)
    {
        if (!isInMap(id))
            return -1;

        return occupancy_buffer_[id(0) + id(1) * grid_size_(0)];
    }

    void VoxelMap::getSurf(std::vector<Eigen::Vector3d> &points) const
    {
        points.clear();
        points.reserve(surf_.size());

        Eigen::Vector3i step(1, grid_size_(0), grid_size_x_multiply_y_);
        Eigen::Vector3d oc(origin_ + Eigen::Vector3d::Constant(0.5 * resolution_));
        Eigen::Vector3d stepScale(step.cast<double>().cwiseInverse() * resolution_);

        for (const Eigen::Vector3i &id : surf_)
        {
            points.push_back(id.cast<double>().cwiseProduct(stepScale) + oc);
        }
    }

} // namespace map_utils