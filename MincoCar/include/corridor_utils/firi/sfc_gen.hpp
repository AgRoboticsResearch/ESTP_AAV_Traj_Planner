/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef _SFC_GEN_HPP_
#define _SFC_GEN_HPP_

#include "firi.hpp"
#include "geo_utils/geo_utils_3d.hpp"
#include "map_utils/voxel_map.hpp"

#include <deque>
#include <memory>
#include <Eigen/Eigen>
#include <cmath>
#include <vector>
#include <utility>

namespace sfc_gen
{
    inline void convexCover(const std::vector<Eigen::Vector3d> &path,
                            const std::vector<std::vector<Eigen::Vector3d>> &vertices_list,
                            const std::vector<Eigen::Vector3d> &points,
                            const Eigen::Vector3d &lowCorner,
                            const Eigen::Vector3d &highCorner,
                            const double &range,
                            std::vector<Eigen::MatrixX4d> &hpolys,
                            const double eps = 1.0e-6)
    {
        hpolys.clear();
        if (path.size() != vertices_list.size())
        {
            return;
        }

        const int n = path.size();
        Eigen::Matrix<double, 6, 4> bd = Eigen::Matrix<double, 6, 4>::Zero();
        bd(0, 0) = 1.0;
        bd(1, 0) = -1.0;
        bd(2, 1) = 1.0;
        bd(3, 1) = -1.0;
        bd(4, 2) = 1.0;
        bd(5, 2) = -1.0;

        Eigen::MatrixX4d hp, gap;
        Eigen::Vector3d a;
        std::vector<Eigen::Vector3d> valid_pc;
        valid_pc.reserve(points.size());
        for (int i = 0; i < n; i++)
        {
            a = path[i];

            bd(0, 3) = -std::min(a(0) + range, highCorner(0));
            bd(1, 3) = +std::max(a(0) - range, lowCorner(0));
            bd(2, 3) = -std::min(a(1) + range, highCorner(1));
            bd(3, 3) = +std::max(a(1) - range, lowCorner(1));
            // bd(4, 3) = -std::min(a(2) + range, highCorner(2));
            // bd(5, 3) = +std::max(a(2) - range, lowCorner(2));
            bd(4, 3) = -highCorner(2);
            bd(5, 3) = lowCorner(2);

            valid_pc.clear();
            for (const Eigen::Vector3d &p : points)
            {
                if ((bd.leftCols<3>() * p + bd.rightCols<1>()).maxCoeff() < 0.0)
                {
                    valid_pc.emplace_back(p);
                }
            }
            Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>> pc(valid_pc[0].data(), 3, valid_pc.size());
            Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>> rv(vertices_list[i][0].data(), 3, vertices_list[i].size());

            firi::firi(bd, pc, rv, a, hp);

            if (hpolys.size() != 0)
            {
                const Eigen::Vector4d ah(a(0), a(1), a(2), 1.0);
                if (3 <= ((hp * ah).array() > -eps).cast<int>().sum() +
                             ((hpolys.back() * ah).array() > -eps).cast<int>().sum())
                {
                    firi::firi(bd, pc, rv, a, gap, 1);
                    hpolys.emplace_back(gap);
                }
            }

            hpolys.emplace_back(hp);
        }
    }

    inline void shortCut(std::vector<Eigen::MatrixX4d> &hpolys)
    {
        std::vector<Eigen::MatrixX4d> htemp = hpolys;
        if (htemp.size() == 1)
        {
            Eigen::MatrixX4d headPoly = htemp.front();
            htemp.insert(htemp.begin(), headPoly);
        }
        hpolys.clear();

        int M = htemp.size();
        Eigen::MatrixX4d hPoly;
        bool overlap;
        std::deque<int> idices;
        idices.push_front(M - 1);
        for (int i = M - 1; i >= 0; i--)
        {
            for (int j = 0; j < i; j++)
            {
                if (j < i - 1)
                {
                    overlap = geo_utils_3d::overlap(htemp[i], htemp[j], 0.01);
                }
                else
                {
                    overlap = true;
                }
                if (overlap)
                {
                    idices.push_front(j);
                    i = j + 1;
                    break;
                }
            }
        }
        for (const auto &ele : idices)
        {
            hpolys.push_back(htemp[ele]);
        }
    }

    /**
     * Generate a sparse flight corridor
     *
     * @param hPolys A vector of polytopes in H-representation
     * @param map_ptr A voxel map pointer
     * @param startGoal A vector of start and goal positions
     * @param route A vector of waypoints
     * @param vertices_list A vector of vehicle vertices
     * @param range Range of the corridor
     */
    void generateSFC(std::vector<Eigen::MatrixX4d> &hPolys,
                     const map_utils::VoxelMap::Ptr &map_ptr,
                     const std::vector<Eigen::Vector3d> &startGoal,
                     const std::vector<Eigen::Vector3d> &route,
                     const std::vector<std::vector<Eigen::Vector3d>> &vertices_list,
                     double range = 10.0)
    {
        hPolys.clear();
        if (route.empty())
            return;

        std::vector<Eigen::Vector3d> pc;
        map_ptr->getSurf(pc);
        convexCover(route,
                    vertices_list,
                    pc,
                    map_ptr->getOrigin(),
                    map_ptr->getCorner(),
                    range,
                    hPolys);

        // if (hPolys.size() > 1)
        // {
        //     shortCut(hPolys);
        // }
    }

    /**
     *  Convert polytopes in H-representation to vertices
     *
     * @param hPolys A vector of polytopes in H-representation
     */
    std::vector<std::vector<Eigen::Vector3d>> convertHPolytoVertex(const std::vector<Eigen::MatrixX4d> &hPolys)
    {
        std::vector<std::vector<Eigen::Vector3d>> groupedVertices;
        for (const auto &hPoly : hPolys)
        {
            Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vPoly;
            geo_utils_3d::enumerateVs(hPoly, vPoly);

            std::vector<Eigen::Vector3d> vertices;
            for (int i = 0; i < vPoly.cols(); i++)
            {
                vertices.emplace_back(vPoly.col(i));
            }

            groupedVertices.emplace_back(vertices);
        }

        return groupedVertices;
    }

    /**
     *  Convert polytopes in H-representation to triangle meshs
     *
     * @param hPolys A vector of polytopes in H-representation
     */
    Eigen::Matrix3Xd convertHPolytoMesh(const std::vector<Eigen::MatrixX4d> &hPolys)
    {
        // Due to the fact that H-representation cannot be directly visualized
        // We first conduct vertex enumeration of them, then apply quickhull
        // to obtain triangle meshs of polyhedra
        Eigen::Matrix3Xd mesh(3, 0), curTris(3, 0), oldTris(3, 0);
        for (size_t id = 0; id < hPolys.size(); id++)
        {
            oldTris = mesh;
            Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vPoly;
            geo_utils_3d::enumerateVs(hPolys[id], vPoly);

            quickhull::QuickHull<double> tinyQH;
            const auto polyHull = tinyQH.getConvexHull(vPoly.data(), vPoly.cols(), false, true);
            const auto &idxBuffer = polyHull.getIndexBuffer();
            int hNum = idxBuffer.size() / 3;

            curTris.resize(3, hNum * 3);
            for (int i = 0; i < hNum * 3; i++)
            {
                curTris.col(i) = vPoly.col(idxBuffer[i]);
            }
            mesh.resize(3, oldTris.cols() + curTris.cols());
            mesh.leftCols(oldTris.cols()) = oldTris;
            mesh.rightCols(curTris.cols()) = curTris;
        }

        return mesh;
    }

} // namespace sfc_gen

#endif // _SFC_GEN_HPP_
