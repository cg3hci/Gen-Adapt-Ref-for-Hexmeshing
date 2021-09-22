/*****************************************************************************************
 *              MIT License                                                              *
 *                                                                                       *
 * Copyright (c) 2021 Luca Pitzalis, Marco Livesu, Gianmarco Cherchi, Enrico Gobbetti    *
 * and Riccardo Scateni                                                                  *  
 *                                                                                       *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this  *
 * software and associated documentation files (the "Software"), to deal in the Software *
 * without restriction, including without limitation the rights to use, copy, modify,    *
 * merge, publish, distribute, sublicense, and/or sell copies of the Software, and to    *
 * permit persons to whom the Software is furnished to do so, subject to the following   *
 * conditions:                                                                           *
 *                                                                                       *
 * The above copyright notice and this permission notice shall be included in all copies *
 * or substantial portions of the Software.                                              *
 *                                                                                       *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,   *
 * INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A         *
 * PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT    *
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION     *
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE        *
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                                *
 *                                                                                       *
 * Authors:                                                                              *
 *      Luca Pitzalis (luca.pitzalis94@unica.it)                                         *
 *      http://cg3hci.dmi.unica.it/lab/en/people/pitzalis.luca/                          *
 *                                                                                       *
 *      Marco Livesu (marco.livesu@ge.imati.cnr.it)                                      *
 *      http://pers.ge.imati.cnr.it/livesu/                                              *
 *                                                                                       *
 *      Gianmarco Cherchi (g.cherchi@unica.it)                                           *
 *      https://people.unica.it/gianmarcocherchi/                                        *
 *                                                                                       *
 *      Encirco Gobbetti (gobbetti@crs4.it)                                              *
 *      https://www.crs4.it/it/peopledetails/8/enrico-gobbetti/                          *
 *                                                                                       *
 *      Riccardo Scateni (riccardo@unica.it)                                             *
 *      https://people.unica.it/riccardoscateni/                                         *
 *                                                                                       *
 * ***************************************************************************************/

#ifndef ADAPTIVE_GRID_H
#define ADAPTIVE_GRID_H

#include "utils.h"
#include <cinolib/octree.h>


enum directions{TNW=7, TNE=6 , TSW=4, TSE=5, BNW=3, BNE=2, BSW=0, BSE=1};

struct AdaptiveGridAttributes : cinolib::Polyhedron_std_attributes{

    int father;
    std::vector<uint> children;
    bool is_leaf;
    bool is_inside_polycube;
};

typedef cinolib::Hexmesh<cinolib::Mesh_std_attributes, cinolib::Vert_std_attributes, cinolib::Edge_std_attributes, cinolib::Polygon_std_attributes, AdaptiveGridAttributes> AdaptiveGrid;


static std::map<cinolib::vec3d, uint, vert_compare> v_map;

void build_empty_bbox(AdaptiveGrid &grid, const cinolib::AABB &bbox);
void build_empty(AdaptiveGrid &grid, double size);
//void build(AdaptiveGrid &grid, const cinolib::Hexmesh<> &m_tree);
void build(AdaptiveGrid &grid, const cinolib::Trimesh<> &m, bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Trimesh<> &), uint min_depth=0, uint max_depth=8);
//void build(AdaptiveGrid &grid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m,  bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Tetmesh<> &, const cinolib::Tetmesh<> &), uint max_depth=8, uint min_depth=3);
void build_polycube(AdaptiveGrid &grid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m,  bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Tetmesh<> &, const cinolib::Tetmesh<> &), uint min_depth = 5, uint max_depth=8);


void split_cell(AdaptiveGrid &grid, uint pid);
std::vector<uint> sibilings(const AdaptiveGrid &grid, uint pid);

void balancing(AdaptiveGrid &grid);
void strong_balancing(AdaptiveGrid &grid);
void pairing(AdaptiveGrid &grid);
bool children_are_paired(const AdaptiveGrid &grid, uint pid);
bool test_pairing(const AdaptiveGrid &grid);
bool test_balancing(const AdaptiveGrid &grid);
bool test_strong_balancing(const AdaptiveGrid &grid);

int get_neighbor_of_greater_or_equal_size(const AdaptiveGrid &grid, uint pid, uint direction);
void extract_submesh_for_refinement(const AdaptiveGrid &grid, cinolib::Hexmesh<> &submesh, uint refinement, std::unordered_map<uint, uint> &poly_g2s, std::unordered_map<uint, uint> &poly_s2g);
void extract_non_conformal_hexmesh(const AdaptiveGrid &grid, cinolib::Hexmesh<> &hm, std::unordered_map<uint, uint> &poly_map);
void refine_minors(AdaptiveGrid &grid, uint pid, const std::vector<uint> &offs, uint target_depth);
int find_adj_hanging_vert(const AdaptiveGrid &grid, uint vid, uint &v1, uint &v2);
bool vert_is_on_border(const AdaptiveGrid &grid, uint vid);
void clear(AdaptiveGrid &grid);
uint num_leaves(const AdaptiveGrid &grid);
std::pair<uint, uint> get_min_and_max_refinement(const AdaptiveGrid &grid);
std::vector<uint> neighbors(const AdaptiveGrid &grid, uint pid);

#include "adaptive_grid.cpp"

#endif // ADAPTIVE_GRID_H
