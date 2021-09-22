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

#include "adaptive_grid.h"

void build_empty(AdaptiveGrid &grid, double size)
{
    double step = size*0.5;
    std::vector<cinolib::vec3d> verts = {cinolib::vec3d(-step, -step, -step),
                                         cinolib::vec3d(step, -step, -step),
                                         cinolib::vec3d(step, step, -step),
                                         cinolib::vec3d(-step, step, -step),
                                         cinolib::vec3d(-step, -step, step),
                                         cinolib::vec3d(step, -step, step),
                                         cinolib::vec3d(step, step, step),
                                         cinolib::vec3d(-step, step, step)};

    std::vector<std::vector<uint>> polys = {{4,5,1,0,7,6,2,3}};

    grid.clear();
    grid.init(verts, polys);
    grid.poly_data(0).label = 0;
    grid.poly_data(0).father = -1;
    grid.poly_data(0).is_leaf = true;

}

void build_empty_bbox(AdaptiveGrid &grid, const cinolib::AABB &bbox){

    std::vector<cinolib::vec3d> verts = bbox.corners();
    std::vector<std::vector<uint>> polys = {{4,5,1,0,7,6,2,3}};

    grid.clear();
    grid.init(verts, polys);
    grid.poly_data(0).label = 0;
    grid.poly_data(0).father = -1;
    grid.poly_data(0).is_leaf = true;

}

void build(AdaptiveGrid &grid, const cinolib::Trimesh<> &m,  bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Trimesh<> &), uint min_depth, uint max_depth){
    std::vector<double> deltas = {m.bbox().delta_x(), m.bbox().delta_y(), m.bbox().delta_z()};
    double max_delta = *std::max_element(deltas.begin(), deltas.end());

    build_empty(grid, max_delta);
    auto bbox_center =  m.bbox().center();
    for(auto &vert : grid.vector_verts()){
        vert +=bbox_center;
    }
    grid.update_bbox();

    std::deque<uint> queue;
    queue.push_back(0);

    while(!queue.empty()){
        uint curr = queue.front();
        queue.pop_front();
        if(	(grid.poly_data(curr).label < static_cast<int>(min_depth) || split_metric(grid, curr, m))
 && grid.poly_data(curr).label < static_cast<int>(max_depth)){
            split_cell(grid, curr);
            for(uint child : grid.poly_data(curr).children){
                queue.push_back(child);
            }
        }
    }
}

/*void build(AdaptiveGrid &grid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m,  bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Tetmesh<> &, const cinolib::Tetmesh<> &), uint max_depth, uint min_depth){
    std::vector<double> deltas = {pc_m.bbox().delta_x(), pc_m.bbox().delta_y(), pc_m.bbox().delta_z()};
    double max_delta = *std::max_element(deltas.begin(), deltas.end());

    build_empty(grid, max_delta);
    auto bbox_center =  pc_m.bbox().center();
    for(auto &vert : grid.vector_verts()){
        vert +=bbox_center;
    }
    grid.update_bbox();
    //build_empty_bbox(grid, pc_m.bbox());

    std::deque<uint> queue;
    queue.push_back(0);

    while(!queue.empty()){
        uint curr = queue.front();
        queue.pop_front();

        if((grid.poly_data(curr).label < static_cast<int>(min_depth) || split_metric(grid, curr, m, pc_m))  && grid.poly_data(curr).label < static_cast<int>(max_depth)){
            split_cell(grid, curr);
            for(uint child : grid.poly_data(curr).children){
                queue.push_back(child);
            }
        }
    }
}*/


void build_polycube(AdaptiveGrid &grid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m,  bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Tetmesh<> &, const cinolib::Tetmesh<> &), uint min_depth, uint max_depth){
    std::vector<double> deltas = {pc_m.bbox().delta_x(), pc_m.bbox().delta_y(), pc_m.bbox().delta_z()};

    double max_delta = *std::max_element(deltas.begin(), deltas.end());
    //double step = max_delta/initial_size;

    build_empty(grid, max_delta);
    const cinolib::vec3d &starting_vert = pc_m.bbox().min;
    auto s = grid.vert(0);
    for(auto &vert : grid.vector_verts()){
        vert -=s;
    }    for(auto &vert : grid.vector_verts()){
        vert +=starting_vert;
    }
    grid.update_bbox();


    std::deque<uint> queue;
    for(uint pid=0; pid<grid.num_polys();pid++) queue.push_back(pid);

    while(!queue.empty()){
        uint curr = queue.front();
        queue.pop_front();

        if((grid.poly_data(curr).label < min_depth || split_metric(grid, curr, m, pc_m))  && grid.poly_data(curr).label < static_cast<int>(max_depth)){
            split_cell(grid, curr);
            for(uint child : grid.poly_data(curr).children){
                queue.push_back(child);
            }
        }
    }

    cinolib::Octree o;
    o.build_from_mesh_polys(pc_m);
    for(uint pid=0; pid<grid.num_polys(); pid++){
        if(grid.poly_data(pid).label == min_depth){
            uint id;
            if(o.contains(grid.poly_centroid(pid), false, id)){
                grid.poly_data(pid).is_inside_polycube=true;
            }
        }
    }
}

void split_cell(AdaptiveGrid &grid, uint pid)
{

    grid.poly_data(pid).is_leaf = false;
    grid.poly_data(pid).children.reserve(8);
    cinolib::vec3d center = grid.poly_aabb(pid).center();

    uint centeroid_vid = grid.vert_add(center);

    for(uint vid : grid.poly_verts_id(pid)){
        std::vector<uint> poly;
        poly.reserve(8);
        poly.push_back(centeroid_vid);
        poly.push_back(vid);

        for(uint eid : grid.adj_v2e(vid)){
            if(grid.poly_contains_edge(pid, eid)){
                cinolib::vec3d middle = grid.edge_sample_at(eid, 0.5);
                auto query = v_map.find(middle);
                if(query != v_map.end()){
                    poly.push_back(query->second);
                }
                else{
                    uint fresh_vid = grid.vert_add(middle);
                    poly.push_back(fresh_vid);
                    v_map[middle] = fresh_vid;
                }
            }
        }
        for(uint fid : grid.adj_v2f(vid)){
            if(grid.poly_contains_face(pid, fid)){
                cinolib::vec3d middle = grid.face_centroid(fid);
                auto query = v_map.find(middle);
                if(query != v_map.end()){
                    poly.push_back(query->second);
                }
                else{
                    uint fresh_vid = grid.vert_add(middle);
                    poly.push_back(fresh_vid);
                    v_map[middle] = fresh_vid;
                }
            }
        }

        poly_vert_ordering(grid.vector_verts(), poly);
        uint fresh_pid = grid.poly_add(poly);
        grid.poly_data(pid).children.push_back(fresh_pid);
        grid.poly_data(fresh_pid).father = pid;
        grid.poly_data(fresh_pid).is_leaf = true;
        grid.poly_data(fresh_pid).label = grid.poly_data(pid).label+1;
    }

}



std::vector<uint> sibilings(const AdaptiveGrid &grid, uint pid)
{
    std::vector<uint> sibilings;
    sibilings.reserve(7);
    for(uint sib : grid.poly_data(grid.poly_data(pid).father).children){
        if(pid != sib) sibilings.push_back(sib);
    }

    return sibilings;
}


int get_neighbor_of_greater_or_equal_size(const AdaptiveGrid &grid, uint pid, uint direction){

    if(grid.poly_data(pid).father == -1){
        return -1;
    }

    switch(direction){
    case 0:
    {
        if(grid.poly_data(grid.poly_data(pid).father).children[TNW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TNE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BSE];

        int tmp = get_neighbor_of_greater_or_equal_size(grid, grid.poly_data(pid).father, direction);
        if(tmp == -1 || grid.poly_data(tmp).is_leaf) return tmp;

        if(grid.poly_data(grid.poly_data(pid).father).children[BNW] == pid)
            return grid.poly_data(tmp).children[TNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNE] == pid)
            return grid.poly_data(tmp).children[TNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSW] == pid)
            return grid.poly_data(tmp).children[TSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSE] == pid)
            return grid.poly_data(tmp).children[TSE];
    }
        break;
    case 1:
    {
        if(grid.poly_data(grid.poly_data(pid).father).children[TNW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TSE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BSE];

        int tmp = get_neighbor_of_greater_or_equal_size(grid, grid.poly_data(pid).father, direction);
        if(tmp == -1 || grid.poly_data(tmp).is_leaf) return tmp;

        if(grid.poly_data(grid.poly_data(pid).father).children[TNE] == pid)
            return grid.poly_data(tmp).children[TNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSE] == pid)
            return grid.poly_data(tmp).children[TSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNE] == pid)
            return grid.poly_data(tmp).children[BNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSE] == pid)
            return grid.poly_data(tmp).children[BSW];
    }
        break;
    case 2:
    {
        if(grid.poly_data(grid.poly_data(pid).father).children[BNW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TSE];

        int tmp = get_neighbor_of_greater_or_equal_size(grid, grid.poly_data(pid).father, direction);
        if(tmp == -1 || grid.poly_data(tmp).is_leaf) return tmp;

        if(grid.poly_data(grid.poly_data(pid).father).children[TNW] == pid)
            return grid.poly_data(tmp).children[BNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TNE] == pid)
            return grid.poly_data(tmp).children[BNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSW] == pid)
            return grid.poly_data(tmp).children[BSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSE] == pid)
            return grid.poly_data(tmp).children[BSE];
    }
        break;
    case 3:
    {
        if(grid.poly_data(grid.poly_data(pid).father).children[TNE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BSW];

        int tmp = get_neighbor_of_greater_or_equal_size(grid, grid.poly_data(pid).father, direction);
        if(tmp == -1 || grid.poly_data(tmp).is_leaf) return tmp;

        if(grid.poly_data(grid.poly_data(pid).father).children[TNW] == pid)
            return grid.poly_data(tmp).children[TNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSW] == pid)
            return grid.poly_data(tmp).children[TSE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNW] == pid)
            return grid.poly_data(tmp).children[BNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSW] == pid)
            return grid.poly_data(tmp).children[BSE];
    }
        break;
    case 4:
    {
        if(grid.poly_data(grid.poly_data(pid).father).children[TNW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TNE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TSE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BSE];

        int tmp = get_neighbor_of_greater_or_equal_size(grid, grid.poly_data(pid).father, direction);
        if(tmp == -1 || grid.poly_data(tmp).is_leaf) return tmp;

        if(grid.poly_data(grid.poly_data(pid).father).children[TSW] == pid)
            return grid.poly_data(tmp).children[TNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSE] == pid)
            return grid.poly_data(tmp).children[TNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSW] == pid)
            return grid.poly_data(tmp).children[BNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSE] == pid)
            return grid.poly_data(tmp).children[BNE];
    }
        break;
    case 5:
    {
        if(grid.poly_data(grid.poly_data(pid).father).children[TSW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TSE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[TNE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSW] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BNW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BSE] == pid)
            return grid.poly_data(grid.poly_data(pid).father).children[BNE];

        int tmp = get_neighbor_of_greater_or_equal_size(grid, grid.poly_data(pid).father, direction);
        if(tmp == -1 || grid.poly_data(tmp).is_leaf) return tmp;

        if(grid.poly_data(grid.poly_data(pid).father).children[TNW] == pid)
            return grid.poly_data(tmp).children[TSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[TNE] == pid)
            return grid.poly_data(tmp).children[TSE];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNW] == pid)
            return grid.poly_data(tmp).children[BSW];
        if(grid.poly_data(grid.poly_data(pid).father).children[BNE] == pid)
            return grid.poly_data(tmp).children[BSE];
    }
        break;
    }

    return -1;
}

void strong_balancing(AdaptiveGrid &grid)
{
    std::deque<uint> queue;
    //queue.push_back(0);
    for(uint pid=0; pid<grid.num_polys();pid++){
        if(grid.poly_data(pid).father == -1) queue.push_back(pid);
        else break;
    }

    while(!queue.empty()){

        uint curr = queue.front();
        queue.pop_front();

        for(uint child : grid.poly_data(curr).children){
            if(grid.poly_data(child).is_leaf){

                std::set<uint> neighbors;
                for(uint vid : grid.poly_verts_id(child)){
                    for(uint neigh : grid.adj_v2p(vid)){
                        if((neigh != child) && grid.poly_data(neigh).is_leaf && (grid.poly_data(neigh).label < grid.poly_data(child).label))
                            neighbors.insert(neigh);
                    }
                }
                for(uint neigh : neighbors){
                    int dc = grid.poly_data(child).label;
                    int dn = grid.poly_data(neigh).label;
                    int diff = std::abs(dc-dn);
                    if(diff > 1){
                        if(dc > dn){
                            split_cell(grid, neigh);
                            queue.push_back(neigh);
                        }
                        else{
                            split_cell(grid, child);
                            queue.push_front(child);
                            break;
                        }
                    }
                }
            }
            else{
                queue.push_back(child);
            }
        }
    }
}

void balancing(AdaptiveGrid &grid)
{
   std::deque<uint> queue;
   //queue.push_back(0);
   for(uint pid=0; pid<grid.num_polys();pid++){
       if(grid.poly_data(pid).father == -1) queue.push_back(pid);
       else break;
   }

   while(!queue.empty()){

       uint curr = queue.front();
       queue.pop_front();

       for(uint child : grid.poly_data(curr).children){
           if(grid.poly_data(child).is_leaf){
               for(uint dir=0; dir<6; dir++){
                   int neigh = get_neighbor_of_greater_or_equal_size(grid, child, dir);
                   if(neigh == -1) continue;

                   int dc = grid.poly_data(child).label;
                   int dn = grid.poly_data(neigh).label;
                   int diff = std::abs(dc-dn);
                   if(diff > 1){
                       //std::cout<<"diff "<<dc<<" "<<dn<<std::endl;
                       if(dc > dn){
                           split_cell(grid, neigh);
                           queue.push_back(neigh);
                       }
                       else{
                           split_cell(grid, child);
                           queue.push_front(child);
                           break;
                       }
                   }
               }
           }
           else{
               queue.push_back(child);
           }
       }
   }
}


void pairing(AdaptiveGrid &grid){

    std::deque<uint> queue;
    queue.push_back(0);
    while(!queue.empty()){
        uint curr = queue.front();
        queue.pop_front();
        for(uint child : grid.poly_data(curr).children){

            if(!grid.poly_data(child).is_leaf){
                for(uint sib : sibilings(grid, child)){
                    if(grid.poly_data(sib).is_leaf)
                        split_cell(grid, sib);
                }
                queue.push_back(child);
            }
        }
    }

}


bool children_are_paired(const AdaptiveGrid &grid, uint pid){
    if(grid.poly_data(pid).is_leaf) return true;
    std::set<bool> splitted;
    for(uint child : grid.poly_data(pid).children){
        splitted.insert(!grid.poly_data(child).is_leaf);
    }
    return splitted.size() == 1;
}

bool test_pairing(const AdaptiveGrid &grid){
    std::deque<uint> queue;
    queue.push_back(0);
    bool paired = true;
    while(!queue.empty()){
        uint curr = queue.front();
        queue.pop_front();
        paired = paired && children_are_paired(grid, curr);
        for(uint child : grid.poly_data(curr).children){
                queue.push_back(child);
        }
    }

    return paired;
}

bool test_balancing(const AdaptiveGrid &grid){

    std::deque<uint> queue;
    queue.push_back(0);

    while(!queue.empty()){

        uint curr = queue.front();
        queue.pop_front();

        for(uint child : grid.poly_data(curr).children){
            if(grid.poly_data(child).is_leaf){
                for(uint dir=0; dir<6; dir++){
                    int neigh = get_neighbor_of_greater_or_equal_size(grid, child, dir);
                    if(neigh == -1) continue;

                    int dc = grid.poly_data(child).label;
                    int dn = grid.poly_data(neigh).label;
                    int diff = std::abs(dc-dn);
                    if(diff > 1){
                        return false;
                    }
                }
            }
            else{
                queue.push_back(child);
            }
        }
    }
    return true;
}

bool test_strong_balancing(const AdaptiveGrid &grid){

    std::deque<uint> queue;
    queue.push_back(0);

    while(!queue.empty()){

        uint curr = queue.front();
        queue.pop_front();

        for(uint child : grid.poly_data(curr).children){
            if(grid.poly_data(child).is_leaf){

                std::set<uint> neighbors;
                for(uint vid : grid.poly_verts_id(child)){
                    for(uint neigh : grid.adj_v2p(vid)){
                        if((neigh != child) && grid.poly_data(neigh).is_leaf && (grid.poly_data(neigh).label < grid.poly_data(child).label))
                            neighbors.insert(neigh);
                    }
                }
                for(uint neigh : neighbors){
                    int dc = grid.poly_data(child).label;
                    int dn = grid.poly_data(neigh).label;
                    int diff = std::abs(dc-dn);
                    if(diff > 1){
                        std::cout<<diff<<std::endl;
                        return false;
                    }
                }
            }
            else{
                queue.push_back(child);
            }
        }
    }
    return true;
}

void extract_non_conformal_hexmesh(const AdaptiveGrid &grid, cinolib::Hexmesh<> &hm, std::unordered_map<uint, uint> &poly_map){

    std::vector<int> poly_labels;
    std::vector<std::vector<uint>> polys;
    for(uint pid=0; pid<grid.num_polys(); pid++){
        if(grid.poly_data(pid).is_leaf){
            poly_map[pid] = polys.size();
            polys.push_back(grid.poly_verts_id(pid));
            poly_labels.push_back(grid.poly_data(pid).label);
        }
    }

    hm.clear();
    hm.init(grid.vector_verts(), polys);
    hm.poly_apply_labels(poly_labels);

}


void extract_submesh_for_refinement(const AdaptiveGrid &grid, cinolib::Hexmesh<> &submesh, uint refinement, std::unordered_map<uint, uint> &poly_g2s, std::unordered_map<uint, uint> &poly_s2g){

    std::vector<std::vector<uint>> polys;
    std::vector<int> labels;

    for(uint pid=0; pid<grid.num_polys(); pid++){
        if(static_cast<uint>(grid.poly_data(pid).label) == refinement){
            poly_g2s[pid] = polys.size();
            poly_s2g[polys.size()] = pid;
            polys.push_back(grid.poly_verts_id(pid));
            if(grid.poly_data(pid).is_leaf){
                labels.push_back(0);
            }
            else{
                labels.push_back(1);
            }
        }
    }

    submesh.clear();
    submesh.init(grid.vector_verts(), polys);
    submesh.poly_apply_labels(labels);
}

void refine_minors(AdaptiveGrid &grid, uint pid, const std::vector<uint> &offs, uint target_depth){

    split_cell(grid, pid);

    for(uint off : offs){
       uint next = grid.poly_data(pid).children[off];
       for(int i=grid.poly_data(next).label; i<static_cast<int>(target_depth); i++){
           split_cell(grid, next);
           next = grid.poly_data(next).children[off];
       }
    }
}

int find_adj_hanging_vert(const AdaptiveGrid &grid, uint vid, uint &v1, uint &v2){

    double eps = 1e-6;
    for(uint pid=0; pid<grid.num_polys(); pid++){
        if(grid.poly_data(pid).is_leaf){

            for(uint eid : grid.adj_p2e(pid)){
                cinolib::vec3d middle = grid.edge_sample_at(eid, 0.5);
                double dist = middle.dist(grid.vert(vid));
                if(dist < eps){
                    v1 = grid.edge_vert_id(eid, 0);
                    v2 = grid.edge_vert_id(eid, 1);
                    return pid;
                }
            }

        }
    }



    return -1;
}

bool vert_is_on_border(const AdaptiveGrid &grid, uint vid){

    const cinolib::vec3d &min = grid.bbox().min;
    const cinolib::vec3d &max = grid.bbox().max;

    if(eps_equal(min.x(), grid.vert(vid).x()) || eps_equal(max.x(), grid.vert(vid).x())) return true;
    if(eps_equal(min.y(), grid.vert(vid).y()) || eps_equal(max.y(), grid.vert(vid).y())) return true;
    if(eps_equal(min.z(), grid.vert(vid).z()) || eps_equal(max.z(), grid.vert(vid).z())) return true;

    return false;

}

std::pair<uint, uint> get_min_and_max_refinement(const AdaptiveGrid &grid){

    int max = 0;
    int min = 1e5;
    for(uint pid=0; pid<grid.num_polys(); pid++){
        if(grid.poly_data(pid).is_leaf){
            if(grid.poly_data(pid).label > max)
                max = grid.poly_data(pid).label;
            if(grid.poly_data(pid).label < min)
                min = grid.poly_data(pid).label;
        }
    }

    return std::make_pair(min, max);
}

std::vector<uint> neighbors(const AdaptiveGrid &grid, uint pid){

    std::vector<uint> neighs;
    for(uint dir=0; dir<6; dir++){
        int neigh = get_neighbor_of_greater_or_equal_size(grid, pid, dir);
        if(neigh != -1 && grid.poly_data(pid).is_leaf){
            neighs.push_back(neigh);
        }
    }
    return neighs;
}

uint num_leaves(const AdaptiveGrid &grid)
{
   uint n_leaves=0;
   for(uint pid=0; pid<grid.num_polys(); pid++){
       if(grid.poly_data(pid).is_leaf){
           n_leaves++;
       }
   }
   return n_leaves;
}

void clear(AdaptiveGrid &grid)
{
   grid.clear();
   v_map.clear();
}


void build(AdaptiveGrid &grid, const cinolib::Hexmesh<> &m_tree)
{
    std::vector<double> deltas = {m_tree.bbox().delta_x(), m_tree.bbox().delta_y(), m_tree.bbox().delta_z()};
    double max_delta = *std::max_element(deltas.begin(), deltas.end());

    build_empty(grid, max_delta);
    auto bbox_center =  m_tree.bbox().center();
    for(auto &vert : grid.vector_verts()){
        vert +=bbox_center;
    }
    grid.update_bbox();

    std::set<cinolib::vec3d, vert_compare> vert_set;
    for(const auto &vert : m_tree.vector_verts()){
        vert_set.insert(vert);
    }

    std::deque<uint> queue;
    queue.push_back(0);

    while(!queue.empty()){
        uint curr = queue.front();
        queue.pop_front();

        cinolib::vec3d poly_centroid = grid.poly_centroid(curr);
        if(cinolib::CONTAINS(vert_set, poly_centroid)){
            split_cell(grid, curr);
            for(uint child : grid.poly_data(curr).children){
                queue.push_back(child);
            }
        }

    }

}
