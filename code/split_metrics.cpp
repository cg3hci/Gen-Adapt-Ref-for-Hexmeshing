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

#include "split_metrics.h"

//SURFACE

bool sdf_split_metric(const AdaptiveGrid &grid, uint pid, const cinolib::Trimesh<> &m){

    static cinolib::Octree tree = [&m](){
            cinolib::Octree o;
            o.build_from_mesh_polys(m);
            return o;
    }();

    static std::vector<double> sdf = [&m](){
        std::vector<double> sdf_values;
        get_sdf_values(m, sdf_values);
        return sdf_values;
    }();

    std::unordered_set<uint> faces_in_cell;
    tree.intersects_box(grid.poly_aabb(pid), faces_in_cell);
    double min_sdf=std::numeric_limits<double>::infinity();
    for(uint fid : faces_in_cell){
        for(uint eid : m.adj_p2e(fid)){
            if(m.edge_length(eid) < 1e-5)
                sdf[fid] = std::numeric_limits<double>::infinity();
        }
        min_sdf = sdf[fid] < min_sdf ? sdf[fid] : min_sdf;
    }

    std::vector<double> deltas = {grid.poly_aabb(pid).delta_x(), grid.poly_aabb(pid).delta_y(), grid.poly_aabb(pid).delta_z()};
    double max_delta = *std::max_element(deltas.begin(), deltas.end());

    return (max_delta > min_sdf*0.5);

}


//POLYCUBE

bool polycube_distortion_split_metric(const AdaptiveGrid &grid, uint pid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m){

    static cinolib::Octree tet_tree = [&pc_m](){
            cinolib::Octree o;
            o.build_from_mesh_polys(pc_m);
            return o;
    }();

    static cinolib::Octree tet_tree_orig = [&m](){
            cinolib::Octree o;
            o.build_from_mesh_polys(m);
            return o;
    }();

    static cinolib::Octree tri_tree = [&pc_m](){
        cinolib::Octree o;
        for(uint fid=0; fid<pc_m.num_faces(); fid++){
            if(pc_m.face_is_on_srf(fid)){
                o.push_triangle(fid, pc_m.face_verts(fid));
            }
        }
        o.build();
        return o;
    }();

    bool cell_is_on_surf = false;
    for(uint vid : grid.poly_verts_id(pid)){
        if(vert_is_on_border(grid, vid)){
            cell_is_on_surf = true;
            break;
        }
    }

    if(!cell_is_on_surf) return false;

    std::vector<cinolib::vec3d> new_coords;
    double tmp_w[4];


    for(uint idv : grid.poly_verts_id(pid)){
        bool not_found=false;
        if(grid.vert_is_on_srf(idv)){						 //GET INTERNAL VERTICES COORDS

            uint query = 0;
            cinolib::vec3d p;
            double dist;

            tri_tree.closest_point(grid.vert(idv), query,p,dist);

            if(true){
                std::vector<uint> v = pc_m.adj_f2v(query);

                cinolib::triangle_barycentric_coords(pc_m.vert(v[0]), pc_m.vert(v[1]), pc_m.vert(v[2]), p, tmp_w);

                std::vector<uint> v_orig = m.adj_f2v(query);
                new_coords.push_back(m.vert(v_orig[0])*tmp_w[0] + m.vert(v_orig[1])*tmp_w[1] + m.vert(v_orig[2])*tmp_w[2]);
            }
            else{

                std::cerr<<"Surface_vtx not found "<<grid.vert(idv)<<std::endl;
                not_found=true;
            }
        }
        if(!grid.vert_is_on_srf(idv) || not_found){

            uint query = 0;
            cinolib::vec3d p;
            double dist;

            tet_tree.closest_point(grid.vert(idv), query, p, dist);

            if(true){
                std::vector<uint> v = pc_m.adj_p2v(query);

                cinolib::tet_barycentric_coords(pc_m.vert(v[0]), pc_m.vert(v[1]), pc_m.vert(v[2]), pc_m.vert(v[3]), p, tmp_w);

                std::vector<uint> v_orig = m.adj_p2v(query);
                new_coords.push_back(m.vert(v_orig[0])*tmp_w[0] + m.vert(v_orig[1])*tmp_w[1] + m.vert(v_orig[2])*tmp_w[2] + m.vert(v_orig[3])*tmp_w[3]);

            }
            else{

                std::cerr<<"Internal vertex not found into tetmesh"<<std::endl;

            }
        }
    }

    std::vector<uint> tets;
    std::vector<cinolib::vec3d> points;
    for(auto &vert : new_coords){
        uint query = 0;
        cinolib::vec3d p;
        double dist;
        tet_tree_orig.closest_point(vert, query,p,dist);
        tets.push_back(query);
        points.push_back(p);
    }

    std::vector<uint> hex_tets;
    cinolib::hex_to_tets(grid.poly_verts_id(pid),hex_tets);



    for(uint vid=0; vid<m.num_verts(); vid++){
        for(uint i=0; i<hex_tets.size(); i+=4){
            std::vector<uint> tet = {hex_tets[i], hex_tets[i+1],hex_tets[i+2],hex_tets[i+3]};
            if(cinolib::point_in_tet(pc_m.vert(vid), grid.vert(tet[0]), grid.vert(tet[1]), grid.vert(tet[2]), grid.vert(tet[3])) != 0){

                for(uint tid : tets){
                    uint shortest = m.num_edges();
                    for(uint v : m.poly_verts_id(tid)){
                        std::vector<uint> path;
                        cinolib::dijkstra(m, vid, v, path);
                        shortest = std::min(static_cast<uint>(path.size()), shortest);
                    }
                    if(shortest > 2) return true;

                }
                break;
            }
        }
    }


   return false;
}
