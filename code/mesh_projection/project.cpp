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

#include "project.h"

void remove_external_polys(const cinolib::Trimesh<> &m, cinolib::Hexmesh<> &hm){

    std::cout<<"Removing external polyhedra...."<<std::endl;
    std::vector<bool> polys_to_delete_bool(hm.num_polys());
    std::vector<uint> polys_to_delete;
    cinolib::PARALLEL_FOR(0, hm.num_polys(),1000, [&polys_to_delete_bool, &hm, &m](uint pid){
        cinolib::vec3d centroid = hm.poly_centroid(pid);
        if(cinolib::winding_number(m,centroid) == 0){
            polys_to_delete_bool[pid] = true;

        }

    });

    for(uint pid=0; pid<hm.num_polys(); pid++)
        if(polys_to_delete_bool[pid]) polys_to_delete.push_back(pid);

    cinolib::Hexmesh<> tmp_hm = hm;
    tmp_hm.polys_remove(polys_to_delete);

    std::map<cinolib::vec3d, uint, vert_comparator> vmap;

    for(uint vid=0; vid<hm.num_verts(); vid++) vmap[hm.vert(vid)] = vid;

    for(uint eid=0; eid<tmp_hm.num_edges(); eid++){
        if(!tmp_hm.edge_is_manifold(eid)){

            uint eid_ = hm.edge_id(vmap[tmp_hm.edge_vert(eid, 0)], vmap[tmp_hm.edge_vert(eid, 1)]);
            for(uint pid : hm.adj_e2p(eid_)){
                polys_to_delete_bool[pid] = false;
            }

        }
    }

    for(uint vid=0; vid<tmp_hm.num_verts(); vid++){
        if(!tmp_hm.vert_is_manifold(vid)){

            uint vid_ = vmap[tmp_hm.vert(vid)];
            for(uint pid : hm.adj_v2p(vid_)){
                polys_to_delete_bool[pid] = false;
            }
        }
    }
    polys_to_delete.clear();
    for(uint pid=0; pid<hm.num_polys(); pid++){
        if(polys_to_delete_bool[pid]){
            polys_to_delete.push_back(pid);
            hm.poly_data(pid).label = 1;
        }
        else hm.poly_data(pid).label = 0;
    }


    hm.polys_remove(polys_to_delete);
}

void project_on_surface_smoother(cinolib::Hexmesh<> &hm, cinolib::Trimesh<> &target, cinolib::Quadmesh<> &peel){


    std::unordered_map<uint, uint> m_map;
    std::unordered_map<uint, uint> v_map;
    cinolib::export_surface(hm, peel, m_map, v_map);
    target.edge_set_flag(cinolib::MARKED, false);
    peel.edge_set_flag(cinolib::MARKED, false);
    cinolib::Quadmesh tmp_peel = peel;

    target.edge_mark_sharp_creases();
    cinolib::Octree tree;
    tree.build_from_mesh_edges(peel);

    uint num_samples = 10;
    for(uint eid=0; eid<target.num_edges(); eid++){
        break;
        if(target.edge_data(eid).flags[cinolib::MARKED]){

            std::vector<cinolib::vec3d> samples;
            for(uint i=0; i<=num_samples; i++){
                samples.push_back(target.edge_sample_at(eid, static_cast<double>(i)/static_cast<double>(num_samples)));
            }


            for(uint i=0; i<samples.size()-1; i++){
                uint id1, id2;
                double dist1, dist2;
                cinolib::vec3d pos1, pos2;

                tree.closest_point(samples[i], id1, pos1, dist1);
                tree.closest_point(samples[i+1], id2, pos2, dist2);

                uint source = 0;
                if(pos1.dist(peel.edge_vert(id1, 0)) < pos1.dist(peel.edge_vert(id1, 1)))
                    source = peel.edge_vert_id(id1, 0);
                else
                    source = peel.edge_vert_id(id1, 1);
                uint dest = 0;
                if(pos2.dist(peel.edge_vert(id2, 0)) < pos2.dist(peel.edge_vert(id2, 1)))
                    dest = peel.edge_vert_id(id2, 0);
                else
                    dest = peel.edge_vert_id(id2, 1);

                std::vector<uint> path;
                cinolib::dijkstra(peel, source, dest, path);
                if(path.size() == 0) continue;
                for(uint j=0; j<path.size()-1; j++){
                    int edge_to_mark = peel.edge_id(path[j], path[j+1]);
                    if(edge_to_mark < 0) std::cerr<<"This edge doen't exist"<<std::endl;
                    else{
                        peel.edge_set_flag(cinolib::MARKED, true, {static_cast<uint>(edge_to_mark)});
                    }
                }

            }


        }
    }
    cinolib::SmootherOptions options;
    options.n_iters = 3;
    options.laplacian_mode = cinolib::UNIFORM;
    cinolib::mesh_smoother(peel, target, options);


    std::cout<<"Projecting..."<<std::endl;
    for(uint vid=0; vid<peel.num_verts(); vid++){
        hm.vert(v_map[vid]) = peel.vert(vid);
    }


}

void project_on_surface_grid_projector(cinolib::Hexmesh<> &hm, cinolib::Trimesh<> &target){

    std::cout<<"Starting projection"<<std::endl;
    std::vector<std::vector<uint>> fn_tri, fn_quad;
    cinolib::feature_network(target, fn_tri);

    cinolib::Quadmesh<> peel;
    std::unordered_map<uint, uint> h2q, q2h;
    cinolib::export_surface(hm, peel, h2q, q2h);

    cinolib::feature_mapping(target, fn_tri, peel, fn_quad);

    for(auto &f : fn_quad){

        for(uint i=1; i<f.size(); ++i){
            uint v0 = q2h.at(f.at(i));
            uint v1 = q2h.at(f.at(i-1));
            int eid = hm.edge_id(v0, v1);

            assert(eid >= 0 );

            hm.edge_data(eid).flags[cinolib::CREASE] = true;
            hm.edge_data(eid).flags[cinolib::MARKED] = true;
        }
    }

    cinolib::GridProjectorOptions options;
    //options.SJ_thresh = -1;
    cinolib::grid_projector(hm, target, options);

    std::cout<<"Projection completed"<<std::endl;

}

void polycube_projector(cinolib::Hexmesh<> &grid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m){
 
    //:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    cinolib::Octree tet_tree(9);
    cinolib::Octree tri_tree(9);
    tet_tree.build_from_mesh_polys(pc_m);
    for(uint fid=0; fid<pc_m.num_faces(); fid++){
        if(pc_m.face_is_on_srf(fid)){
            tri_tree.push_triangle(fid, pc_m.face_verts(fid));
        }
    }
    tri_tree.build();


    std::vector<cinolib::vec3d> new_coords(grid.num_verts());
    double tmp_w[4];


    for(uint idv=0; idv<grid.num_verts(); idv++){
        if(grid.vert_is_on_srf(idv)){

            uint query = 0;
            cinolib::vec3d p;
            double dist;

            tri_tree.closest_point(grid.vert(idv), query,p,dist);

            std::vector<uint> v = pc_m.adj_f2v(query);

            cinolib::triangle_barycentric_coords(pc_m.vert(v[0]), pc_m.vert(v[1]), pc_m.vert(v[2]), p, tmp_w);

            std::vector<uint> v_orig = m.adj_f2v(query);
            new_coords[idv] = m.vert(v_orig[0])*tmp_w[0] + m.vert(v_orig[1])*tmp_w[1] + m.vert(v_orig[2])*tmp_w[2];
        }
        else{

            uint query = 0;
            cinolib::vec3d p;
            double dist;

            tet_tree.closest_point(grid.vert(idv), query,p,dist);

            std::vector<uint> v = pc_m.adj_p2v(query);

            cinolib::tet_barycentric_coords(pc_m.vert(v[0]), pc_m.vert(v[1]), pc_m.vert(v[2]), pc_m.vert(v[3]), p, tmp_w);

            std::vector<uint> v_orig = m.adj_p2v(query);
            new_coords[idv] = m.vert(v_orig[0])*tmp_w[0] + m.vert(v_orig[1])*tmp_w[1] + m.vert(v_orig[2])*tmp_w[2] + m.vert(v_orig[3])*tmp_w[3];

            if(new_coords[idv].is_inf()){
                new_coords[idv] = m.poly_centroid(query);
            }


        }

    }


    grid.vector_verts() = new_coords;

}


void project_mesh(cinolib::Hexmesh<> &grid, cinolib::Trimesh<> &target, bool use_grid_projector){
    
    remove_external_polys(target, grid);
    cinolib::padding(grid, true);

    cinolib::Quadmesh<> peel;

    if(!use_grid_projector) project_on_surface_smoother(grid, target, peel);
    else project_on_surface_grid_projector(grid, target);
}

void project_mesh(cinolib::Hexmesh<> &grid, const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m){
    
    cinolib::padding(grid, true);

    polycube_projector(grid, m, pc_m);

}