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

#include "find_transition_verts.h"


void extract_submesh_for_refinement(cinolib::AbstractPolyhedralMesh<> &m, cinolib::Hexmesh<> &sub_m, int ref){

    std::vector<std::vector<uint>> polys;

    for(uint pid=0; pid<m.num_polys(); pid++){
        if(m.poly_data(pid).label == ref)
            polys.push_back(m.poly_verts_id(pid));
    }

    sub_m.clear();
    std::vector<int> poly_labels(m.num_polys(), 0);
    sub_m.init(m.vector_verts(), polys, m.vector_vert_labels(), poly_labels);
}

void extract_quadmesh_for_refinement(const cinolib::AbstractPolyhedralMesh<> &m, uint ref, cinolib::Quadmesh<> &qm){

    std::vector<std::vector<uint>> polys;
    std::map<cinolib::vec3d, uint, vert_comparator> v_map;
    for(uint vid=0; vid<m.num_verts(); vid++) v_map[m.vert(vid)] = vid;

    for(uint pid=0; pid<m.num_polys();pid++){

        if(m.poly_data(pid).label == static_cast<int>(ref)){
            for(uint fid : m.adj_p2f(pid)){
                if(m.face_is_on_srf(fid)){
                    auto centroid = m.face_centroid(fid);
                    if(v_map.find(centroid) != v_map.end()){
                        polys.push_back(m.face_verts_id(fid));
                    }
                }
            }
        }
    }
    qm.clear();
    qm.init(m.vector_verts(), polys);

}

void find_t_vertices(cinolib::AbstractPolyhedralMesh<> &m, uint div_factor){
    
    std::vector<uint> solutions(m.num_verts());
    auto labels = m.vector_poly_labels();

    const std::vector<int> &refs = m.vector_poly_labels();
    int max_ref = *std::max_element(refs.begin(), refs.end());
    int min_ref = *std::min_element(refs.begin(), refs.end());

    std::cout<<max_ref<<" "<<min_ref<<std::endl;
    for(int ref = max_ref-1; ref>=min_ref; --ref){
        std::vector<std::pair<cinolib::Quadmesh<>, uint>> to_check;
        cinolib::Quadmesh<> submesh;
        extract_quadmesh_for_refinement(m, ref, submesh);


    try
    {
        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 0);

        GRBModel model = GRBModel(env);

        // variables
        //std::vector<GRBVar> vert_variables(m.num_verts());
        GRBVar  *vert_variables   = model.addVars(submesh.num_verts(), GRB_BINARY);


        /*for(uint vid = 0; vid < m.num_verts(); vid++){
            vert_variables[vid] = model.addVar(0, 1, 1,GRB_BINARY);
        }*/

        // objective function
        GRBLinExpr obj = 0;


        for(uint vid=0; vid<submesh.num_verts(); vid++){
            double mult = 1;

            if(!submesh.vert_is_boundary(vid)){
                mult = submesh.num_verts();
            }
            else{
                if(submesh.adj_v2e(vid).size() == 3){
                    mult = static_cast<double>(submesh.num_verts())/2;
                }
                else{
                    mult = static_cast<double>(submesh.num_verts())/div_factor;
                }
            }

            obj += vert_variables[vid] * mult;
        }


        model.setObjective(obj, GRB_MAXIMIZE);


        for(uint vid=0; vid<submesh.num_verts(); vid++){

            if(submesh.vert_valence(vid) == 0){
                model.addConstr(vert_variables[vid] == 0);
            }
            else{

                std::set<cinolib::vec3d, vert_comparator> normals;
                for(uint fid : submesh.adj_v2p(vid)){
                    normals.insert(submesh.poly_data(fid).normal);
                }
                if(normals.size() >= 2){
                     model.addConstr(vert_variables[vid] == 0);
                     continue;
                }

                std::set<uint> vs;
                for(uint fid : submesh.adj_v2p(vid)){
                    for(uint v : submesh.poly_verts_id(fid)){
                        if(vid != v)
                            vs.insert(v);
                    }
                }
                for(uint v : vs){
                    model.addConstr(vert_variables[vid] + vert_variables[v] <= 1);
                }
            }


        }

        for(uint pid=0; pid<submesh.num_polys(); ++pid)
        {
            for(uint nbr : submesh.adj_p2p(pid))
            {
                uint eid  = submesh.edge_shared(pid,nbr);
                if(!submesh.edge_is_manifold(eid)){
                    const auto &v_list_pid = submesh.poly_verts_id(pid, true);
                    const auto &v_list_nbr = submesh.poly_verts_id(nbr, true);

                    uint fid_pid = m.face_id(v_list_pid);
                    uint fid_nbr = m.face_id(v_list_nbr);

                    auto adj_pids1 = m.adj_f2p(fid_pid);
                    auto adj_pids2 = m.adj_f2p(fid_nbr);

                    std::sort(adj_pids1.begin(), adj_pids1.end());
                    std::sort(adj_pids2.begin(), adj_pids2.end());

                    std::vector<uint> intersection;

                    std::set_intersection(adj_pids1.begin(), adj_pids1.end(), adj_pids2.begin(), adj_pids2.end(), std::back_inserter(intersection));

                    if(intersection.size() == 0)
                        continue;
                }
                uint vid0 = submesh.edge_vert_id(eid,0);
                uint vid1 = submesh.edge_vert_id(eid,1);
                uint off_pid  = submesh.poly_vert_offset(pid, vid0);
                if(submesh.poly_vert_id(pid,(off_pid+1)%4)==vid1) off_pid = submesh.poly_vert_offset(pid, vid1);
                uint off_nbr  = submesh.poly_vert_offset(nbr, vid0);
                if(submesh.poly_vert_id(nbr,(off_nbr+1)%4)==vid1) off_nbr = submesh.poly_vert_offset(nbr, vid1);
                uint pid_next = submesh.poly_vert_id(pid,(off_pid+1)%4);
                uint nbr_next = submesh.poly_vert_id(nbr,(off_nbr+1)%4);
                model.addConstr(vert_variables[pid_next] + vert_variables[nbr_next] <= 1);
                pid_next = submesh.poly_vert_id(pid,(off_pid+2)%4);
                nbr_next = submesh.poly_vert_id(nbr,(off_nbr+2)%4);
                model.addConstr(vert_variables[pid_next] + vert_variables[nbr_next] <= 1);

            }
        }



        model.optimize();

        // solution parsing
        for(uint vid = 0; vid < submesh.num_verts(); vid++){

            if(submesh.vert_valence(vid)==0) continue;

            double var = std::round(vert_variables[vid].get(GRB_DoubleAttr_X));
            //submesh.vert_data(vid).label += static_cast<uint>(var);
            solutions[vid] = static_cast<uint>(var);
        }



    }
    catch(GRBException e)
    {
        std::cerr << "Error code" << e.getErrorCode() <<std::endl;
        std::cerr << e.getMessage().c_str() <<std::endl;
    }
    catch(...)
    {
        std::cerr << "Exception during optimization" << std::endl;
    }

    }
    for(uint vid = 0; vid < m.num_verts(); vid++){
        if(solutions[vid] == 1)
            m.vert_data(vid).label = 2;
    }

}