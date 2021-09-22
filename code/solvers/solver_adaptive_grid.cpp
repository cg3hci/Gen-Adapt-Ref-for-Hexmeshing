#include "solver_adaptive_grid.h"

bool weak_balancing = true;

void solve(AdaptiveGrid &grid, bool use_weak_balancing){
    weak_balancing = use_weak_balancing;
    auto min_max = get_min_and_max_refinement(grid);

    std::cout<<min_max.first<<" "<<min_max.second<<std::endl;
    if(weak_balancing) balancing(grid);
    else strong_balancing(grid);
    for(int ref=static_cast<int>(min_max.second)-1; ref>=static_cast<int>(min_max.first); --ref){

        //std::cerr<<"REF: "<<ref<<std::endl;
        cinolib::Hexmesh<> submesh;
        std::unordered_map<uint, uint> poly_g2s;
        std::unordered_map<uint, uint> poly_s2g;
        extract_submesh_for_refinement(grid, submesh, ref, poly_g2s, poly_s2g);

        perform_pairing(grid, submesh, ref, poly_g2s, poly_s2g);
        if(weak_balancing) balancing(grid);
        else strong_balancing(grid);
    }
}

void perform_pairing(AdaptiveGrid &grid, const cinolib::Hexmesh<> &submesh, uint refinement, std::unordered_map<uint, uint> &poly_g2s, std::unordered_map<uint, uint> &poly_s2g){

    try{
        auto start_construction = startChrono();

        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 0);

        GRBModel model = GRBModel(env);

        // variables
        GRBVar* vert_variables = model.addVars(submesh.num_verts(), GRB_BINARY);

        // objective function
        std::vector<double> dists(grid.num_verts(), 0);
        //compute_weights_per_vertex(grid, submesh, dists);
        GRBLinExpr obj = 0;
        for(uint pid = 0; pid < submesh.num_polys(); pid++)
        {
            GRBLinExpr sub_obj = 0;
            int min = std::min(submesh.poly_data(pid).label, 1);

            for(uint vid : submesh.poly_verts_id(pid))
            {
                sub_obj += vert_variables[vid];

            }

            obj += sub_obj - min;
        }
        GRBLinExpr sum = 0;
        for(uint vid=0; vid<submesh.num_verts(); vid++){
            break;
            if(submesh.vert_valence(vid) == 0) continue;

            uint mult=0;
            uint abs_G_v = 0;
            std::unordered_set<uint> disjoint_elements;
            bool vert_polys_touch_srf = false;
            for(uint pid : submesh.adj_v2p(vid)){
                if(submesh.poly_is_on_surf(pid)) {
                    vert_polys_touch_srf = true;
                    break;
                }
            }
            if(vert_polys_touch_srf){

                for(uint adj_sub : submesh.adj_v2p(vid)){
                    for(uint adj_global : neighbors(grid, poly_s2g[adj_sub])){
                        auto query = poly_g2s.find(adj_global);
                        if(query == poly_g2s.end() && grid.poly_data(adj_global).is_leaf)
                            disjoint_elements.insert(adj_global);
                    }
                }

                abs_G_v = disjoint_elements.size();
                if(abs_G_v != 0){
                    mult = abs_G_v + dists[vid];
                }
            }

            sum += vert_variables[vid] * mult;
        }
        //obj+=sum;
        model.setObjective(obj);


        //Constraints

        //Grant minimum refinement
        for(uint pid = 0; pid < submesh.num_polys(); pid++)
        {
            int min = std::min(submesh.poly_data(pid).label, 1);
            GRBLinExpr expr = 0;
            for(uint vid : submesh.poly_verts_id(pid))
            {
                expr += vert_variables[vid];
            }
            model.addConstr(expr >= min);
        }

        //Avoid misaligned multi resolution concaves
        if(weak_balancing){
            for(uint vid=0; vid<submesh.num_verts(); vid++){
                if(submesh.vert_valence(vid) == 0 || (submesh.vert_is_on_srf(vid) && !vert_is_on_border(grid, vid))) continue;

                std::set<uint> verts_to_check;
                for(uint fid : submesh.adj_v2f(vid)){
                    for(uint v : submesh.face_verts_id(fid)){
                        if( v!= vid){
                            if(submesh.edge_id(vid, v) == -1 && submesh.vert_is_on_srf(v)){
                                verts_to_check.insert(v);
                                break;
                            }
                        }
                    }
                }

                for(uint adj_v : verts_to_check){
                    if(submesh.vert_is_on_srf(adj_v)){
                        for(uint pid : grid.adj_v2p(adj_v)){
                            if(grid.poly_data(pid).is_leaf){
                                if(grid.poly_data(pid).label == static_cast<int>(refinement)-1){
                                    model.addConstr(vert_variables[vid] == 0);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }


        // conflict constraints
        std::set< std::pair<uint, uint> > vert_constraints;
        globalVertexConstraints(submesh, vert_constraints);


        for(auto &p : vert_constraints)
        {
            uint v0_id = p.first;
            uint v1_id = p.second;

            model.addConstr(vert_variables[v0_id] + vert_variables[v1_id] <= 1);
        }



        model.optimize();



        //Parse solution
        std::vector<uint> v_refs(submesh.num_verts(),0);
        std::vector<int> h_ref = submesh.vector_poly_labels();

        for(uint vid=0; vid<submesh.num_verts();vid++){
            if(submesh.vert_valence(vid) == 0) continue;

            double var = std::round(vert_variables[vid].get(GRB_DoubleAttr_X));
            int v_ref = static_cast<uint>(var);
            v_refs[vid] = v_ref;
            grid.vert_data(vid).label = v_ref;
            for(uint pid : submesh.adj_v2p(vid)){
                h_ref[pid] = std::max(h_ref[pid], v_ref);
            }
        }

        std::unordered_set<uint> refined;
        for(uint pid=0; pid<submesh.num_polys(); pid++){
            if(h_ref[pid] == 1 && submesh.poly_data(pid).label == 0){
                split_cell(grid, poly_s2g[pid]);
                refined.insert(pid);
            }
        }

        //Refine minor
        std::unordered_map<uint, std::vector<uint>> to_refine;
        for(uint vid=0; vid<submesh.num_verts(); vid++){
            if(v_refs[vid] == 0 || !submesh.vert_is_on_srf(vid)) continue;

            bool should_refine=false;
            for(uint pid : submesh.adj_v2p(vid)){
                if(refined.find(pid) != refined.end()){
                    should_refine = true;
                    break;
                }
            }

            if(should_refine){
                uint num_nodes_to_refine = 0;
                for(uint pid : grid.adj_v2p(vid)){
                    if(grid.poly_data(pid).is_leaf && grid.poly_data(pid).label < static_cast<int>(refinement)){
                        to_refine[pid].push_back(grid.poly_vert_offset(pid, vid));
                        grid.vert_data(vid).label = 1;
                        num_nodes_to_refine++;
                    }
                }

                if(num_nodes_to_refine == 0 && !vert_is_on_border(grid, vid)){ //Hanging vert
                    uint v1, v2;
                    int pid = find_adj_hanging_vert(grid, vid, v1,v2);
                    if(pid == -1){
                        std::cerr<<"CELL NOT FOUND"<<std::endl;
                        continue;
                    }
                    to_refine[pid].push_back(grid.poly_vert_offset(pid, v1));
                    to_refine[pid].push_back(grid.poly_vert_offset(pid, v2));
                }
            }
        }
        for(auto &p : to_refine){
            refine_minors(grid, p.first, p.second, refinement+1);
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
