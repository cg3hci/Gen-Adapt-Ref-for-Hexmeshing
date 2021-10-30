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

#include "install.h"

inline bool polys_share_edge(const cinolib::AbstractPolyhedralMesh<> &m, const std::vector<uint> &polys)
{
    for(uint i = 0; i < polys.size(); i++)
    {
        for(uint j = 0; j < polys.size(); j++)
        {
            auto verts_h1 = m.poly_verts_id(polys[i]);
            auto verts_h2 = m.poly_verts_id(polys[j]);

            std::set<uint> s1(verts_h1.begin(), verts_h1.end());
            std::set<uint> s2(verts_h2.begin(), verts_h2.end());

            std::vector<uint> inters;
            std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), std::back_inserter(inters));

            if(inters.size() < 2)
                return false;
        }
    }
    return true;
}

void mark_candidates(cinolib::AbstractPolyhedralMesh<> &hm){

    for(uint vid=0; vid<hm.num_verts(); vid++){

        std::vector<uint> refinements;
        std::set<uint> refs_set;

        for(uint pid : hm.adj_v2p(vid)){
            refs_set.insert(hm.poly_data(pid).label);
            refinements.push_back(hm.poly_data(pid).label);
        }

        if(refs_set.size() == 2 && hm.adj_v2p(vid).size()==8){ //INSIDE
            std::sort(refinements.begin(), refinements.end());
            if(refinements[0] == refinements[3] && refinements[4] == refinements[7]){

                int ref = *refs_set.begin();


                std::vector<uint> cluster1;
                std::vector<uint> cluster2;

                for(uint pid : hm.adj_v2p(vid)){
                    if(hm.poly_data(pid).label == ref){
                        cluster1.push_back(pid);
                    }
                    else{
                        cluster2.push_back(pid);
                    }
                }

                if(polys_share_edge(hm, cluster1) && polys_share_edge(hm, cluster2))
                    hm.vert_data(vid).label = 1;
            }
        }
        else if(refs_set.size() == 2 && hm.adj_v2p(vid).size()==4 && hm.vert_is_on_srf(vid)){ //SURFACE
            std::sort(refinements.begin(), refinements.end());
            if(refinements[0] == refinements[1] && refinements[2] == refinements[3]){
                int ref = *refs_set.begin();


                std::vector<uint> cluster1;
                std::vector<uint> cluster2;

                for(uint pid : hm.adj_v2p(vid)){
                    if(hm.poly_data(pid).label == ref){
                        cluster1.push_back(pid);
                    }
                    else{
                        cluster2.push_back(pid);
                    }
                }

                if(polys_share_edge(hm, cluster1) && polys_share_edge(hm, cluster2))
                    hm.vert_data(vid).label = 1;
            }
        }
        else if(refs_set.size() == 2 && hm.adj_v2p(vid).size()==2 && hm.vert_is_on_srf(vid)){ //EDGE

            hm.vert_data(vid).label = 1;

        }


    }
}

bool install_schemes_to_grid(cinolib::Polyhedralmesh<> &input_grid, cinolib::Hexmesh<> &output_grid, bool clip_cells)
{
    mark_candidates(input_grid);
    bool passed=true;
    cinolib::Polyhedralmesh<> out;

    for(uint i : {4,3,2,5}){
        
        find_t_vertices(input_grid, i);

        std::vector<bool> transition_verts(input_grid.num_verts());
        for(uint vid=0; vid<input_grid.num_verts(); vid++){
            if(input_grid.vert_data(vid).label == 2){
                transition_verts[vid] = true;
            }
        }

        cinolib::hex_transition_install(input_grid, transition_verts, out);

        cinolib::Polygonmesh<> peel;
        cinolib::export_surface(out, peel);

        if(connected_components(peel) != 1){
            std::cerr<<"Something went wrong, retrying with different paramenters"<<std::endl;
            passed = false;
        }
        else{
            passed = true;
            break;
        }

    }

    if(passed) std::cout<<"Schemes installed successfully :)"<<std::endl;
    else{
        std::cerr<<"Could not install schemes on this grid"<<std::endl;
        return false;
    }

    std::cout<<"Making dual mesh..."<<std::endl;
    //DUAL::::::::::::::::::::::::::::::::::::::::::::::::::
    out.edge_mark_sharp_creases();
    cinolib::Polyhedralmesh<> dual;
    cinolib::dual_mesh(out, dual, clip_cells);
    std::vector<std::vector<uint>> polys(dual.num_polys());
    for(uint pid=0; pid<dual.num_polys();pid++){
        polys[pid] = dual.poly_verts_id(pid);
    }

    cinolib::Hexmesh<> dual_full_hexa(dual.vector_verts(), polys);

    for(uint pid=0; pid<dual_full_hexa.num_polys(); pid++){
        auto poly_verts = dual_full_hexa.poly_verts(pid);

        if(cinolib::hex_volume(poly_verts[0],poly_verts[1],poly_verts[2],poly_verts[3],poly_verts[4],poly_verts[5],poly_verts[6],poly_verts[7]) < 0){
            dual_full_hexa.poly_flip_winding(pid);
            dual_full_hexa.poly_reorder_p2v(pid);
        }
    }

    output_grid = dual_full_hexa;


    return true;
}
