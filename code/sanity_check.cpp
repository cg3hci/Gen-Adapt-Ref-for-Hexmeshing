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

#include "sanity_check.h"

bool sanityCheckSingleMesh(cinolib::Quadmesh<> &m)
{
    // label == 0 ignore, 1 valid for check, 2 already visited

    std::queue<uint> corners;

    // flag border edges
    for(uint e_id = 0; e_id < m.num_edges(); e_id++)
    {
        if(m.edge_is_boundary(e_id)) continue; // skyp boundary edges

        std::set<int> normals;

        for(auto &f : m.adj_e2p(e_id))
            normals.insert(normal2int(m.poly_data(f).normal));

        if(normals.size() == 2)
            m.edge_data(e_id).label = 1; // corner
    }

    // flag corner vertices
    for(uint v_id = 0; v_id < m.num_verts(); v_id++)
    {
        if(m.adj_v2p(v_id).size() == 4 || m.adj_v2p(v_id).size() == 8) continue;
        std::set<int> normals;

        for(auto &e : m.adj_v2e(v_id))
        {
            cinolib::vec3d n = (m.edge_vert(e, 0) - m.edge_vert(e, 1));
            n.normalize();
            normals.insert(normal2int(n));
        }

        if(normals.size() == 3 && !m.vert_is_boundary(v_id))
        {
            m.vert_data(v_id).label = 1; // corner
            corners.push(v_id);
        }
    }

    // sanity check
    while(!corners.empty())
    {
        uint start_v = corners.front();
        corners.pop();

        for(auto &e_id : m.adj_v2e(start_v))
        {
            if(m.edge_data(e_id).label == 1) // valid border edge
            {
                uint step = 0;
                int next_v = m.edge_vert_id(e_id, 0) == start_v ? static_cast<int>(m.edge_vert_id(e_id, 1)) : static_cast<int>(m.edge_vert_id(e_id, 0));
                m.edge_data(e_id).label = 2;
                step++;
                int prev_v;
                do
                {
                    prev_v = next_v;
                    next_v = nextVert(m, static_cast<uint>(next_v));
                    if(next_v != -1) step++;

                } while(next_v != -1); //corner

                if(step % 2 > 0 && m.vert_data(prev_v).label == 1)
                {
                    m.edge_data(e_id).label = 3;
                    return false;
                }
            }
        }
    }
    return true;
}


bool sanityCheck(std::vector<std::pair<cinolib::Quadmesh<>, uint>> &ref_meshes)
{
    for(uint r = 0; r < ref_meshes.size(); r++)
    {
        if(!sanityCheckSingleMesh(ref_meshes[r].first))
        {
            std::cerr << "sanity check FAILS at ref " << ref_meshes[r].second << std::endl;
            return false;
        }
    }

    std::cerr << "sanity check PASSED :)" << std::endl;
    return true;
}

void extract_quadmesh_for_refinement(const cinolib::Hexmesh<> &m, uint ref, cinolib::Quadmesh<> &qm){

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
    remove_unreferenced_vertices(qm);
    //std::string path = "/Users/lucapitzalis/Desktop/ref_"+std::to_string(ref)+".obj";
    //qm.save(path.c_str());

}

bool sanity_check(const cinolib::Hexmesh<> &m){

    const std::vector<int> &ref = m.vector_poly_labels();
    int max_ref = *std::max_element(ref.begin(), ref.end());
    int min_ref = *std::min_element(ref.begin(), ref.end());

    std::cout<<max_ref<<" "<<min_ref<<std::endl;
    std::vector<uint> failures;
    bool all_passed = true;
    for(int ref = max_ref-1; ref>=min_ref; --ref){
        std::vector<std::pair<cinolib::Quadmesh<>, uint>> to_check;
        cinolib::Quadmesh<> qm;
        extract_quadmesh_for_refinement(m, ref, qm);

        to_check.push_back(std::make_pair(qm, ref));
        bool passed = sanityCheck(to_check);
        all_passed = all_passed && passed;
    }

    return all_passed;

}
