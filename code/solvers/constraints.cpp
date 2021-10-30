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

#include "constraints.h"

using namespace cinolib;

static std::vector<vec3d> confl_verts = {vec3d(-1, -2, -2),vec3d(1, -2, -2), vec3d(-2, -1, -2), vec3d(-1, -1, -2), vec3d(0, -1, -2), vec3d(1, -1, -2), vec3d(2, -1, -2),  vec3d(-1, 0, -2),  vec3d(1, 0, -2),  vec3d(-2, 1, -2),
                                         vec3d(-1, 1, -2), vec3d(0, 1, -2),  vec3d(1, 1, -2),   vec3d(2, 1, -2),   vec3d(-1, 2, -2), vec3d(1, 2, -2),  vec3d(-2, -2, -1), vec3d(-1, -2, -1), vec3d(0, -2, -1), vec3d(1, -2, -1),
                                         vec3d(2, -2, -1), vec3d(-2, -1, -1),vec3d(-1, -1, -1), vec3d(0, -1, -1),  vec3d(1, -1, -1), vec3d(2, -1, -1), vec3d(-2, 0, -1),  vec3d(-1, 0, -1),  vec3d(0, 0, -1),  vec3d(1, 0, -1),
                                         vec3d(2, 0, -1),  vec3d(-2, 1, -1), vec3d(-1, 1, -1),  vec3d(0, 1, -1),   vec3d(1, 1, -1),  vec3d(2, 1, -1),  vec3d(-2, 2, -1),  vec3d(-1, 2, -1),  vec3d(0, 2, -1),  vec3d(1, 2, -1),
                                         vec3d(2, 2, -1),  vec3d(-1, -2, 0), vec3d(1, -2, 0),   vec3d(-2, -1, 0),  vec3d(-1, -1, 0), vec3d(0, -1, 0),  vec3d(1, -1, 0),   vec3d(2, -1, 0),   vec3d(-1, 0, 0),  vec3d(1, 0, 0),
                                         vec3d(-2, 1, 0),  vec3d(-1, 1, 0),  vec3d(0, 1, 0),    vec3d(1, 1, 0),    vec3d(2, 1, 0),   vec3d(-1, 2, 0),  vec3d(1, 2, 0),    vec3d(-2, -2, 1),  vec3d(-1, -2, 1), vec3d(0, -2, 1),
                                         vec3d(1, -2, 1),  vec3d(2, -2, 1),  vec3d(-2, -1, 1),  vec3d(-1, -1, 1),  vec3d(0, -1, 1),  vec3d(1, -1, 1),  vec3d(2, -1, 1),   vec3d(-2, 0, 1),   vec3d(-1, 0, 1),  vec3d(0, 0, 1),
                                         vec3d(1, 0, 1),   vec3d(2, 0, 1),   vec3d(-2, 1, 1),   vec3d(-1, 1, 1),   vec3d(0, 1, 1),   vec3d(1, 1, 1),   vec3d(2, 1, 1),    vec3d(-2, 2, 1),   vec3d(-1, 2, 1),  vec3d(0, 2, 1),
                                         vec3d(1, 2, 1),   vec3d(2, 2, 1),   vec3d(-1, -2, 2),  vec3d(1, -2, 2),   vec3d(-2, -1, 2), vec3d(-1, -1, 2), vec3d(0, -1, 2),   vec3d(1, -1, 2),   vec3d(2, -1, 2),  vec3d(-1, 0, 2),
                                         vec3d(1, 0, 2),   vec3d(-2, 1, 2),  vec3d(-1, 1, 2),   vec3d(0, 1, 2),    vec3d(1, 1, 2),   vec3d(2, 1, 2),   vec3d(-1, 2, 2),   vec3d(1, 2, 2),
                                         vec3d(-3, -3, -2),   vec3d(-3, -3, -1),   vec3d(-3, -3, 0),   vec3d(-3, -3, 1),   vec3d(-3, -3, 2),   vec3d(-3, -2, -3),   vec3d(-3, -2, -2),   vec3d(-3, -2, -1),   vec3d(-3, -2, 0),   vec3d(-3, -2, 1),
                                         vec3d(-3, -2, 2),   vec3d(-3, -2, 3),   vec3d(-3, -1, -3),   vec3d(-3, -1, -2),   vec3d(-3, -1, -1),   vec3d(-3, -1, 0),   vec3d(-3, -1, 1),   vec3d(-3, -1, 2),   vec3d(-3, -1, 3),   vec3d(-3, 0, -3),
                                         vec3d(-3, 0, -2),   vec3d(-3, 0, -1),   vec3d(-3, 0, 0),   vec3d(-3, 0, 1),   vec3d(-3, 0, 2),   vec3d(-3, 0, 3),   vec3d(-3, 1, -3),   vec3d(-3, 1, -2),   vec3d(-3, 1, -1),   vec3d(-3, 1, 0),
                                         vec3d(-3, 1, 1),   vec3d(-3, 1, 2),   vec3d(-3, 1, 3),   vec3d(-3, 2, -3),   vec3d(-3, 2, -2),   vec3d(-3, 2, -1),   vec3d(-3, 2, 0),   vec3d(-3, 2, 1),   vec3d(-3, 2, 2),   vec3d(-3, 2, 3),
                                         vec3d(-3, 3, -2),   vec3d(-3, 3, -1),   vec3d(-3, 3, 0),   vec3d(-3, 3, 1),   vec3d(-3, 3, 2),   vec3d(-2, -3, -3),   vec3d(-2, -3, -2),   vec3d(-2, -3, -1),   vec3d(-2, -3, 0),   vec3d(-2, -3, 1),
                                         vec3d(-2, -3, 2),   vec3d(-2, -3, 3),   vec3d(-2, -2, -3),   vec3d(-2, -2, 3),   vec3d(-2, -1, -3),   vec3d(-2, -1, 3),   vec3d(-2, 0, -3),   vec3d(-2, 0, 3),   vec3d(-2, 1, -3),   vec3d(-2, 1, 3),
                                         vec3d(-2, 2, -3),   vec3d(-2, 2, 3),   vec3d(-2, 3, -3),   vec3d(-2, 3, -2),   vec3d(-2, 3, -1),   vec3d(-2, 3, 0),   vec3d(-2, 3, 1),   vec3d(-2, 3, 2),   vec3d(-2, 3, 3),   vec3d(-1, -3, -3),
                                         vec3d(-1, -3, -2),   vec3d(-1, -3, -1),   vec3d(-1, -3, 0),   vec3d(-1, -3, 1),   vec3d(-1, -3, 2),   vec3d(-1, -3, 3),   vec3d(-1, -2, -3),   vec3d(-1, -2, 3),   vec3d(-1, -1, -3),
                                         vec3d(-1, -1, 3),   vec3d(-1, 0, -3),   vec3d(-1, 0, 3),   vec3d(-1, 1, -3),   vec3d(-1, 1, 3),   vec3d(-1, 2, -3),   vec3d(-1, 2, 3),   vec3d(-1, 3, -3),   vec3d(-1, 3, -2),   vec3d(-1, 3, -1),
                                         vec3d(-1, 3, 0),   vec3d(-1, 3, 1),   vec3d(-1, 3, 2),   vec3d(-1, 3, 3),   vec3d(0, -3, -3),   vec3d(0, -3, -2),   vec3d(0, -3, -1),   vec3d(0, -3, 0),   vec3d(0, -3, 1),   vec3d(0, -3, 2),
                                         vec3d(0, -3, 3),   vec3d(0, -2, -3),   vec3d(0, -2, 3),   vec3d(0, -1, -3),   vec3d(0, -1, 3),   vec3d(0, 0, -3),   vec3d(0, 0, 3),   vec3d(0, 1, -3),   vec3d(0, 1, 3),   vec3d(0, 2, -3),
                                         vec3d(0, 2, 3),   vec3d(0, 3, -3),   vec3d(0, 3, -2),   vec3d(0, 3, -1),   vec3d(0, 3, 0),   vec3d(0, 3, 1),   vec3d(0, 3, 2),   vec3d(0, 3, 3),   vec3d(1, -3, -3),   vec3d(1, -3, -2),
                                         vec3d(1, -3, -1),   vec3d(1, -3, 0),   vec3d(1, -3, 1),   vec3d(1, -3, 2),   vec3d(1, -3, 3),   vec3d(1, -2, -3),   vec3d(1, -2, 3),   vec3d(1, -1, -3),   vec3d(1, -1, 3),   vec3d(1, 0, -3),
                                         vec3d(1, 0, 3),   vec3d(1, 1, -3),   vec3d(1, 1, 3),   vec3d(1, 2, -3),   vec3d(1, 2, 3),   vec3d(1, 3, -3),   vec3d(1, 3, -2),   vec3d(1, 3, -1),   vec3d(1, 3, 0),   vec3d(1, 3, 1),
                                         vec3d(1, 3, 2),   vec3d(1, 3, 3),   vec3d(2, -3, -3),   vec3d(2, -3, -2),   vec3d(2, -3, -1),   vec3d(2, -3, 0),   vec3d(2, -3, 1),   vec3d(2, -3, 2),   vec3d(2, -3, 3),   vec3d(2, -2, -3),
                                         vec3d(2, -2, 3),   vec3d(2, -1, -3),   vec3d(2, -1, 3),   vec3d(2, 0, -3),   vec3d(2, 0, 3),   vec3d(2, 1, -3),   vec3d(2, 1, 3),   vec3d(2, 2, -3),   vec3d(2, 2, 3),   vec3d(2, 3, -3),
                                         vec3d(2, 3, -2),   vec3d(2, 3, -1),   vec3d(2, 3, 0),   vec3d(2, 3, 1),   vec3d(2, 3, 2),   vec3d(2, 3, 3),   vec3d(3, -3, -2),   vec3d(3, -3, -1),   vec3d(3, -3, 0),   vec3d(3, -3, 1),
                                         vec3d(3, -3, 2),   vec3d(3, -2, -3),   vec3d(3, -2, -2),   vec3d(3, -2, -1),   vec3d(3, -2, 0),   vec3d(3, -2, 1),   vec3d(3, -2, 2),   vec3d(3, -2, 3),   vec3d(3, -1, -3),   vec3d(3, -1, -2),
                                         vec3d(3, -1, -1),   vec3d(3, -1, 0),   vec3d(3, -1, 1),   vec3d(3, -1, 2),   vec3d(3, -1, 3),   vec3d(3, 0, -3),   vec3d(3, 0, -2),   vec3d(3, 0, -1),   vec3d(3, 0, 0),   vec3d(3, 0, 1),
                                         vec3d(3, 0, 2),   vec3d(3, 0, 3),   vec3d(3, 1, -3),   vec3d(3, 1, -2),   vec3d(3, 1, -1),   vec3d(3, 1, 0),   vec3d(3, 1, 1),   vec3d(3, 1, 2),   vec3d(3, 1, 3),   vec3d(3, 2, -3),
                                         vec3d(3, 2, -2),   vec3d(3, 2, -1),   vec3d(3, 2, 0),   vec3d(3, 2, 1),   vec3d(3, 2, 2),   vec3d(3, 2, 3),   vec3d(3, 3, -2),   vec3d(3, 3, -1),   vec3d(3, 3, 0),   vec3d(3, 3, 1),
                                         vec3d(3, 3, 2)};


static std::vector<uint> overlappd_cells = {22,23,24,27,28,29,32,33,34,44,45,46,48,49,51,52,53,63,64,65,68,69,70,73,74,75};
static std::vector<uint> overlappd_faces = {3,4,5,7,8,10,11,12,17,18,19,21,25,26,30,31,35,37,38,39,41,42,43,47,50,54,55,56,58,59,60,62,66,67,71,72,76,78,79,80,85,86,87,89,90,92,93,94};
static std::vector<uint> overlappd_edges = {0,20,40,1,61,81,91,2,82,13,83,14,84,15,95,6,16,36,96,57,77,97,88,9};



void vertexConstraints(const cinolib::Hexmesh<> &m, uint start_v, double edge_l, std::map<cinolib::vec3d, uint, vert_comparator> &vmap, std::vector<uint> &constr, bool verbose)
{

    vec3d ref_v = m.vert(start_v);

    for(uint v_id = 0; v_id < confl_verts.size(); v_id++)
    {
        vec3d query_v = (confl_verts[v_id] * edge_l) + ref_v;

        uint count3=0;
        uint count0=0;
        for(uint i=0;i<3;i++){
            if(abs(confl_verts[v_id][i]) == 3) count3++;
            if(abs(confl_verts[v_id][i]) == 0 || abs(confl_verts[v_id][i]) == 2) count0++;
        }
        //if(count3 == 2 && count0 == 1) continue;

        auto find = vmap.find(query_v);

        if(find != vmap.end())
            constr.push_back(find->second);
    }

    if(verbose){
        if(m.vert_is_on_srf(start_v)){
            std::cout<<"ON SURFACE: ";
        }
        std::cout << "Num TOT constr fo vid "<<start_v<<": "<< constr.size() << std::endl;
    }

}

void globalVertexConstraints(const cinolib::Hexmesh<> &m, std::set< std::pair<uint, uint> > &global_constr)
{
    std::map<cinolib::vec3d, uint,vert_comparator> v_map;
    for(uint vid=0; vid<m.num_verts();vid++)
        v_map[m.vert(vid)] = vid;

    double edge_length = m.edge_length(0);
    uint max = 0;
    for(uint v_id = 0; v_id < m.num_verts(); v_id++)
    {
        if(m.vert_valence(v_id)==0) continue;

        std::vector<uint> v_constr;
        vertexConstraints(m, v_id, edge_length, v_map, v_constr, false);
        max = std::max(max, static_cast<uint>(v_constr.size()));

        for(uint vc : v_constr)
        {
            assert(vc != v_id);
            if(v_id < vc)
                global_constr.insert(std::make_pair(v_id, vc));
            else
                global_constr.insert(std::make_pair(vc, v_id));

        }
    }
    //std::cerr<<"N constr: "<<max<<std::endl;

}
