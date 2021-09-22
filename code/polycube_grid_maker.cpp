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

#include "polycube_grid_maker.h"

void solve_grid_polycube(AdaptiveGrid &grid, bool octree, bool weak_balancing)
{
    if(octree)
    {
        if(weak_balancing){
            while(!(test_balancing(grid) && test_pairing(grid))){
                balancing(grid);
                pairing(grid);
            }
        }
        else{
            while(!(test_strong_balancing(grid) && test_pairing(grid))){
                strong_balancing(grid);
                pairing(grid);
            }
        }

    }
    else{
        solve(grid, weak_balancing);
    }

}

void mark_external_polys(const AdaptiveGrid &grid, std::unordered_map<uint, uint> &poly_map, std::vector<uint> &external_polys){
    
    external_polys.clear();
    std::vector<uint> externals_tmp;

    auto min_max = get_min_and_max_refinement(grid);
    for(uint pid=0; pid<grid.num_polys(); pid++){
        if(!grid.poly_data(pid).is_inside_polycube && grid.poly_data(pid).label==min_max.first){
            if(grid.poly_data(pid).is_leaf)
                externals_tmp.push_back(pid);
            std::deque<uint> q;
            q.push_back(pid);

            while(!q.empty()){
                auto curr = q.front();
                q.pop_front();
                for(auto child : grid.poly_data(curr).children){
                    if(!grid.poly_data(child).is_leaf){
                        q.push_back(child);
                    }
                    else{
                        externals_tmp.push_back(child);
                    }
                }
            }
        }
    }

    std::sort(externals_tmp.begin(), externals_tmp.end(), std::greater());
    for(uint pid : externals_tmp){
        external_polys.push_back(poly_map[pid]);
    }

    
}



void make_grid(const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m, cinolib::Hexmesh<> &grid, uint min_refinement, uint max_refinement, bool weak_balancing, bool use_octree, bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Tetmesh<> &, const cinolib::Tetmesh<> &), bool perform_sanity_check){

    std::string external_polys_path = grid.mesh_data().filename;

    AdaptiveGrid adaptive_grid;
    build_polycube(adaptive_grid,m,pc_m, split_metric, min_refinement, max_refinement);
    solve_grid(adaptive_grid, use_octree, weak_balancing);
    std::unordered_map<uint, uint> poly_map;
    extract_non_conformal_hexmesh(adaptive_grid, grid, poly_map);
    if(perform_sanity_check) sanity_check(grid);
    std::vector<uint> external_polys;

    mark_external_polys(adaptive_grid, poly_map, external_polys);
    std::ofstream f;
    f.open(external_polys_path);
    assert(f.is_open());
    for(uint pid : external_polys)
    {
        f <<pid<< "\n";
    }
    f.close();
}
