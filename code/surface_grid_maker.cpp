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

#include "surface_grid_maker.h"


void solve_grid(AdaptiveGrid &grid, bool octree, bool weak_balancing)
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


void make_grid(const cinolib::Trimesh<> &m, cinolib::Hexmesh<> &grid, uint min_refinement, uint max_refinement, bool weak_balancing, bool use_octree, bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Trimesh<> &), bool perform_sanity_check){

    AdaptiveGrid adaptive_grid;
    build(adaptive_grid, m, split_metric, min_refinement, max_refinement);
    solve_grid(adaptive_grid, use_octree, weak_balancing);
    std::unordered_map<uint, uint> poly_map;
    extract_non_conformal_hexmesh(adaptive_grid, grid, poly_map);
    if(perform_sanity_check) sanity_check(grid);

}
