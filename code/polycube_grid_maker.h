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

#ifndef POLYCUBE_GRID_MAKER_H
#define POLYCUBE_GRID_MAKER_H

#include "adaptive_grid.h"
#include "split_metrics.h"
#include "sanity_check.h"
#include "solvers/solver_adaptive_grid.h"


void solve_grid_polycube(AdaptiveGrid &grid, bool octree, bool weak_balancing);

/**
*This function performs all the steps necessary to build a balanced and paired adaptive grid 
*from an input tetrahedral mesh and its polycube map.
*
*@param m: the input mesh
*@param pc_m: the input polycube mesh
*@param grid: the output balanced and paired grid
*@param min_refinement: minimum depth of the tree representing the grid 
*@param max_refinement: maximum depth of the tree representing the grid 
*@param weak_balancing: if false use strong balancing instead of weak balancing
*@param split_metric: function that defines the split policy of a grid cell
*@param perform_sanity_check: if false the correcteness of the pairing is not checked at the end of the process
*/
void make_grid(const cinolib::Tetmesh<> &m, const cinolib::Tetmesh<> &pc_m, cinolib::Hexmesh<> &grid, 
               uint min_refinement, uint max_refinement, 
               bool weak_balancing, bool use_octree,  
               bool(*split_metric)(const AdaptiveGrid &, uint, const cinolib::Tetmesh<> &, const cinolib::Tetmesh<> &) = polycube_distortion_split_metric, 
               bool perform_sanity_check = true);

#include "polycube_grid_maker.cpp"

#endif
