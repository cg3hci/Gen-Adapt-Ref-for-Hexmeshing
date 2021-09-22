#ifndef SOLVER_ADAPTIVE_GRID_H
#define SOLVER_ADAPTIVE_GRID_H

#include <gurobi_c++.h>
#include "../adaptive_grid.h"
#include <cinolib/dijkstra.h>
#include "constraints.h"


void solve(AdaptiveGrid &grid, bool use_weak_balancing=true);
void perform_pairing(AdaptiveGrid &grid, const cinolib::Hexmesh<> &submesh, uint refinement, std::unordered_map<uint, uint> &poly_g2s, std::unordered_map<uint, uint> &poly_s2g);

#include "solver_adaptive_grid.cpp"

#endif // SOLVER_ADAPTIVE_GRID_H
