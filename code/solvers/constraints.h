#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <cinolib/meshes/meshes.h>
#include "utils.h"

void vertexConstraints(const cinolib::Hexmesh<> &m, uint start_v, double edge_l, std::map<cinolib::vec3d, uint, vert_compare> &vmap, std::vector<uint> &constr, bool verbose);
void globalVertexConstraints(const cinolib::Hexmesh<> &m, std::set< std::pair<uint, uint> > &global_constr);

#include "constraints.cpp"

#endif // CONSTRAINTS_H
