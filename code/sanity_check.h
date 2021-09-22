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

#ifndef SANITY_CHECK_H
#define SANITY_CHECK_H

#include <cinolib/meshes/meshes.h>
#include "utils.h"

inline int normal2int(const cinolib::vec3d &n)
{
    double eps = 0.0001;
    if(std::abs(n.y()) < eps && std::abs(n.z()) < eps) return 0;
    if(std::abs(n.x()) < eps && std::abs(n.z()) < eps) return 1;
    return 2;
}

inline int nextVert(cinolib::Quadmesh<> &m, uint start_v)
{
    if(m.vert_data(start_v).label == 1) return -1; //corner

    for(auto &e : m.adj_v2e(start_v))
        if(m.edge_data(e).label == 1)
        {
            m.edge_data(e).label = 2; //visited

            if(m.edge_vert_id(e, 0) == start_v)
                 return static_cast<int>(m.edge_vert_id(e, 1));
            else return static_cast<int>(m.edge_vert_id(e, 0));
        }

    return -1;
}

bool sanityCheckSingleMesh(cinolib::Quadmesh<> &m);
bool sanityCheck(std::vector<std::pair<cinolib::Quadmesh<>, uint>> &ref_meshes);
bool sanity_check(const cinolib::Hexmesh<> &m);

#include "sanity_check.cpp"

#endif // SANITY_CHECK_H
