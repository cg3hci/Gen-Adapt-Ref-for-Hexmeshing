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

#include "utils.h"

struct Vert{
    cinolib::vec3d coords;
    uint idx;
};

struct order_x
{
    inline bool operator() (const Vert& a, const Vert& b)
    {
        return (a.coords.x() < b.coords.x());
    }
};
struct order_y
{
    inline bool operator() (const Vert& a, const Vert& b)
    {
        return (a.coords.y() < b.coords.y());
    }
};
struct order_z
{
    inline bool operator() (const Vert& a, const Vert& b)
    {
        return (a.coords.z() < b.coords.z());
    }
};

void poly_vert_ordering(const std::vector<cinolib::vec3d> &vertices, std::vector<uint> &poly){

    std::vector<Vert> tmp(8);
    for(uint i=0; i<8; i++){
        tmp[i] = {vertices[poly[i]], poly[i]};
    }
    std::sort(tmp.begin(), tmp.end(), order_y());
    std::sort(tmp.begin(), tmp.begin()+4, order_x());
    std::sort(tmp.begin()+4, tmp.end(), order_x());
    std::sort(tmp.begin(), tmp.begin()+2, order_z());
    std::sort(tmp.begin()+2, tmp.begin()+4, order_z());
    std::sort(tmp.begin()+4, tmp.begin()+6, order_z());
    std::sort(tmp.begin()+6, tmp.end(), order_z());
    std::swap(tmp[0], tmp[3]);
    std::swap(tmp[0], tmp[1]);
    std::swap(tmp[4], tmp[7]);
    std::swap(tmp[4], tmp[5]);

    for(uint i=0; i<8; i++)
        poly[i] = tmp[i].idx;


}

void remove_unreferenced_vertices(cinolib::AbstractPolygonMesh<> &m){

    std::vector<uint> to_remove;
    for(uint vid=0; vid<m.num_verts(); vid++){
       if(m.vert_valence(vid) == 0){
           to_remove.push_back(vid);
       }
    }
    std::sort(to_remove.begin(), to_remove.end(), std::greater<uint>());
    for(uint vid : to_remove){
        m.vert_remove_unreferenced(vid);
    }
}

void remove_unreferenced_vertices(cinolib::AbstractPolyhedralMesh<> &m){

    std::vector<uint> to_remove;
    for(uint vid=0; vid<m.num_verts(); vid++){
       if(m.vert_valence(vid) == 0){
           to_remove.push_back(vid);
       }
    }
    std::sort(to_remove.begin(), to_remove.end(), std::greater<uint>());
    for(uint vid : to_remove){
        m.vert_remove_unreferenced(vid);
    }

}

std::vector<std::string> split_string(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string substring;
    while(std::getline(ss, substring, delim)) {
        elems.push_back(substring);
    }
    return elems;
}