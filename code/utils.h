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

#ifndef UTILS_H
#define UTILS_H

#include <cinolib/meshes/meshes.h>

struct vert_compare{
    bool operator()(const cinolib::vec3d& a, const cinolib::vec3d& b) const {

       double eps = 1e-6;
       if(a.x()-b.x() < 0.0 && abs(a.x()-b.x()) > eps){
           return true;
       }
       else if(abs(a.x()-b.x()) < eps){
           if(a.y()-b.y() < 0.0 && abs(a.y()-b.y()) > eps){
               return true;
           }
           else if(abs(a.y()-b.y()) < eps){
               if(a.z()-b.z() < 0.0 && abs(a.z()-b.z()) > eps){
                   return true;
               }
           }
       }

       return false;
    }
};

inline bool eps_equal(const double &a, const double &b)
{
    return std::abs(a - b) < 1e-5;
}

inline std::chrono::time_point<std::chrono::system_clock> startChrono()
{
    return std::chrono::system_clock::now();
}

inline double stopChrono(std::chrono::time_point<std::chrono::system_clock> &start)
{
    auto time = std::chrono::system_clock::now() - start;
    return std::chrono::duration <double, std::milli> (time).count() / 1000;
}

void remove_unreferenced_vertices(cinolib::AbstractPolyhedralMesh<> &m);
void remove_unreferenced_vertices(cinolib::AbstractPolygonMesh<> &m);
void poly_vert_ordering(const std::vector<cinolib::vec3d> &vertices, std::vector<uint> &poly);
std::vector<std::string> split_string(const std::string &s, char delim);

#include "utils.cpp"

#endif // UTILS_H