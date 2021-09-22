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

#include "cino_cgal_interface.h"

template <class HDS>
class BuildMesh : public CGAL::Modifier_base<HDS> {
public:
    cinolib::Trimesh<> mesh;
    BuildMesh(const cinolib::Trimesh<> &mesh) {this->mesh=mesh;}
    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.

        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( mesh.num_verts(), mesh.num_polys(), 0);
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        for(const cinolib::vec3d &vert : mesh.vector_verts()){
            B.add_vertex( Point( vert.x(), vert.y(), vert.z()));
        }
        for(const std::vector<uint> &poly: mesh.vector_polys()){
            B.begin_facet();
            for(uint pid : poly){
                B.add_vertex_to_facet(pid);
            }
            B.end_facet();
        }
        B.end_surface();
    }
};
void get_sdf_values(const cinolib::Trimesh<> &mesh, std::vector<double> &sdf)
{
    Polyhedron P;
    BuildMesh<HalfedgeDS> cgal_mesh(mesh);
    P.delegate(cgal_mesh);

    typedef std::map<face_descriptor, double> Face_double_map;
    Face_double_map internal_map;
    boost::associative_property_map<Face_double_map> sdf_property_map(internal_map);
    // compute SDF values
    //std::pair<double, double> min_max_sdf = CGAL::sdf_values(P, sdf_property_map);
    // It is possible to compute the raw SDF values and post-process them using
    // the following lines:
     const std::size_t number_of_rays = 25;  // cast 25 rays per face
     const double cone_angle = 2.0 / 3.0 * CGAL_PI; // set cone opening-angle
     std::pair<double, double> min_max_sdf = CGAL::sdf_values(P, sdf_property_map, cone_angle, number_of_rays, false);
    // std::pair<double, double> min_max_sdf =
    //  CGAL::sdf_values_postprocessing(mesh, sdf_property_map);
    // print minimum & maximum SDF values
    std::cout << "minimum SDF: " << min_max_sdf.first
              << " maximum SDF: " << min_max_sdf.second << std::endl;
    // print SDF values
    for(face_descriptor f : faces(P)) {
        sdf.push_back(sdf_property_map[f]);
    }
}

