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

#include <stdlib.h>
#include <iostream>
#include <regex>
#include "surface_grid_maker.h"
#include "polycube_grid_maker.h"
#include "hex_transitions_install/install.h"
#include "mesh_projection/project.h"

int main(int argc, char *argv[]){

    std::regex regex_sw_balancing("--strong_balancing|--weak_balancing");
    std::regex regex_use_octree("--use_octree");
    std::regex regex_max_refinement("--max_refinement=[0-9]+");
    std::regex regex_mesh_path("--input_mesh_path=.*");
    std::regex regex_pc_mesh_path("--input_pc_mesh_path=.*");
    std::regex regex_grid_path("--output_grid_path=.*.mesh");
    std::regex regex_sanity_check("--sanity_check=(true|false)");
    std::regex regex_install_schemes("--install_schemes=(true|false)");
    std::regex regex_project_mesh("--project_mesh=(true|false)");

    if(argc > 3){

        if(std::string(argv[1]) == "--surface"){

            std::string mesh_path;
            std::string grid_save_path;
            uint min_refinement = 0;
            uint max_refinement = 8;
            bool weak_balancing = true;
            bool use_octree = false;
            bool perform_sanity_check = true;
            bool install_schemes = false;
            bool project_mesh_ = false;

            bool mesh_found, grid_found;

            for(int i=2;i<argc;i++){

                std::string argument(argv[i]);

                if(std::regex_match(argument, regex_sw_balancing)){
                    if(argument == "--strong_balancing") weak_balancing = false;
                }
                else if(std::regex_match(argument, regex_use_octree)) use_octree = true;
                else if(std::regex_match(argument, regex_max_refinement)){
                    max_refinement = std::stoi(argument.erase(0, argument.find("=") + 1));
                }
                else if(std::regex_match(argument, regex_mesh_path)){
                    mesh_path = argument.erase(0, argument.find("=") + 1);
                    mesh_found = true;
                }
                else if(std::regex_match(argument, regex_grid_path)){
                    grid_save_path = argument.erase(0, argument.find("=") + 1);
                    grid_found = true;
                }
                else if(std::regex_match(argument, regex_sanity_check)){
                    std::string bool_string = argument.erase(0, argument.find("=") + 1);
                    perform_sanity_check = bool_string.front() == 't' ? true : false;
                }
                else if(std::regex_match(argument, regex_install_schemes)){
                    std::string bool_string = argument.erase(0, argument.find("=") + 1);
                    install_schemes = bool_string.front() == 't' ? true : false;
                }
                else if(std::regex_match(argument, regex_project_mesh)){
                    std::string bool_string = argument.erase(0, argument.find("=") + 1);
                    project_mesh_ = bool_string.front() == 't' ? true : false;
                }
            }

            if(!(mesh_found && grid_found)){
                std::cerr<<"Required argument --input_mesh_path or --output_grid_path not found"<<std::endl;
                return -1;
            }

            cinolib::Trimesh<> m(mesh_path.c_str());
            cinolib::Hexmesh<> grid;

            make_grid(m, grid, min_refinement, max_refinement, weak_balancing, use_octree, sdf_split_metric, perform_sanity_check);
            grid.save(grid_save_path.c_str());

            if(install_schemes || project_mesh_){
                cinolib::Polyhedralmesh<> pm(grid_save_path.c_str());
                cinolib::Hexmesh<> conforming_grid;
                bool result = install_schemes_to_grid(pm, conforming_grid, true);
                std::vector<std::string> split_extension= split_string(grid_save_path, '.');
                if(result) conforming_grid.save((split_extension[0]+"_conforming.mesh").c_str());

                if(project_mesh_){
                    project_mesh(conforming_grid, m);
                    conforming_grid.save((split_extension[0]+"_projected.mesh").c_str()); 
                }
            }
            
        }


        if(std::string(argv[1]) == "--polycube"){

            std::string mesh_path;
            std::string pc_mesh_path;
            std::string grid_save_path;
            uint min_refinement = 5;
            uint max_refinement = 8;
            bool weak_balancing = true;
            bool use_octree = false;
            bool perform_sanity_check = true;
            bool install_schemes = false;
            bool project_mesh_ = false;
            bool mesh_found, pc_mesh_found, grid_found;

            for(int i=2;i<argc;i++){

                std::string argument(argv[i]);

                if(std::regex_match(argument, regex_sw_balancing)){
                    if(argument == "--strong_balancing") weak_balancing = false;
                }
                else if(std::regex_match(argument, regex_use_octree)) use_octree = true;
                else if(std::regex_match(argument, regex_max_refinement)){
                    max_refinement = std::stoi(argument.erase(0, argument.find("=") + 1));
                }
                else if(std::regex_match(argument, regex_mesh_path)){
                    mesh_path = argument.erase(0, argument.find("=") + 1);
                    mesh_found = true;
                }
                else if(std::regex_match(argument, regex_pc_mesh_path)){
                    pc_mesh_path = argument.erase(0, argument.find("=") + 1);
                    pc_mesh_found = true;
                }
                else if(std::regex_match(argument, regex_grid_path)){
                    grid_save_path = argument.erase(0, argument.find("=") + 1);
                    grid_found = true;
                }
                else if(std::regex_match(argument, regex_sanity_check)){
                    std::string bool_string = argument.erase(0, argument.find("=") + 1);
                    perform_sanity_check = bool_string.front() == 't' ? true : false;
                }
                else if(std::regex_match(argument, regex_install_schemes)){
                    std::string bool_string = argument.erase(0, argument.find("=") + 1);
                    install_schemes = bool_string.front() == 't' ? true : false;
                }
                else if(std::regex_match(argument, regex_project_mesh)){
                    std::string bool_string = argument.erase(0, argument.find("=") + 1);
                    project_mesh_ = bool_string.front() == 't' ? true : false;
                }
            }

            if(!(mesh_found && grid_found && pc_mesh_found)){
                std::cerr<<"Required argument --input_mesh_path or --input_pc_mesh_path or --output_grid_path not found"<<std::endl;
                return -1;
            }

            cinolib::Tetmesh<> m(mesh_path.c_str());
            cinolib::Tetmesh<> pc_m(pc_mesh_path.c_str());
            cinolib::Hexmesh<> grid;

            make_grid(m, pc_m, grid, min_refinement, max_refinement, weak_balancing, false, polycube_distortion_split_metric, perform_sanity_check);
            grid.save(grid_save_path.c_str());

            if(install_schemes || project_mesh_){
                cinolib::Polyhedralmesh<> pm(grid_save_path.c_str());
                cinolib::Hexmesh<> conforming_grid;
                bool result = install_schemes_to_grid(pm, conforming_grid, false);
                std::vector<std::string> split_extension= split_string(grid_save_path, '.');
                if(result) conforming_grid.save((split_extension[0]+"_conforming.mesh").c_str());

                if(project_mesh_){
                    project_mesh(conforming_grid, m, pc_m);
                    conforming_grid.save((split_extension[0]+"_projected.mesh").c_str()); 
                }
            }

        }


    }

    else if(argc == 2 && std::string(argv[1]) == "--help"){
        std::cout<<"usage: ./grid_maker (--surface | --polycube) --input_mesh_path=MESH_PATH --output_grid_path=GRID_PATH [Options]"<<std::endl;
        std::cout<<"Options:"<<std::endl;
        std::cout<<"--input_pc_mesh_path=PATH (required for polycube pipeline). Specify the path of the polycube map"<<std::endl;
        std::cout<<"--min_refinement=VALUE (optional, default 0[5 for polycube])"<<std::endl;
        std::cout<<"--max_refinement=VALUE (optional, default 8)"<<std::endl;
        std::cout<<"--use_octree (optional). Use the algorithmic pairing process (for surface pipeline only)"<<std::endl;
        std::cout<<"--weak_balancing | --strong_balancing (optional, default weak_balancing)"<<std::endl;
        std::cout<<"--sanity_check=BOOL (optional, default true). Test if the final mesh is paired correctly"<<std::endl;
        std::cout<<"--install_schemes=BOOL (optional, default false). Install the transition schemes to get a conforming all-hexa grid"<<std::endl;
        std::cout<<"--project_mesh=BOOL (optional, default false). Project the grid on the target mesh"<<std::endl;

    }
    else{
        std::cerr<<"Wrong number of parameter. Run with --help for more information";
    }

    return 0;
}
