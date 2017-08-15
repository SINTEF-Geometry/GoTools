/*
* Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
* Applied Mathematics, Norway.
*
* Contact information: E-mail: tor.dokken@sintef.no                      
* SINTEF ICT, Department of Applied Mathematics,                         
* P.O. Box 124 Blindern,                                                 
* 0314 Oslo, Norway.                                                     
*
* This file is part of GoTools.
*
* GoTools is free software: you can redistribute it and/or modify
* it under the terms of the GNU Affero General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version. 
*
* GoTools is distributed in the hope that it will be useful,        
* but WITHOUT ANY WARRANTY; without even the implied warranty of         
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
* GNU Affero General Public License for more details.
*
* You should have received a copy of the GNU Affero General Public
* License along with GoTools. If not, see
* <http://www.gnu.org/licenses/>.
*
* In accordance with Section 7(b) of the GNU Affero General Public
* License, a covered work must retain the producer line in every data
* file that is created or manipulated using GoTools.
*
* Other Usage
* You can be released from the requirements of the license by purchasing
* a commercial license. Buying such a license is mandatory as soon as you
* develop commercial activities involving the GoTools library without
* disclosing the source code of your own applications.
*
* This file may be used in accordance with the terms contained in a
* written agreement between you and SINTEF ICT. 
*/

#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/LineCloud.h"

#include <vector>
#include <iostream>
#include <fstream>
#include "json/json.h"
#include <string>

using std::vector;
using std::string;


// The indices denote the triangles, each triplet constituting a triangle.
void readData(std::istream& filein_json,
	      vector<double>& vertices,
	      vector<double>& normals,
	      vector<int>& indices)
{
    Json::Value root;
    Json::Reader reader;
//    std::string json_string(filein_json);
//    std::istream* filein_json_is = static_cast<std::ifstream*>(&filein_json);
    bool parsedSuccess = reader.parse(filein_json,
				      root,
				      false);

    std::cout << "root.size(): " << root.size() << std::endl;

    Json::Value json_vert = root.get("vertices", "ASCII");
    std::cout << "vert_size: " << json_vert.size() << std::endl;
    vertices.resize(json_vert.size());
    for (int ki = 0; ki < vertices.size(); ++ki)
    {
	vertices[ki] = json_vert[ki].asDouble();
    }

    Json::Value json_norm = root.get("normals", "ASCII");
    std::cout << "norm_size: " << json_norm.size() << std::endl;
    normals.resize(json_norm.size());
    for (int ki = 0; ki < normals.size(); ++ki)
    {
	normals[ki] = json_norm[ki].asDouble();
    }

    Json::Value json_ind = root.get("indices", "ASCII");
    std::cout << "ind_size: " << json_ind.size() << std::endl;
    indices.resize(json_ind.size());
    for (int ki = 0; ki < indices.size(); ++ki)
    {
	indices[ki] = json_ind[ki].asDouble();
    }

#if 0
    for (size_t ki = 0; ki < root.size(); ki++)
    {
	double vert = root[ki].get("vertices", "ASCII").asDouble();
	std::cout << "vert: " << vert << std::endl;
	// etc.
    }
#endif
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
	// The output
	std::cout << "Usage: triangulation.json vert_and_edges.g2" << std::endl;
	return 1;
    }

    std::ifstream filein_json(argv[1]);
    std::ofstream fileout(argv[2]);

    // We expect the input file to consist of vertices, normals and indices. We do not care about the normals.
    vector<double> vertices;
    vector<double> normals;
    vector<int> indices;
    readData(filein_json,
    	     vertices, normals, indices);

    const int dim = 3;
    const int num_vert = vertices.size()/dim;
    Go::PointCloud3D pt_cl(vertices.begin(), num_vert);
    pt_cl.writeStandardHeader(fileout);
    pt_cl.write(fileout);

    std::cout << "indices.size(): " << indices.size() << std::endl;
    int num_triang = indices.size()/3;
    int num_lines = num_triang*3;
    vector<double> triang_vert;
    for (size_t ki = 0; ki < indices.size(); ki += 3)
    {
	vector<double> node1(vertices.begin() + indices[ki]*dim,
			     vertices.begin() + (indices[ki]+1)*dim);
	vector<double> node2(vertices.begin() + indices[ki+1]*dim,
			     vertices.begin() + (indices[ki+1]+1)*dim);
	vector<double> node3(vertices.begin() + indices[ki+2]*dim,
			     vertices.begin() + (indices[ki+2]+1)*dim);
	triang_vert.insert(triang_vert.end(), node1.begin(), node1.end());
	triang_vert.insert(triang_vert.end(), node2.begin(), node2.end());
	triang_vert.insert(triang_vert.end(), node2.begin(), node2.end());
	triang_vert.insert(triang_vert.end(), node3.begin(), node3.end());
	triang_vert.insert(triang_vert.end(), node3.begin(), node3.end());
	triang_vert.insert(triang_vert.end(), node1.begin(), node1.end());
    }
    Go::LineCloud line_cl(triang_vert.begin(), num_lines);
    line_cl.writeStandardHeader(fileout);
    line_cl.write(fileout);

}
