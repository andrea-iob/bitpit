/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

//Standard Template Library
# include <ctime>
# include <chrono>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_CG.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"
# include "bitpit_levelset.hpp"

/*!
* Subtest 001
*
* Testing basic features of a 3D levelset.
*/
int subtest_001()
{
    int dimensions(3) ;


    double distance;
    std::array<double, 3> projection;

    std::array<double, 3> point = {{56.925, 50.715, -100.567}};
    std::array<std::array<double, 3>, 3> vertexCoordinates;

    std::array<double, 3> lambda;


    vertexCoordinates[0] = {{55.2, 50.6, -100.841}};
    vertexCoordinates[1] = {{55.2, 50.6, -101.2}};
    vertexCoordinates[2] = {{55.2, 51.2504, -100.93}};

    projection = bitpit::CGElem::projectPointTriangle(point, vertexCoordinates[0], vertexCoordinates[1], vertexCoordinates[2]);
    distance   = norm2(point - projection);

    std::cout << "    point > " << point << std::endl;
    std::cout << "    V0 > " << vertexCoordinates[0] << std::endl;
    std::cout << "    V1 > " << vertexCoordinates[1] << std::endl;
    std::cout << "    V2 > " << vertexCoordinates[2] << std::endl;
    std::cout << "    lambda > " << lambda << std::endl;
    std::cout << "    distance > " << distance << std::endl;
    std::cout << "    projection > " << projection << std::endl;



    vertexCoordinates[0] = {{55.2, 50.6, -100.841}};
    vertexCoordinates[1] = {{55.2, 50.6, -101.2}};
    vertexCoordinates[2] = {{55.2, 51.2504, -100.93}};

    distance = bitpit::CGElem::distancePointTriangle(point, vertexCoordinates[0], vertexCoordinates[1], vertexCoordinates[2], lambda);
    projection = bitpit::CGElem::reconstructPointFromBarycentricTriangle(vertexCoordinates[0], vertexCoordinates[1], vertexCoordinates[2], lambda);


    std::cout << "    point > " << point << std::endl;
    std::cout << "    V0 > " << vertexCoordinates[0] << std::endl;
    std::cout << "    V1 > " << vertexCoordinates[1] << std::endl;
    std::cout << "    V2 > " << vertexCoordinates[2] << std::endl;
    std::cout << "    lambda > " << lambda << std::endl;
    std::cout << "    distance > " << distance << std::endl;
    std::cout << "    projection > " << projection << std::endl;















    vertexCoordinates[0] = {{55.2, 50.6, -100.841}};
    vertexCoordinates[1] = {{55.2, 51.2504, -100.93}};
    vertexCoordinates[2] = {{54.3115, 51.4885, -101.2}};

    projection = bitpit::CGElem::projectPointTriangle(point, vertexCoordinates[0], vertexCoordinates[1], vertexCoordinates[2]);
    distance   = norm2(point - projection);

    std::cout << "    point > " << point << std::endl;
    std::cout << "    V0 > " << vertexCoordinates[0] << std::endl;
    std::cout << "    V1 > " << vertexCoordinates[1] << std::endl;
    std::cout << "    V2 > " << vertexCoordinates[2] << std::endl;
    std::cout << "    lambda > " << lambda << std::endl;
    std::cout << "    distance > " << distance << std::endl;
    std::cout << "    projection > " << projection << std::endl;



    vertexCoordinates[0] = {{55.2, 50.6, -100.841}};
    vertexCoordinates[1] = {{55.2, 51.2504, -100.93}};
    vertexCoordinates[2] = {{54.3115, 51.4885, -101.2}};


    distance = bitpit::CGElem::distancePointTriangle(point, vertexCoordinates[0], vertexCoordinates[1], vertexCoordinates[2], lambda);
    projection = bitpit::CGElem::reconstructPointFromBarycentricTriangle(vertexCoordinates[0], vertexCoordinates[1], vertexCoordinates[2], lambda);


    std::cout << "    point > " << point << std::endl;
    std::cout << "    V0 > " << vertexCoordinates[0] << std::endl;
    std::cout << "    V1 > " << vertexCoordinates[1] << std::endl;
    std::cout << "    V2 > " << vertexCoordinates[2] << std::endl;
    std::cout << "    lambda > " << lambda << std::endl;
    std::cout << "    distance > " << distance << std::endl;
    std::cout << "    projection > " << projection << std::endl;



//     exit(0);



















    // Input geometry
#if BITPIT_ENABLE_MPI
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(2, 3, MPI_COMM_NULL) );
#else
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(2, 3) );
#endif

    bitpit::log::cout()<< " - Loading stl geometry" << std::endl;

    STL->importSTL("./data/extracted.2.ascii.stl", bitpit::STLReader::FormatUnknown, true);

    STL->initializeAdjacencies() ;

    STL->getVTK().setName("geometry_002") ;
    STL->write() ;

    bitpit::log::cout()<< "n. vertex: " << STL->getVertexCount() << std::endl;
    bitpit::log::cout()<< "n. simplex: " << STL->getCellCount() << std::endl;


    // create cartesian mesh around geometry 
    bitpit::log::cout()<< " - Setting mesh" << std::endl;
    std::array<double,3> meshMin, meshMax, delta ;
    std::array<int,3> nc = {{64, 64, 64}} ;

    STL->getBoundingBox( meshMin, meshMax ) ;

    delta = meshMax -meshMin ;
    meshMin -=  0.1*delta ;
    meshMax +=  0.1*delta ;

    delta = meshMax -meshMin ;

    std::array<double, 3> origin = {{42., 12., 76.}};
    std::array<double, 3> length = {{18., 18., 18.}};
    std::array<int, 3> nnnn = {{18, 18, 18}};

    bitpit::VolCartesian mesh( 1, dimensions, origin, length, nnnn);
    mesh.update() ;
    mesh.initializeAdjacencies() ;

    // Compute level set  in narrow band
    bitpit::LevelSet levelset ;
    std::chrono::time_point<std::chrono::system_clock> start, end;

    levelset.setMesh(&mesh) ;
    int id0 = levelset.addObject( std::move(STL), BITPIT_PI/8. ) ;

    levelset.setPropagateSign(true) ;
    start = std::chrono::system_clock::now();
    levelset.compute( ) ;
    end = std::chrono::system_clock::now();

    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout()<< "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    bitpit::log::cout()<< " - Exporting data" << std::endl;
    mesh.getVTK().setName("levelset_002") ;
    bitpit::LevelSetObject &object = levelset.getObject(id0);

    object.enableVTKOutput( bitpit::LevelSetWriteField::ALL);

    mesh.write() ;

    return 0;
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
	MPI_Init(&argc,&argv);
#else
	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
#endif

	// Initialize the logger
	bitpit::log::manager().initialize(bitpit::log::COMBINED);

	// Run the subtests
	bitpit::log::cout() << "Testing basic levelset features" << std::endl;

	int status;
	try {
		status = subtest_001();
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		bitpit::log::cout() << exception.what();
		exit(1);
	}

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif
}
