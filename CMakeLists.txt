cmake_minimum_required(VERSION 3.10)

#project name
project(make_grid)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)

#specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

find_package(CGAL)
find_package(GUROBI REQUIRED)

message("BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
message("GUROBI_HOME: " $ENV{GUROBI_HOME})
message("GUROBI_INCLUDE: "${GUROBI_INCLUDE_DIRS} )
message("GUROBI_LIBS: "${GUROBI_LIBRARY} " " ${GUROBI_CXX_LIBRARY})

add_executable(${PROJECT_NAME} main.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC
	./
	code/
	external/Cinolib/external/eigen/
	external/Cinolib/include/
	${GUROBI_INCLUDE_DIRS}
)

if(UNIX AND NOT APPLE)
	target_link_libraries(${CMAKE_PROJECT_NAME} 
		$ENV{GUROBI_HOME}/lib/libgurobi_g++5.2.a
	)
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} 
	${GUROBI_CXX_LIBRARY}
	${GUROBI_LIBRARY}
	CGAL::CGAL
)


##################################################################################################