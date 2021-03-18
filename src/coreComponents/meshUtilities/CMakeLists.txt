
#
# Specify all headers
#
set(meshUtilities_headers
    ComputationalGeometry.hpp
    MeshManager.hpp
    MeshGeneratorBase.hpp
    InternalMeshGenerator.hpp
    InternalWellGenerator.hpp
    PerforationData.hpp
    Perforation.hpp
    MeshUtilities.hpp
    SimpleGeometricObjects/GeometricObjectManager.hpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.hpp
    SimpleGeometricObjects/Box.hpp
    SimpleGeometricObjects/ThickPlane.hpp
    SimpleGeometricObjects/BoundedPlane.hpp
    StructuredGridUtilities.hpp
   )

#
# Specify all sources
#
set(meshUtilities_sources
    ComputationalGeometry.cpp
    MeshManager.cpp
    MeshGeneratorBase.cpp
    InternalMeshGenerator.cpp
    InternalWellGenerator.cpp
    PerforationData.cpp
    Perforation.cpp
    MeshUtilities.cpp
    SimpleGeometricObjects/GeometricObjectManager.cpp
    SimpleGeometricObjects/SimpleGeometricObjectBase.cpp
    SimpleGeometricObjects/Box.cpp
    SimpleGeometricObjects/Cylinder.cpp
    SimpleGeometricObjects/ThickPlane.cpp
    SimpleGeometricObjects/BoundedPlane.cpp
   )

if( ENABLE_PAMELA )
    message(STATUS "Adding PAMELAMeshGenerator sources and headers")
    set( meshUtilities_headers ${meshUtilities_headers} PAMELAMeshGenerator.hpp )
    set( meshUtilities_sources ${meshUtilities_sources} PAMELAMeshGenerator.cpp)
    set( dependencyList ${dependencyList} PAMELA )
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

convert_filenames_to_full_paths(meshUtilities_sources)
convert_filenames_to_full_paths(meshUtilities_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${meshUtilities_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${meshUtilities_headers} PARENT_SCOPE)
set(CORE_DEPENDENCIES ${CORE_DEPENDENCIES} ${dependencyList} PARENT_SCOPE)

add_subdirectory( unitTests )

geosx_add_code_checks(PREFIX meshUtilities )
