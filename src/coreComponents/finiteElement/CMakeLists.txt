#
# Specify all headers
#
set( finiteElement_headers
     FiniteElementDiscretization.hpp
     FiniteElementDiscretizationManager.hpp
     FiniteElementDispatch.hpp
     Kinematics.h
     kernelInterface/KernelBase.hpp
     kernelInterface/ImplicitKernelBase.hpp
     FiniteElementDiscretizationManager.hpp
     FiniteElementDispatch.hpp
     elementFormulations/FiniteElementBase.hpp
     elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp
     elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp
     elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp
     elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp
     elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp
     elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp
     elementFormulations/LagrangeBasis1.hpp
     elementFormulations/LagrangeBasis2.hpp
   )
#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     FiniteElementDiscretizationManager.cpp
   )

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

convert_filenames_to_full_paths(finiteElement_sources)
convert_filenames_to_full_paths(finiteElement_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${finiteElement_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${finiteElement_headers} PARENT_SCOPE)
set(CORE_DEPENDENCIES ${CORE_DEPENDENCIES} ${dependencyList} PARENT_SCOPE)

geosx_add_code_checks( PREFIX finiteElement )

add_subdirectory( unitTests )

