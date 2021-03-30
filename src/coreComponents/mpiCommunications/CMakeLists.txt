
#
# Specify all headers
#
set(mpiCommunications_headers
    CommunicationTools.hpp
    MpiWrapper.hpp
    NeighborCommunicator.hpp
    PartitionBase.hpp
    SpatialPartition.hpp
    NeighborData.hpp
   )


#
# Specify all sources
#
set(mpiCommunications_sources
    CommunicationTools.cpp
    MpiWrapper.cpp
    NeighborCommunicator.cpp
    PartitionBase.cpp
    SpatialPartition.cpp
   )

set( dependencyList trilinos )

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if( ENABLE_MPI )
    set( dependencyList ${dependencyList} mpi )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()
               
convert_filenames_to_full_paths(mpiCommunications_sources)
convert_filenames_to_full_paths(mpiCommunications_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${mpiCommunications_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${mpiCommunications_headers} PARENT_SCOPE)
set(CORE_DEPENDENCIES ${CORE_DEPENDENCIES} ${dependencyList} PARENT_SCOPE)

geosx_add_code_checks(PREFIX mpiCommunications )
