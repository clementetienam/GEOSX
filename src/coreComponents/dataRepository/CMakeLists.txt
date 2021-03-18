#
# Specify all headers
#
set(dataRepository_headers
    BufferOps_inline.hpp
    BufferOps.hpp
    BufferOpsDevice.hpp
    ExecutableGroup.hpp
    Group.hpp
    InputFlags.hpp
    KeyIndexT.hpp
    MappedVector.hpp
    ObjectCatalog.hpp
    ReferenceWrapper.hpp
    RestartFlags.hpp
    ConduitRestart.hpp
    wrapperHelpers.hpp
    Wrapper.hpp
    WrapperBase.hpp
    xmlWrapper.hpp
    )

#
# Specify all sources
#
set(dataRepository_sources
    BufferOpsDevice.cpp
    ExecutableGroup.cpp
    Group.cpp
    ConduitRestart.cpp
    Wrapper.cpp
    WrapperBase.cpp
    xmlWrapper.cpp
    )

set( extraComponentsLinkList "")

if( ENABLE_MPI )
  set( extraComponentsLinkList ${extraComponentsLinkList} mpi )
endif()

set (dependencyList pugixml )

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

set( dependencyList ${dependencyList} conduit )

if( ENABLE_PYGEOSX )
  set( dependencyList ${dependencyList} pylvarray )
endif()

if( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

convert_filenames_to_full_paths(dataRepository_sources)
convert_filenames_to_full_paths(dataRepository_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${dataRepository_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${dataRepository_headers} PARENT_SCOPE)
set(CORE_DEPENDENCIES ${CORE_DEPENDENCIES} ${dependencyList} PARENT_SCOPE)

geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
