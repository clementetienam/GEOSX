#
# Specify all headers
#
set( virtualElement_headers
  ConformingVirtualElementOrder1.hpp
  VirtualElementBase.hpp
  )
#
# Specify all sources
#
set( virtualElement_sources
  ConformingVirtualElementOrder1.cpp
  )

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

convert_filenames_to_full_paths(virtualElement_sources)
convert_filenames_to_full_paths(virtualElement_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${virtualElement_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${virtualElement_headers} PARENT_SCOPE)
set(CORE_DEPENDENCIES ${CORE_DEPENDENCIES} ${dependencyList} PARENT_SCOPE)

geosx_add_code_checks( PREFIX virtualElement )
