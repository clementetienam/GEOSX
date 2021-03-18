# Specify all headers
set(mesh_headers
    AggregateElementSubRegion.hpp
    BufferOps.hpp
    CellBlockManager.hpp
    CellBlock.hpp
    CellElementRegion.hpp
    CellElementSubRegion.hpp
    EdgeManager.hpp
    ElementSubRegionBase.hpp
    ElementRegionBase.hpp
    ElementRegionManager.hpp
    SurfaceElementSubRegion.hpp
    FaceElementSubRegion.hpp
    SurfaceElementRegion.hpp
    EmbeddedSurfaceSubRegion.hpp
    WellElementSubRegion.hpp
    WellElementRegion.hpp
    ExtrinsicMeshData.hpp
    FaceManager.hpp
    InterObjectRelation.hpp
    MeshBody.hpp
    MeshLevel.hpp
    NodeManager.hpp
    ToElementRelation.hpp
  )

# Specify all sources
set(mesh_sources
    AggregateElementSubRegion.cpp
    BufferOps.cpp
    CellBlockManager.cpp
    CellBlock.cpp
    CellElementRegion.cpp
    CellElementSubRegion.cpp
    EdgeManager.cpp
    ElementSubRegionBase.cpp
    ElementRegionBase.cpp
    ElementRegionManager.cpp
    SurfaceElementRegion.cpp
    SurfaceElementSubRegion.cpp
    FaceElementSubRegion.cpp
    EmbeddedSurfaceSubRegion.cpp
    WellElementRegion.cpp
    WellElementSubRegion.cpp
    FaceManager.cpp
    MeshBody.cpp
    MeshLevel.cpp
    NodeManager.cpp
    ToElementRelation.cpp
   )

set(dependencyList
    metis
    )


if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

convert_filenames_to_full_paths(mesh_sources)
convert_filenames_to_full_paths(mesh_headers)

set(CORE_SOURCES ${CORE_SOURCES} ${mesh_sources} PARENT_SCOPE)
set(CORE_HEADERS ${CORE_HEADERS} ${mesh_headers} PARENT_SCOPE)
set(CORE_DEPENDENCIES ${CORE_DEPENDENCIES} ${dependencyList} PARENT_SCOPE)

geosx_add_code_checks(PREFIX mesh )
