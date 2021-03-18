/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MultiscaleMeshUtils.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHUTILS_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHUTILS_HPP

#include "common/DataTypes.hpp"
#include "mesh/ObjectManagerBase.hpp"

namespace geosx
{
namespace multiscale
{

class MultiscaleMeshObjectManager;

namespace meshUtils
{

template< typename T >
void filterArray( arrayView1d< T const > const & src,
                  arrayView1d< T const > const & map,
                  array1d< T > & dst )
{
  dst.reserve( src.size() );
  for( T const & val : src )
  {
    T const newVal = map[val];
    if( newVal >= 0 )
    {
      dst.emplace_back( newVal );
    }
  }
}

template< typename T >
void filterArrayUnique( arrayView1d< T const > const & src,
                        arrayView1d< T const > const & map,
                        array1d< T > & dst )
{
  SortedArray< T > values;
  values.reserve( src.size() );
  for( T const & val : src )
  {
    T const newVal = map[val];
    if( newVal >= 0 )
    {
      values.insert( newVal );
    }
  }
  for( T const & val : values )
  {
    dst.emplace_back( val );
  }
}

template< typename T >
void filterSet( SortedArrayView< T const > const & src,
                arrayView1d< T const > const & map,
                SortedArray< T > & dst )
{
  dst.reserve( src.size() );
  for( T const & val : src )
  {
    T const newVal = map[val];
    if( newVal >= 0 )
    {
      dst.insert( newVal );
    }
  }
}

template< integer MAX_NEIGHBORS, typename L2C_MAP_TYPE, typename C2L_MAP_TYPE, typename LAMBDA >
void forUniqueNeighbors( localIndex const locIdx,
                         L2C_MAP_TYPE const & locToConn,
                         C2L_MAP_TYPE const & connToLoc,
                         integer const minCommonConn,
                         LAMBDA lambda )
{
  localIndex neighbors[MAX_NEIGHBORS];
  integer numNeighbors = 0;
  for( localIndex connIdx : locToConn[locIdx] )
  {
    for( localIndex nbrIdx : connToLoc[connIdx] )
    {
      GEOSX_ERROR_IF_GE_MSG( numNeighbors, MAX_NEIGHBORS, "Too many neighbors, need to increase stack limit" );
      neighbors[numNeighbors++] = nbrIdx;
    }
  }

  if( numNeighbors == 0 )
  {
    return;
  }
  LvArray::sortedArrayManipulation::makeSorted( neighbors, neighbors + numNeighbors );

  integer numConn = 0;
  for( integer i = 0; i < numNeighbors - 1; ++i )
  {
    ++numConn;
    if( neighbors[i + 1] != neighbors[i] )
    {
      if( numConn >= minCommonConn )
      {
        lambda( neighbors[i] );
      }
      numConn = 0;
    }
  }
  if( ++numConn >= minCommonConn )
  {
    lambda( neighbors[numNeighbors - 1] );
  }
}

template< typename FUNC >
void copyNeighborData( ObjectManagerBase const & srcManager,
                       string const & mapKey,
                       std::vector< integer > const & ranks,
                       ObjectManagerBase & dstManager,
                       FUNC copyFunc )
{
  arrayView1d< localIndex const > const map = srcManager.getReference< array1d< localIndex > >( mapKey );
  for( integer const rank : ranks )
  {
    NeighborData const & srcData = srcManager.getNeighborData( rank );
    NeighborData & dstData = dstManager.getNeighborData( rank );
    copyFunc( srcData.ghostsToSend(), map, dstData.ghostsToSend() );
    copyFunc( srcData.ghostsToReceive(), map, dstData.ghostsToReceive() );
    copyFunc( srcData.adjacencyList(), map, dstData.adjacencyList() );
    copyFunc( srcData.matchedPartitionBoundary(), map, dstData.matchedPartitionBoundary() );
  }
}

void copySets( ObjectManagerBase const & srcManager,
               string const & mapKey,
               ObjectManagerBase & dstManager );

} // namespace meshUtils
} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHUTILS_HPP
