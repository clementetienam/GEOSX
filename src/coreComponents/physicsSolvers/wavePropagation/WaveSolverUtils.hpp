/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file WaveSolverUtils.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_

#include "mesh/utilities/ComputationalGeometry.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

struct WaveSolverUtils
{
  static constexpr real64 epsilonLoc = 1e-8;
  static constexpr real64 eps64 = std::numeric_limits< real64 >::epsilon();
  static constexpr real32 eps32 = std::numeric_limits< real32 >::epsilon();

  using EXEC_POLICY = parallelDevicePolicy< >;
  using wsCoordType = real32;

  GEOS_HOST_DEVICE
  static real32 evaluateRicker( real64 const time_n, real32 const f0, real32 const t0, localIndex const order )
  {
    real32 const delay = t0 > 0 ? t0 : 1 / f0;
    real32 pulse = 0.0;
    if( time_n <= -0.9 * delay || time_n >= 2.9 * delay )
    {
      return pulse;
    }

    real32 const alpha = -pow( f0 * M_PI, 2 );
    real32 const time_d = time_n - delay;
    real32 const gaussian = exp( alpha * pow( time_d, 2 ));
    localIndex const sgn = pow( -1, order + 1 );

    switch( order )
    {
      case 0:
        pulse = sgn * gaussian;
        break;
      case 1:
        pulse = sgn * (2 * alpha * time_d) * gaussian;
        break;
      case 2:
        pulse = sgn * (2 * alpha + 4 * pow( alpha, 2 ) * pow( time_d, 2 )) * gaussian;
        break;
      case 3:
        pulse = sgn * (12 * pow( alpha, 2 ) * time_d + 8 * pow( alpha, 3 ) * pow( time_d, 3 )) * gaussian;
        break;
      case 4:
        pulse = sgn * (12 * pow( alpha, 2 ) + 48 * pow( alpha, 3 ) * pow( time_d, 2 ) + 16 * pow( alpha, 4 ) * pow( time_d, 4 )) * gaussian;
        break;
      default:
        GEOS_ERROR( "This option is not supported yet, rickerOrder must be in range {0:4}" );
    }

    return pulse;
  }

  /**
   * @brief Initialize (clear) the trace file.
   */
  static void initTrace( char const * prefix,
                         string const & name,
                         bool const outputSeismoTrace,
                         localIndex const nReceivers,
                         arrayView1d< localIndex const > const receiverIsLocal )
  {
    if( !outputSeismoTrace ) return;

    string const outputDir = OutputBase::getOutputDirectory();
    RAJA::ReduceSum< ReducePolicy< serialPolicy >, localIndex > count( 0 );

    forAll< serialPolicy >( nReceivers, [=] ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        count += 1;
        string const fn = joinPath( outputDir, GEOS_FMT( "{}_{}_{:03}.txt", prefix, name, ircv ) );
        std::ofstream f( fn, std::ios::out | std::ios::trunc );
      }
    } );

    localIndex const total = MpiWrapper::sum( count.get() );
    GEOS_ERROR_IF( nReceivers != total, GEOS_FMT( ": Invalid distribution of receivers: nReceivers={} != MPI::sum={}.", nReceivers, total ) );
  }

  /**
   * @brief Convenient helper for 3D vectors calling 3 times the scalar version with only the sampled variable argument changed.
   */
  static void writeSeismoTraceVector( char const * prefix,
                                      string const & name,
                                      bool const outputSeismoTrace,
                                      localIndex const nReceivers,
                                      arrayView1d< localIndex const > const receiverIsLocal,
                                      localIndex const nsamplesSeismoTrace,
                                      arrayView2d< real32 const > const varAtReceiversx,
                                      arrayView2d< real32 const > const varAtReceiversy,
                                      arrayView2d< real32 const > const varAtReceiversz )
  {
    writeSeismoTrace( prefix, name, outputSeismoTrace, nReceivers, receiverIsLocal, nsamplesSeismoTrace, varAtReceiversx );
    writeSeismoTrace( prefix, name, outputSeismoTrace, nReceivers, receiverIsLocal, nsamplesSeismoTrace, varAtReceiversy );
    writeSeismoTrace( prefix, name, outputSeismoTrace, nReceivers, receiverIsLocal, nsamplesSeismoTrace, varAtReceiversz );
  }

  /**
   * @brief Write the seismo traces to a file.
   */
  static void writeSeismoTrace( char const * prefix,
                                string const & name,
                                bool const outputSeismoTrace,
                                localIndex const nReceivers,
                                arrayView1d< localIndex const > const receiverIsLocal,
                                localIndex const nsamplesSeismoTrace,
                                arrayView2d< real32 const > const varAtReceivers )
  {
    if( !outputSeismoTrace ) return;

    string const outputDir = OutputBase::getOutputDirectory();
    forAll< serialPolicy >( nReceivers, [=] ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        string const fn = joinPath( outputDir, GEOS_FMT( "{}_{}_{:03}.txt", prefix, name, ircv ) );
        std::ofstream f( fn, std::ios::app );
        if( f )
        {
          GEOS_LOG_RANK( GEOS_FMT( "Append to seismo trace file {}", fn ) );
          for( localIndex iSample = 0; iSample < nsamplesSeismoTrace; ++iSample )
          {
            // index - time - value
            f << iSample << " " << varAtReceivers[iSample][nReceivers] << " " << varAtReceivers[iSample][ircv] << std::endl;
          }
          f.close();
        }
        else
        {
          GEOS_WARNING( GEOS_FMT( "Failed to open output file {}", fn ) );
        }
      }
    } );
  }

  /**
   * @brief Compute the seismo traces.
   */
  static void computeSeismoTrace( real64 const time_n,
                                  real64 const dt,
                                  real64 const timeSeismo,
                                  localIndex const iSeismo,
                                  arrayView2d< localIndex const > const receiverNodeIds,
                                  arrayView2d< real64 const > const receiverConstants,
                                  arrayView1d< localIndex const > const receiverIsLocal,
                                  arrayView1d< real32 const > const var_np1,
                                  arrayView1d< real32 const > const var_n,
                                  arrayView2d< real32 > varAtReceivers )
  {
    real64 const time_np1 = time_n + dt;

    real32 const a1 = abs( dt ) < epsilonLoc ? 1.0 : (time_np1 - timeSeismo) / dt;
    real32 const a2 = 1.0 - a1;

    localIndex const nReceivers = receiverConstants.size( 0 );

    forAll< EXEC_POLICY >( nReceivers, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        real32 vtmp_np1 = 0.0, vtmp_n = 0.0;
        for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
        {
          vtmp_np1 += var_np1[receiverNodeIds( ircv, inode )] * receiverConstants( ircv, inode );
          vtmp_n += var_n[receiverNodeIds( ircv, inode )] * receiverConstants( ircv, inode );
        }
        // linear interpolation between the pressure value at time_n and time_{n+1}
        varAtReceivers( iSeismo, ircv ) = a1 * vtmp_n + a2 * vtmp_np1;
        // NOTE: varAtReceivers has size(1) = numReceiversGlobal + 1, this does not OOB
        // left in the forAll loop for sync issues since the following does not depend on `ircv`
        varAtReceivers( iSeismo, nReceivers ) = a1 * time_n + a2 * time_np1;
      }
    } );
  }

  static void compute2dVariableSeismoTrace( real64 const time_n,
                                            real64 const dt,
                                            localIndex const regionIndex,
                                            arrayView1d< localIndex const > const receiverRegion,
                                            real64 const timeSeismo,
                                            localIndex const iSeismo,
                                            arrayView1d< localIndex const > const rcvElem,
                                            arrayView2d< real64 const > const receiverConstants,
                                            arrayView1d< localIndex const > const receiverIsLocal,
                                            arrayView2d< real32 const > const var_np1,
                                            arrayView2d< real32 const > const var_n,
                                            arrayView2d< real32 > varAtReceivers )
  {
    real64 const time_np1 = time_n + dt;

    real32 const a1 = dt < epsilonLoc ? 1.0 : (time_np1 - timeSeismo) / dt;
    real32 const a2 = 1.0 - a1;

    localIndex const nReceivers = receiverConstants.size( 0 );

    forAll< EXEC_POLICY >( nReceivers, [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        if( receiverRegion[ircv] == regionIndex )
        {
          real32 vtmp_np1 = 0.0, vtmp_n = 0.0;
          for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
          {
            vtmp_np1 += var_np1( rcvElem[ircv], inode ) * receiverConstants( ircv, inode );
            vtmp_n += var_n( rcvElem[ircv], inode ) * receiverConstants( ircv, inode );
          }
          // linear interpolation between the pressure value at time_n and time_{n+1}
          varAtReceivers( iSeismo, ircv ) = a1 * vtmp_n + a2 * vtmp_np1;
          // NOTE: varAtReceivers has size(1) = numReceiversGlobal + 1, this does not OOB
          // left in the forAll loop for sync issues since the following does not depend on `ircv`
          varAtReceivers( iSeismo, nReceivers ) = a1 * time_n + a2 * time_np1;
        }
      }
    } );
  }

  /**
   * @brief Check if the source point is inside an element or not
   * @param numFacesPerElem number of face on an element
   * @param elemCenter array containing the center of the elements
   * @param faceNormal array containing the normal of all faces
   * @param faceCenter array containing the center of all faces
   * @param elemsToFaces map to get the global faces from element index and local face index
   * @param coords coordinate of the point
   * @return true if coords is inside the element
   */
  GEOS_HOST_DEVICE
  static bool
  locateSourceElement( real64 const numFacesPerElem,
                       real64 const (&elemCenter)[3],
                       arrayView2d< real64 const > const faceNormal,
                       arrayView2d< real64 const > const faceCenter,
                       arraySlice1d< localIndex const > const elemsToFaces,
                       real64 const (&coords)[3] )
  {
    //Loop over the element faces
    real64 tmpVector[3]{};
    for( localIndex kfe = 0; kfe < numFacesPerElem; ++kfe )
    {

      localIndex const iface = elemsToFaces[kfe];
      real64 faceCenterOnFace[3] = {faceCenter[iface][0],
                                    faceCenter[iface][1],
                                    faceCenter[iface][2]};
      real64 faceNormalOnFace[3] = {faceNormal[iface][0],
                                    faceNormal[iface][1],
                                    faceNormal[iface][2]};

      //Test to make sure if the normal is outwardly directed
      LvArray::tensorOps::copy< 3 >( tmpVector, faceCenterOnFace );
      LvArray::tensorOps::subtract< 3 >( tmpVector, elemCenter );
      if( LvArray::tensorOps::AiBi< 3 >( tmpVector, faceNormalOnFace ) < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormalOnFace, -1 );
      }

      // compute the vector face center to query point
      LvArray::tensorOps::subtract< 3 >( faceCenterOnFace, coords );
      localIndex const s = computationalGeometry::sign( LvArray::tensorOps::AiBi< 3 >( faceNormalOnFace, faceCenterOnFace ));

      // all dot products should be non-negative (we enforce outward normals)
      if( s < 0 ) return false;

    }
    return true;
  }

  /**
   * @brief Convert a mesh element point coordinate into a coordinate on the reference element
   * @tparam FE_TYPE finite element type
   * @param[in] coords coordinate of the point
   * @param[in] elemsToNodes map to obtaint global nodes from element index
   * @param[in] nodeCoords array of mesh nodes coordinates
   * @param[out] coordsOnRefElem to contain the coordinate computed in the reference element
   */
  template< typename FE_TYPE >
  GEOS_HOST_DEVICE
  static void
  computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                        arraySlice1d< localIndex const, cells::NODE_MAP_USD - 1 > const elemsToNodes,
                                        arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                                        real64 (& coordsOnRefElem)[3] )
  {
    real64 xLocal[FE_TYPE::numNodes][3]{};
    for( localIndex a = 0; a < FE_TYPE::numNodes; ++a )
    {
      LvArray::tensorOps::copy< 3 >( xLocal[a], nodeCoords[ elemsToNodes[a] ] );
    }
    // coordsOnRefElem = invJ*(coords-coordsNode_0)
    real64 invJ[3][3]{};
    FE_TYPE::invJacobianTransformation( 0, xLocal, invJ );
    for( localIndex i = 0; i < 3; ++i )
    {
      // init at (-1,-1,-1) as the origin of the referential elem
      coordsOnRefElem[i] = -1.0;
      for( localIndex j = 0; j < 3; ++j )
      {
        coordsOnRefElem[i] += invJ[i][j] * (coords[j] - xLocal[0][j]);
      }
    }
  }

};

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERUTILS_HPP_ */
