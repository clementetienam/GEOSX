/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// using some utility classes from the following unit test
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "common/DataTypes.hpp"
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/wavePropagation/WaveSolverBase.hpp"
#include "physicsSolvers/wavePropagation/WaveSolverBaseFields.hpp"
#include "physicsSolvers/wavePropagation/AcousticFirstOrderWaveEquationSEM.hpp"

#include <gtest/gtest.h>

using namespace geos;
using namespace geos::dataRepository;
using namespace geos::testing;

CommandLineOptions g_commandLineOptions;

// This unit test checks the interpolation done to extract seismic traces from a wavefield.
// It computes a seismogram at a receiver co-located with the source and compares it to the surrounding receivers.
char const * xmlInput =
  R"xml(
  <?xml version="1.0" ?>
  <Problem>
    <Solvers>
      <AcousticFirstOrderSEM
        name="acousticFirstOrderSolver"
        cflFactor="0.25"
        discretization="FE1"
        targetRegions="{ Region }"
        sourceCoordinates="{ { 30, 30, 30 } }"
        timeSourceFrequency="2"
        receiverCoordinates="{ { 0.1, 0.1, 0.1 }, { 0.1, 0.1, 99.9 }, { 0.1, 99.9, 0.1 }, { 0.1, 99.9, 99.9 },
                                { 99.9, 0.1, 0.1 }, { 99.9, 0.1, 99.9 }, { 99.9, 99.9, 0.1 }, { 99.9, 99.9, 99.9 },
                                { 50, 50, 50 } }"
        outputSeismoTrace="0"
        dtSeismoTrace="0.05"
        rickerOrder="1"/>
    </Solvers>
    <Mesh>
      <InternalMesh
        name="mesh"
        elementTypes="{ C3D8 }"
        xCoords="{ 0, 100 }"
        yCoords="{ 0, 100 }"
        zCoords="{ 0, 100 }"
        nx="{ 1 }"
        ny="{ 1 }"
        nz="{ 1 }"
        cellBlockNames="{ cb }"/>
    </Mesh>
    <Events
      maxTime="1">
      <PeriodicEvent
        name="solverApplications"
        forceDt="0.05"
        targetExactStartStop="0"
        targetExactTimestep="0"
        target="/Solvers/acousticFirstOrderSolver"/>
      <PeriodicEvent
        name="waveFieldUxCollection"
        timeFrequency="0.05"
        targetExactTimestep="0"
        target="/Tasks/waveFieldUxCollection" />
      <PeriodicEvent
        name="waveFieldUyCollection"
        timeFrequency="0.05"
        targetExactTimestep="0"
        target="/Tasks/waveFieldUyCollection" />
      <PeriodicEvent
        name="waveFieldUzCollection"
        timeFrequency="0.05"
        targetExactTimestep="0"
        target="/Tasks/waveFieldUzCollection" />
      <PeriodicEvent
        name="waveFieldPressureCollection"
        timeFrequency="0.05"
        targetExactTimestep="0"
        target="/Tasks/waveFieldPressureCollection" />
    </Events>
    <NumericalMethods>
      <FiniteElements>
        <FiniteElementSpace
          name="FE1"
          order="1"
          formulation="SEM" />
      </FiniteElements>
    </NumericalMethods>
    <ElementRegions>
      <CellElementRegion
        name="Region"
        cellBlocks="{ cb }"
        materialList="{ nullModel }"/>
    </ElementRegions>
    <Constitutive>
      <NullModel
        name="nullModel"/>
    </Constitutive>
    <FieldSpecifications>
      <FieldSpecification
        name="cellVelocity"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="mediumVelocity"
        scale="1500"
        setNames="{ all }"/>
      <FieldSpecification
        name="cellDensity"
        initialCondition="1"
        objectPath="ElementRegions/Region/cb"
        fieldName="mediumDensity"
        scale="1"
        setNames="{ all }"/>
    </FieldSpecifications>
    <Tasks>
      <PackCollection
        name="waveFieldPressureCollection"
        objectPath="nodeManager"
        fieldName="pressure_np1"/>
      <PackCollection
        name="waveFieldUxCollection"
        objectPath="ElementRegions/Region/cb"
        fieldName="velocity_x"/>
      <PackCollection
        name="waveFieldUyCollection"
        objectPath="ElementRegions/Region/cb"
        fieldName="velocity_y"/>
      <PackCollection
        name="waveFieldUzCollection"
        objectPath="ElementRegions/Region/cb"
        fieldName="velocity_z"/>
    </Tasks>
  </Problem>
  )xml";

class AcousticFirstOrderWaveEquationSEMTest : public ::testing::Test
{
public:

  AcousticFirstOrderWaveEquationSEMTest():
    state( std::make_unique< CommandLineOptions >( g_commandLineOptions ) )
  {}

protected:

  void SetUp() override
  {
    setupProblemFromXML( state.getProblemManager(), xmlInput );
  }

  static real64 constexpr time = 0.0;
  static real64 constexpr dt = 5e-2;
  static real64 constexpr eps = std::numeric_limits< real64 >::epsilon();

  GeosxState state;
  AcousticFirstOrderWaveEquationSEM * propagator;
};

real64 constexpr AcousticFirstOrderWaveEquationSEMTest::time;
real64 constexpr AcousticFirstOrderWaveEquationSEMTest::dt;
real64 constexpr AcousticFirstOrderWaveEquationSEMTest::eps;

TEST_F( AcousticFirstOrderWaveEquationSEMTest, SeismoTrace )
{

  DomainPartition & domain = state.getProblemManager().getDomainPartition();
  propagator = &state.getProblemManager().getPhysicsSolverManager().getGroup< AcousticFirstOrderWaveEquationSEM >( "acousticFirstOrderSolver" );
  real64 time_n = time;
  // run for 1s (20 steps)
  for( int i=0; i<20; i++ )
  {
    propagator->solverStep( time_n, dt, i, domain );
    time_n += dt;
  }
  // cleanup (triggers calculation of the remaining seismograms data points)
  propagator->cleanup( 1.0, 20, 0, 0, domain );

  // retrieve seismo
  arrayView2d< real32 > const pReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::pressureNp1AtReceiversString() ).toView();
  arrayView2d< real32 > const uxReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::uxNp1AtReceiversString() ).toView();
  arrayView2d< real32 > const uyReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::uyNp1AtReceiversString() ).toView();
  arrayView2d< real32 > const uzReceivers = propagator->getReference< array2d< real32 > >( AcousticFirstOrderWaveEquationSEM::viewKeyStruct::uzNp1AtReceiversString() ).toView();

  // move it to CPU, if needed
  pReceivers.move( LvArray::MemorySpace::host, false );
  uxReceivers.move( LvArray::MemorySpace::host, false );
  uyReceivers.move( LvArray::MemorySpace::host, false );
  uzReceivers.move( LvArray::MemorySpace::host, false );

  // check number of seismos and trace length
  ASSERT_EQ( pReceivers.size( 1 ), 10 );
  ASSERT_EQ( pReceivers.size( 0 ), 21 );
  ASSERT_EQ( uxReceivers.size( 1 ), 10 );
  ASSERT_EQ( uxReceivers.size( 0 ), 21 );
  ASSERT_EQ( uyReceivers.size( 1 ), 10 );
  ASSERT_EQ( uyReceivers.size( 0 ), 21 );
  ASSERT_EQ( uzReceivers.size( 1 ), 10 );
  ASSERT_EQ( uzReceivers.size( 0 ), 21 );

  // check seismo content. The pressure and velocity values cannot be directly checked as the problem is too small.
  // Since the basis is linear, check that the seismograms are nonzero (for t>0) and the seismogram at the center is equal
  // to the average of the others.
  for( int i = 0; i < 21; i++ )
  {
    if( i > 0 )
    {
      ASSERT_TRUE( std::abs( pReceivers[i][8] ) > 0 );
      ASSERT_TRUE( std::abs( uxReceivers[i][8] ) > 0 );
      ASSERT_TRUE( std::abs( uyReceivers[i][8] ) > 0 );
      ASSERT_TRUE( std::abs( uzReceivers[i][8] ) > 0 );
    }
    double avgP = 0;
    double avgUx = 0;
    double avgUy = 0;
    double avgUz = 0;
    for( int r=0; r<8; r++ )
    {
      avgP += pReceivers[i][r];
      avgUx += uxReceivers[i][r];
      avgUy += uyReceivers[i][r];
      avgUz += uzReceivers[i][r];
    }
    avgP /= 8.0;
    avgUx /= 8.0;
    avgUy /= 8.0;
    avgUz /= 8.0;
    ASSERT_TRUE( std::abs( pReceivers[i][8] - avgP ) < 0.00001 );
    ASSERT_TRUE( std::abs( uxReceivers[i][8] - avgUx ) < 0.00001 );
    ASSERT_TRUE( std::abs( uyReceivers[i][8] - avgUy ) < 0.00001 );
    ASSERT_TRUE( std::abs( uzReceivers[i][8] - avgUz ) < 0.00001 );
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geos::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::basicCleanup();
  return result;
}
