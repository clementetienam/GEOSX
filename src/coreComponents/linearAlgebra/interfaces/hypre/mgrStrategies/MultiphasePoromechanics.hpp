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
 * @file MultiphasePoromechanics.hpp
 */

#ifndef GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICS_HPP_
#define GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICS_HPP_

#include "linearAlgebra/interfaces/hypre/HypreMGR.hpp"

namespace geos
{

namespace hypre
{

namespace mgr
{

/**
 * @brief MultiphasePoromechanics strategy.
 *
 * Labels description stored in point_marker_array
 *   - dofLabel: 0             = nodal displacement, x-component
 *   - dofLabel: 1             = nodal displacement, y-component
 *   - dofLabel: 2             = nodal displacement, z-component
 *   - dofLabel: 3             = pressure
 *   - dofLabel: 4             = density
 *             ...             = densities
 *   - dofLabel: numLabels - 1 = density
 *
 * 3-level MGR reduction strategy based on CompositionalMultiphaseFVM
 *   - 1st level: eliminate displacements (0,1,2)
 *   - 2nd level: eliminate the reservoir density associated with the volume constraint (numLabels - 1)
 *   - 3nd level: eliminate the other reservoir densities
 *   - The coarse grid is solved with BoomerAMG.
 *
 */
class MultiphasePoromechanics : public MGRStrategyBase< 3 >
{
public:
  /**
   * @brief Constructor.
   * @param numComponentsPerField array with number of components for each field
   */
  explicit MultiphasePoromechanics( arrayView1d< int const > const & numComponentsPerField )
    : MGRStrategyBase( LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] + numComponentsPerField[1] ) )
  {
    // Level 0: eliminate displacement degrees of freedom
    m_labels[0].resize( m_numBlocks - 3 );
    std::iota( m_labels[0].begin(), m_labels[0].end(), 3 );
    // Level 1: eliminate last density which corresponds to the volume constraint equation
    m_labels[1].resize( m_numBlocks - 4 );
    std::iota( m_labels[1].begin(), m_labels[1].end(), 3 );
    // Level 2: eliminate the remaining reservoir densities
    m_labels[2].push_back( 3 ); // pressure

    setupLabels();

    // Level 0
    m_levelFRelaxType[0]         = MGRFRelaxationType::amgVCycle;
    m_levelFRelaxIters[0]        = 1;
    m_levelInterpType[0]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[0]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[0]   = MGRCoarseGridMethod::nonGalerkin;
    m_levelGlobalSmootherType[0] = MGRGlobalSmootherType::none;

    // Level 1
    m_levelFRelaxType[1]         = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[1]        = 1;
    m_levelInterpType[1]         = MGRInterpolationType::jacobi;
    m_levelRestrictType[1]       = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[1]   = MGRCoarseGridMethod::galerkin;
    m_levelGlobalSmootherType[1] = MGRGlobalSmootherType::none;

    // Level 2
    m_levelFRelaxType[2]          = MGRFRelaxationType::jacobi;
    m_levelFRelaxIters[2]         = 1;
    m_levelInterpType[2]          = MGRInterpolationType::injection;
    m_levelRestrictType[2]        = MGRRestrictionType::injection;
    m_levelCoarseGridMethod[2]    = MGRCoarseGridMethod::cprLikeBlockDiag;
    m_levelGlobalSmootherType[2]  = MGRGlobalSmootherType::ilu0;
    m_levelGlobalSmootherIters[2] = 1;
  }

  /**
   * @brief Setup the MGR strategy.
   * @param precond preconditioner wrapper
   * @param mgrData auxiliary MGR data
   */
  void setup( LinearSolverParameters::MGR const &,
              HyprePrecWrapper & precond,
              HypreMGRData & mgrData )
  {
    setReduction( precond, mgrData );

    GEOS_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( precond.ptr, 0 ));

    // Configure the BoomerAMG solver used as F-relaxation for the first level
    setMechanicsFSolver( precond, mgrData );

    // Configure the BoomerAMG solver used as mgr coarse solver for the pressure reduced system
    setPressureAMG( mgrData.coarseSolver );
  }
};

} // namespace mgr

} // namespace hypre

} // namespace geos

#endif /*GEOS_LINEARALGEBRA_INTERFACES_HYPREMGRMULTIPHASEPOROMECHANICS_HPP_*/
