/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "constitutive/relativePermeability/layouts.hpp"
#include "constitutive/capillaryPressure/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFVMKernels
{

using namespace constitutive;

/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from density, viscosity and relperm
 */
struct PhaseMobilityKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice2d< real64 const, compflow::USD_COMP_PROP_DC-1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE_PROP-2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE_PROP-2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_PROP_DC-2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE_PROP-2 > const & phaseVisc,
           arraySlice1d< real64 const, multifluid::USD_PHASE_PROP-2 > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_PROP_DC-2 > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const, relperm::USD_RELPERM-2 > const & phaseRelPerm,
           arraySlice2d< real64 const, relperm::USD_RELPERM_DS-2 > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE_PROP-1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_PROP_DC-1 > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE_PROP-1 > const & phaseMob,
           arraySlice1d< real64, compflow::USD_PHASE_PROP-1 > const & dPhaseMob_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_PROP_DC-1 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView3d< real64 const, compflow::USD_COMP_PROP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_PROP_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & dPhaseVisc_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_PROP_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE_PROP > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_PROP_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE_PROP > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE_PROP > const & dPhaseMob_dPres,
          arrayView3d< real64, compflow::USD_PHASE_PROP_DC > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const, compflow::USD_COMP_PROP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_PROP_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE_PROP > const & dPhaseVisc_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_PROP_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE_PROP > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_PROP_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE_PROP > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE_PROP > const & dPhaseMob_dPres,
          arrayView3d< real64, compflow::USD_PHASE_PROP_DC > const & dPhaseMob_dComp );
};


/******************************** FluxKernel ********************************/

/**
 * @brief Functions to assemble flux term contributions to residual and Jacobian
 */
struct FluxKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const stencilSize,
           localIndex const numPhases,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           arraySlice1d< real64 const > const stencilWeights,
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE_PROP > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE_PROP > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_PROP_DC > > const & dPhaseMob_dComp,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE_PROP > > const & dPhaseVolFrac_dPres,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_PROP_DC > > const & dPhaseVolFrac_dComp,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_PROP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE_PROP > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE_PROP > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_PROP_DC > > const & dPhaseMassDens_dComp,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dComp,
           ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d< real64 > const localFlux,
           arraySlice2d< real64 > const localFluxJacobian );

  template< localIndex NC, typename STENCIL_TYPE >
  static void
  launch( localIndex const numPhases,
          STENCIL_TYPE const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE_PROP > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE_PROP > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_PROP_DC > > const & dPhaseMob_dComp,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE_PROP > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_PROP_DC > > const & dPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_PROP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE_PROP > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE_PROP > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_PROP_DC > > const & dPhaseMassDens_dComp,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const, cappres::USD_CAPPRES > > const & phaseCapPressure,
          ElementViewConst< arrayView4d< real64 const, cappres::USD_CAPPRES_DS > > const & dPhaseCapPressure_dPhaseVolFrac,
          integer const capPressureFlag,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

} // namespace CompositionalMultiphaseFVMKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFVMKERNELS_HPP
