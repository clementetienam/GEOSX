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
 * @file layouts.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_LAYOUTS_HPP
#define GEOSX_CONSTITUTIVE_FLUID_LAYOUTS_HPP

#include "common/GeosxConfig.hpp"

#include "LvArray/src/typeManipulation.hpp"
#include "RAJA/RAJA.hpp"

namespace geosx
{
namespace constitutive
{
namespace multifluid
{

#if defined( GEOSX_USE_CUDA )

/// Constitutive model phase property array layout
using LAYOUT_PHASE_PROP = RAJA::PERM_JKI;
/// Constitutive model phase property compositional derivative array layout
using LAYOUT_PHASE_PROP_DC = RAJA::PERM_JKLI;

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_JKLI;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_JKLMI;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID_PROP = RAJA::PERM_JI;
/// Constitutive model fluid property compositional derivative array layout
using LAYOUT_FLUID_PROP_DC = RAJA::PERM_JKI;

#else

/// Constitutive model phase property array layout
using LAYOUT_PHASE_PROP = RAJA::PERM_IJK;
/// Constitutive model phase property compositional derivative array layout
using LAYOUT_PHASE_PROP_DC = RAJA::PERM_IJKL;

/// Constitutive model phase composition array layout
using LAYOUT_PHASE_COMP = RAJA::PERM_IJKL;
/// Constitutive model phase composition compositional derivative array layout
using LAYOUT_PHASE_COMP_DC = RAJA::PERM_IJKLM;

/// Constitutive model fluid property array layout
using LAYOUT_FLUID_PROP = RAJA::PERM_IJ;
/// Constitutive model fluid property compositional derivative array layout
using LAYOUT_FLUID_PROP_DC = RAJA::PERM_IJK;

#endif

/// Constitutive model phase property unit stride dimension
static constexpr int USD_PHASE_PROP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_PROP{} );
/// Constitutive model phase property compositional derivative unit stride dimension
static constexpr int USD_PHASE_PROP_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_PROP_DC{} );

/// Constitutive model phase composition unit stride dimension
static constexpr int USD_PHASE_COMP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP{} );
/// Constitutive model phase composition compositional derivative unit stride dimension
static constexpr int USD_PHASE_COMP_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_PHASE_COMP_DC{} );

/// Constitutive model fluid property unit stride dimension
static constexpr int USD_FLUID_PROP = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_PROP{} );
/// Constitutive model fluid property compositional derivative unit stride dimension
static constexpr int USD_FLUID_PROP_DC = LvArray::typeManipulation::getStrideOneDimension( LAYOUT_FLUID_PROP_DC{} );

} // namespace multifluid
} // namespace constitutive
} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_LAYOUTS_HPP
