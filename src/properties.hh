// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The properties for the bulk domain in the rinchards facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_RICHARDS_MATRIX_PROPERTIES_HH
#define DUMUX_TEST_FACETCOUPLING_RICHARDS_MATRIX_PROPERTIES_HH

// #include <dumux/material/components/constant.hh>
// #include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

// Matrix sub-problem
// the spatial parameters (permeabilities, material parameters etc.)
#include "matrixspatialparams.hh"
// we use alu grid for the discretization of the matrix domain
#include <dune/alugrid/grid.hh>

// LowDim sub-problem
// the spatial parameters (permeabilities, material parameters etc.)
#include "lowdimspatialparams.hh"
// we use foam grid for the discretization of the fracture domain
// as this grid manager is able to represent network/surface grids
#include <dune/foamgrid/foamgrid.hh>


#include "problem_matrix.hh"
#include "problem_lowdim.hh"

// default for the bulk grid type
#ifndef MATRIXGRIDTYPE
#define MATRIXGRIDTYPE Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>
#endif

#ifndef LOWDIMGRIDTYPE
#define LOWDIMGRIDTYPE Dune::FoamGrid<1, 2>
#endif

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct RichardsMatrix { using InheritsFrom = std::tuple<Richards>; };
struct RichardsMatrixTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, RichardsMatrix>; };
struct RichardsMatrixMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, RichardsMatrix>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsMatrix> { using type = MATRIXGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsMatrix> { using type = RichardsMatrixProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsMatrix>
{
    using type = MatrixSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// create the type tag nodes
namespace TTag {
struct RichardsLowDim { using InheritsFrom = std::tuple<Richards>; };
struct RichardsLowDimTpfa { using InheritsFrom = std::tuple<RichardsLowDim, CCTpfaModel>; };

// we need an additional type tag for the test using mpfa in the bulk domain
struct RichardsLowDimMpfa { using InheritsFrom = std::tuple<RichardsLowDim, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::RichardsLowDim> { using type = LOWDIMGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::RichardsLowDim> { using type = RichardsLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RichardsLowDim>
{
    using type = LowDimSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// obtain/define some types to be used below in the property definitions and in main
template< class MatrixTypeTag, class LowDimTypeTag >
class TestTraits
{
    using MatrixFVGridGeometry = GetPropType<MatrixTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<MatrixTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<MatrixFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set cm property in the sub-problems
using TpfaTraits = TestTraits<TTag::RichardsMatrixTpfa, TTag::RichardsLowDimTpfa>;
using MpfaTraits = TestTraits<TTag::RichardsMatrixMpfa, TTag::RichardsLowDimMpfa>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::RichardsMatrixTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::RichardsLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::RichardsMatrixMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::RichardsLowDimMpfa> { using type = typename MpfaTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
