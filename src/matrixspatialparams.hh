// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The spatial parameters for the richards facet coupling test.
 */
#ifndef DUMUX_TEST_MATRIX_RICHARDS_SPATIALPARAMS_HH
#define DUMUX_TEST_MATRIX_RICHARDS_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the richards facet coupling test.
 */
template< class GridGeometry, class Scalar >
class MatrixSpatialParams
: public FVPorousMediumFlowSpatialParamsMP< GridGeometry, Scalar, MatrixSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = MatrixSpatialParams< GridGeometry, Scalar >;
    using ParentType = FVPorousMediumFlowSpatialParamsMP< GridGeometry, Scalar, ThisType >;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // use a van-genuchten fluid matrix interaction relationship
    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    //! Export the type used for permeabilities
    using PermeabilityType = Scalar;

    MatrixSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup)
    : ParentType(gridGeometry)
    , pcKrSwCurve_("Matrix.SpatialParams")
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

    //! Returns the porosity
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return porosity_; }

    /*!
     * \brief Returns the fluid-matrix interaction law at a given location
     * \param globalPos The global coordinates for the given location
     */
    auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
    {
        return makeFluidMatrixInteraction(pcKrSwCurve_);
    }

    /*!
     * \brief Returns the temperature at the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 283.15;
    }


private:
    PermeabilityType permeability_;
    Scalar porosity_;
    const PcKrSwCurve pcKrSwCurve_;
};

} // end namespace Dumux

#endif
