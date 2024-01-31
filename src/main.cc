// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief Test for the richards facet coupling model.
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include "properties.hh"

// main program
int main(int argc, char** argv)
{
    using namespace Dumux;

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////
    // try to create the grids (from the given grid file)
    //////////////////////////////////////////////////////
    using MatrixProblemTypeTag = Properties::TTag::MATRIXTYPETAG;
    using LowDimProblemTypeTag = Properties::TTag::LOWDIMTYPETAG;
    using MatrixGrid = GetPropType<MatrixProblemTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimProblemTypeTag, Properties::Grid>;

    using GridManager = FacetCouplingGridManager<MatrixGrid, LowDimGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary, non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid views
    const auto& matrixGridView = gridManager.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridManager.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    using MatrixFVGridGeometry = GetPropType<MatrixProblemTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimProblemTypeTag, Properties::GridGeometry>;
    auto matrixFvGridGeometry = std::make_shared<MatrixFVGridGeometry>(matrixGridView);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);

    // the coupling mapper
    using TestTraits = Properties::TestTraits<MatrixProblemTypeTag, LowDimProblemTypeTag>;
    auto couplingMapper = std::make_shared<typename TestTraits::CouplingMapper>();
    couplingMapper->update(*matrixFvGridGeometry, *lowDimFvGridGeometry, gridManager.getEmbeddings());

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MatrixProblem = GetPropType<MatrixProblemTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimProblemTypeTag, Properties::Problem>;
    auto matrixSpatialParams = std::make_shared<typename MatrixProblem::SpatialParams>(matrixFvGridGeometry, "Matrix");
    auto matrixProblem = std::make_shared<MatrixProblem>(matrixFvGridGeometry, matrixSpatialParams, couplingManager, "Matrix");
    auto lowDimSpatialParams = std::make_shared<typename LowDimProblem::SpatialParams>(lowDimFvGridGeometry, "LowDim");
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, "LowDim");

    // the solution vector
    using MDTraits = typename TestTraits::MDTraits;
    using SolutionVector = typename MDTraits::SolutionVector;
    SolutionVector x;

    static const auto matrixId = typename MDTraits::template SubDomain<0>::Index();
    static const auto lowDimId = typename MDTraits::template SubDomain<1>::Index();
    x[matrixId].resize(matrixFvGridGeometry->numDofs());
    x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
    matrixProblem->applyInitialSolution(x[matrixId]);
    lowDimProblem->applyInitialSolution(x[lowDimId]);

    // initialize coupling manager
    couplingManager->init(matrixProblem, lowDimProblem, couplingMapper, x);

    // the grid variables
    using MatrixGridVariables = GetPropType<MatrixProblemTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimProblemTypeTag, Properties::GridVariables>;
    auto matrixGridVariables = std::make_shared<MatrixGridVariables>(matrixProblem, matrixFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    matrixGridVariables->init(x[matrixId]);
    lowDimGridVariables->init(x[lowDimId]);

    // initialize the vtk output module
    using MatrixSolutionVector = std::decay_t<decltype(x[matrixId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<MatrixGridVariables, MatrixSolutionVector> matrixVtkWriter(*matrixGridVariables, x[matrixId], matrixProblem->name(), "Matrix");
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim");

    // Add model specific output fields
    using MatrixIOFields = GetPropType<MatrixProblemTypeTag, Properties::IOFields>;
    using LowIOFields = GetPropType<LowDimProblemTypeTag, Properties::IOFields>;
    MatrixIOFields::initOutputModule(matrixVtkWriter);
    LowIOFields::initOutputModule(lowDimVtkWriter);

    // write initial solution
    matrixVtkWriter.write(0.0);
    lowDimVtkWriter.write(0.0);

    // the assembler
    using Assembler = MultiDomainFVAssembler<MDTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(matrixProblem, lowDimProblem),
                                                  std::make_tuple(matrixFvGridGeometry, lowDimFvGridGeometry),
                                                  std::make_tuple(matrixGridVariables, lowDimGridVariables),
                                                  couplingManager);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // linearize & solve
    newtonSolver->solve(x);

    // update grid variables for output
    matrixGridVariables->update(x[matrixId]);
    lowDimGridVariables->update(x[lowDimId]);

    // write vtk output
    matrixVtkWriter.write(1.0);
    lowDimVtkWriter.write(1.0);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    // output some
    newtonSolver->report();

    // output used run time parameters
    Parameters::print();

    return 0;
}
