#include "dune/biw4-67-phasefield/biw4-67-phasefield.hh"
#include "dune/istl/bvector.hh"
#include <config.h>

#include <array>
#include <dune/common/indices.hh>
#include <dune/functions/backends/istlvectorbackend.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/compositebasis.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/functions/functionspacebases/taylorhoodbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/gridviewfunction.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <vector>
using namespace Dune;

template <class LocalView, class LocalDispFunctionView,
          class LocalPhaseFunctionView, class Material>
void assembleElementStiffnessMatrix(
    const LocalView &localView, //
    const LocalDispFunctionView &localDispFunction, //
    const LocalPhaseFunctionView &localPhaseFunction, //
    Matrix<double> &elementMatrix, //
    const Material &material
    ) {

  using Element = typename LocalView::Element;
  const Element element = localView.element();
  constexpr int dim = Element::dimension;
  auto geometry = element.geometry();

  elementMatrix.setSize(localView.size(), localView.size());
  elementMatrix = 0;

  using namespace Indices;
  const auto &displacementLocalFiniteElement =
      localView.tree().child(_0, 0).finiteElement();
  const auto &phasefieldLocalFiniteElement =
      localView.tree().child(_1).finiteElement();
  auto num_nodes = displacementLocalFiniteElement.localBasis().size();
  int order =
      2 * (dim * displacementLocalFiniteElement.localBasis().order() - 1);
  const auto &quad = QuadratureRules<double, dim>::rule(element.type(), order);

  // Displacement Gradient
  auto dispDerivative = derivative(localDispFunction);
  // Loop over all quadrature points
  for (auto [position, weight] : quad) {
    const auto jacobianInverseTransposed =
        geometry.jacobianInverseTransposed(position);
    // The multiplicative factor in the integral transformation formula
    const auto integrationElement = geometry.integrationElement(position);

    // calculate relevant function values
    auto phasefieldFunctionValue = localPhaseFunction(position)[0];

    auto dispGradient = dispDerivative(position);
    Dune::FieldMatrix<double, dim, dim> strains(0);

    //std::cout << dispGradient << std::endl;

    //for (int i = 0; i < dim; i++) {
    //  for (int j = 0; j < dim; j++) {
    //    strains[i][j] = 0.5 * ( dispGradient[i][j] + dispGradient[j][i] );
    //  }
    //}

    //////////////////////////////////////////////////////////////////
    // The gradients of the shape functions on the reference element
    std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
    displacementLocalFiniteElement.localBasis().evaluateJacobian(
        position, referenceGradients);
    // Compute the shape function gradients on the grid element
    std::vector<FieldVector<double, dim>> gradients(referenceGradients.size());
    for (size_t i = 0; i < gradients.size(); i++)
      jacobianInverseTransposed.mv(referenceGradients[i][0], gradients[i]);

    std::vector<std::array<FieldMatrix<double, dim, dim>, dim>> deltaLinStrain(
        num_nodes);
    // Loop over the Dofs
    for (size_t i = 0; i < num_nodes; i++) {
      for (size_t k = 0; k < dim; k++) {
        //
        FieldMatrix<double, dim, dim> displacementGradiente(0);
        displacementGradiente[k] = gradients[i];
        // elementMatrix[row][col] += ( gradients[i] * gradients[j] )*
        // quadPoint.weight() * integrationElement;

        FieldMatrix<double, dim, dim> linearisedStrains(0);

        for (int k = 0; k < dim; k++) {
          for (int l = 0; l < dim; l++) {
            linearisedStrains[k][l] = 0.5 * (displacementGradiente[k][l] +
                                             displacementGradiente[l][k]);
          }
        }
        // std::cout << linearisedStrains << std::endl;

        deltaLinStrain[i][k] = linearisedStrains;
      }
    }

    /////////////////////////////////////////////////////////////////
    // The values of the Phasefield shape functions
    std::vector<FieldVector<double, 1>> phaseFieldValue;
    phasefieldLocalFiniteElement.localBasis().evaluateFunction(position,
                                                               phaseFieldValue);

    std::vector<FieldMatrix<double, 1, dim>> referencePhaseGradients;
    phasefieldLocalFiniteElement.localBasis().evaluateJacobian(
        position, referencePhaseGradients);

    std::vector<FieldVector<double, dim>> phaseFieldGradients(
        referencePhaseGradients.size());

    for (size_t i = 0; i < gradients.size(); i++)
      jacobianInverseTransposed.mv(referencePhaseGradients[i][0],
                                   phaseFieldGradients[i]);


    ///////////////////////////////////////////////////////////////////////
    // Displacements--Displacements coupling
    ///////////////////////////////////////////////////////////////////////
    // TODO degradation
    double degradation = material.degradationFunction(phasefieldFunctionValue);
    // Compute the actual matrix entries
    // \int g(\phi) * C(\Delta \eps) : \virt \eps \diff \Omega
    FieldMatrix<double, dim, dim> stressTensor(0);
    for (size_t i = 0; i < displacementLocalFiniteElement.size(); i++)
      for (size_t j = 0; j < displacementLocalFiniteElement.size(); j++)
        for (size_t k = 0; k < dim; k++) {
          stressTensor = 0.0;
          size_t row = localView.tree().child(_0, k).localIndex(i);
          size_t col = localView.tree().child(_0, k).localIndex(j);
          auto inkreStrain = deltaLinStrain[row][i];
          auto virtStrain = deltaLinStrain[col][j];
          material.stresses(inkreStrain, stressTensor);

          double value = 0.0;
          for (int l = 0; l < dim; l++)
            for (int o = 0; o < dim; o++) {
              value += virtStrain[l][o] * stressTensor[l][o];
            }

          elementMatrix[row][col] +=
              degradation * value * weight * integrationElement;
        }

    ///////////////////////////////////////////////////////////////////////
    // Displacements--Phasefield coupling
    ///////////////////////////////////////////////////////////////////////
    

    // Compute the actual matrix entries

    // TODO: Current Stresses, degradation function derivative
    // -\int 2(1 - \phi) \virt\phi \cdot \sigma_0(\eps) : \delt\eps \diff \Omega 
    for (size_t i = 0; i < displacementLocalFiniteElement.size(); i++)
      for (size_t k = 0; k < dim; k++)
        for (size_t j = 0; j < phasefieldLocalFiniteElement.size(); j++) {
          size_t vIndex = localView.tree().child(_0, k).localIndex(i);
          size_t pIndex = localView.tree().child(_1).localIndex(j);
          auto value = - 2.0* (1.0 - phasefieldFunctionValue) * phaseFieldValue[j] * weight *
                       integrationElement;
          elementMatrix[vIndex][pIndex] += value;
          elementMatrix[pIndex][vIndex] += value;
        }

    ///////////////////////////////////////////////////////////////////////
    // Phasefield--Phasefield coupling
    ///////////////////////////////////////////////////////////////////////

    // Compute the actual matrix entries
    // \int 2 \Delta \phi\virt \phi \cdot \psi_0 \diff \Omega 
    // + G_c \int \frac{1}{l} (\Delta \phi \vdiff \phi + l^2 \grad \Delta \phi \cdot \grad \vdiff \phi) \diff \Omega
    for (size_t i = 0; i < phasefieldLocalFiniteElement.size(); i++)
      for (size_t j = 0; j < phasefieldLocalFiniteElement.size(); j++) {
        size_t vIndex = localView.tree().child(_1).localIndex(i);
        size_t pIndex = localView.tree().child(_1).localIndex(j);
        //
        auto valueFirst = 2.0 * phaseFieldValue[i] * phaseFieldValue[j] * weight *
                          integrationElement;
        auto valueSecond = material._regularisationParameter *
                           material._regularisationParameter * weight *
                           integrationElement;
        elementMatrix[vIndex][pIndex] += valueFirst;
      }
  }
}

// Set the occupation pattern of the stiffness matrix
template <class Basis, class Matrix>
void setOccupationPattern(const Basis &basis, Matrix &matrix) {
  enum { dim = Basis::GridView::dimension };
  // MatrixIndexSets store the occupation pattern of a sparse matrix.
  // They are not particularly efficient, but simple to use.
  std::array<std::array<MatrixIndexSet, 2>, 2> nb;
  // Set sizes of the 2x2 submatrices
  for (size_t i = 0; i < 2; i++)
    for (size_t j = 0; j < 2; j++)
      nb[i][j].resize(basis.size({i}), basis.size({j}));
  // A view on the FE basis on a single element
  auto localView = basis.localView();
  // Loop over all leaf elements
  for (const auto &element : elements(basis.gridView())) {
    // Bind the local view to the current element
    localView.bind(element);
    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i = 0; i < localView.size(); i++) {
      // Global index of the i-th local degree of freedom of the current element
      auto row = localView.index(i);
      for (size_t j = 0; j < localView.size(); j++) {
        // Global index of the j-th local degree of freedom of the current
        // element
        auto col = localView.index(j);
        nb[row[0]][col[0]].add(row[1], col[1]);
      }
    }
  }
  // Give the matrix the occupation pattern we want.
  using namespace Indices;
  nb[0][0].exportIdx(matrix[_0][_0]);
  nb[0][1].exportIdx(matrix[_0][_1]);
  nb[1][0].exportIdx(matrix[_1][_0]);
  nb[1][1].exportIdx(matrix[_1][_1]);
}

template <class Matrix, class MultiIndex>
decltype(auto) matrixEntry(Matrix &matrix, const MultiIndex &row,
                           const MultiIndex &col) {
  using namespace Indices;
  if ((row[0] == 0) && (col[0] == 0))
    return matrix[_0][_0][row[1]][col[1]][row[2]][col[2]];
  if ((row[0] == 0) && (col[0] == 1))
    return matrix[_0][_1][row[1]][col[1]][row[2]][0];
  if ((row[0] == 1) && (col[0] == 0))
    return matrix[_1][_0][row[1]][col[1]][0][col[2]];
  return matrix[_1][_1][row[1]][col[1]];
}

// Assemble the Laplace stiffness matrix on the given grid view
template <class Basis, class displacementFunction, class phasefieldFunction,
          class Matrix>
void assemblePhasefieldMatrix(const Basis &basis,
                              const displacementFunction &dispFunction,
                              const phasefieldFunction &phaseFieldFunction,
                              Matrix &matrix) {
  // Set all entries to zero
  matrix = 0;
  // A view on the FE basis on a single element
  auto localView = basis.localView();
  auto localDispFunctionView = localFunction(dispFunction);
  auto localPhaseFunctionView = localFunction(phaseFieldFunction);
  // A loop over all elements of the grid
  auto material = Dune::Biw467::material<2>(0.01, 100, 200, 0.1);
  for (const auto &element : elements(basis.gridView())) {
    // Bind the local FE basis view to the current element
    localView.bind(element);
    localDispFunctionView.bind(element);
    localPhaseFunctionView.bind(element);
    // Now let’s get the element stiffness matrix
    // A dense matrix is used for the element stiffness matrix
    Dune::Matrix<double> elementMatrix;
    Dune::BlockVector<double> elementResidualVector;
    assembleElementStiffnessMatrix(localView, localDispFunctionView,
                                   localPhaseFunctionView, elementMatrix,
                                   material);
    // Add element stiffness matrix onto the global stiffness matrix
    for (size_t i = 0; i < elementMatrix.N(); i++) {
      // The global index of the i-th local degree of freedom
      // of the current element
      auto row = localView.index(i);
      for (size_t j = 0; j < elementMatrix.M(); j++) {
        // The global index of the j-th local degree of freedom
        // of the current element
        auto col = localView.index(j);
        matrixEntry(matrix, row, col) += elementMatrix[i][j];
      }
    }
    // { accumulate_global_matrix_end }
  }
}

int main(int argc, char *argv[]) {
  // Set up MPI, if available
  MPIHelper::instance(argc, argv);
  ///////////////////////////////////
  // Generate the grid
  ///////////////////////////////////
  constexpr int dim = 2;
  using Grid = YaspGrid<dim>;
  FieldVector<double, dim> upperRight = {1, 1};
  std::array<int, dim> nElements = {4, 4};
  Grid grid(upperRight, nElements);
  using GridView = typename Grid::LeafGridView;
  GridView gridView = grid.leafGridView();
  /////////////////////////////////////////////////////////
  // Choose a finite element space
  /////////////////////////////////////////////////////////
  using namespace Functions::BasisFactory;
  constexpr std::size_t p = 1; // Pressure order for Taylor-Hood
  auto dispPhase = makeBasis(gridView,
                             composite(power<dim>(                //
                                           lagrange<p>(),         //
                                           blockedInterleaved()), //
                                       lagrange<p>()              //
                                       ));
  using namespace Indices;
  auto displacementBasis = Functions::subspaceBasis(dispPhase, _0);
  auto phasefieldBasis = Functions::subspaceBasis(dispPhase, _1);
  /////////////////////////////////////////////////////////
  // Stiffness matrix and right hand side vector
  /////////////////////////////////////////////////////////
  using DisplacmentVector = BlockVector<FieldVector<double, dim>>;
  using PhasefieldVector = BlockVector<double>;
  using Vector = MultiTypeBlockVector<DisplacmentVector, PhasefieldVector>;
  using Matrix00 = BCRSMatrix<FieldMatrix<double, dim, dim>>;
  using Matrix01 = BCRSMatrix<FieldMatrix<double, dim, 1>>;
  using Matrix10 = BCRSMatrix<FieldMatrix<double, 1, dim>>;
  using Matrix11 = BCRSMatrix<double>;
  using MatrixRow0 = MultiTypeBlockVector<Matrix00, Matrix01>;
  using MatrixRow1 = MultiTypeBlockVector<Matrix10, Matrix11>;
  using Matrix = MultiTypeBlockMatrix<MatrixRow0, MatrixRow1>;
  /////////////////////////////////////////////////////////
  // Assemble the system
  /////////////////////////////////////////////////////////
  Vector rhs;
  Vector sol;
  Vector solInkrement;
  auto rhsBackend = Functions::istlVectorBackend(rhs);
  rhsBackend.resize(dispPhase);
  // initialisation with copy.
  rhs = 0;
  sol = rhs;
  solInkrement = rhs;

  // make a solution function
  auto displacementFunction =
      Functions::makeDiscreteGlobalBasisFunction<DisplacmentVector>(
          displacementBasis, sol);
  auto phasefieldFunction =
      Functions::makeDiscreteGlobalBasisFunction<PhasefieldVector>(
          phasefieldBasis, sol);

  Matrix stiffnessMatrix;
  // Set matrix size and occupation pattern
  setOccupationPattern(dispPhase, stiffnessMatrix);

  assemblePhasefieldMatrix(dispPhase, displacementFunction, phasefieldFunction,
                           stiffnessMatrix);
  /////////////////////////////////////////////////////////
  // Set Dirichlet values.
  // Only velocity components have Dirichlet boundary values
  /////////////////////////////////////////////////////////
  using VelocityBitVector = std::vector<std::array<char, dim>>;
  using PressureBitVector = std::vector<char>;
  using BitVector = TupleVector<VelocityBitVector, PressureBitVector>;
  BitVector isBoundary;
  auto isBoundaryBackend = Functions::istlVectorBackend(isBoundary);
  isBoundaryBackend.resize(dispPhase);

  for (auto &&b0i : isBoundary[_0])
    for (std::size_t j = 0; j < b0i.size(); ++j)
      b0i[j] = false;
  std::fill(isBoundary[_1].begin(), isBoundary[_1].end(), false);
  Functions::forEachBoundaryDOF(
      Functions::subspaceBasis(dispPhase, _0),
      [&](auto &&index) { isBoundaryBackend[index] = true; });

  using Coordinate = GridView::Codim<0>::Geometry::GlobalCoordinate;
  using DisplacementRange = FieldVector<double, dim>;
  auto &&g = [](Coordinate x) {
    return DisplacementRange{0.0, (x[0] < 1e-8) ? 1.0 : 0.0};
  };
  Functions::interpolate(Functions::subspaceBasis(dispPhase, _0), rhs, g,
                         isBoundary);

  ////////////////////////////////////////////
  // Modify Dirichlet rows
  ////////////////////////////////////////////

  auto localView = dispPhase.localView();
  for (const auto &element : elements(gridView)) {
    localView.bind(element);
    for (size_t i = 0; i < localView.size(); ++i) {
      auto row = localView.index(i);
      // If row corresponds to a boundary entry,
      // modify it to be an identity matrix row.
      if (isBoundaryBackend[row])
        for (size_t j = 0; j < localView.size(); ++j) {
          auto col = localView.index(j);
          matrixEntry(stiffnessMatrix, row, col) = (i == j) ? 1 : 0;
        }
    }
  }

  ////////////////////////////
  // Compute solution
  ////////////////////////////
  // Initial iterate: Start from the rhs vector,
  // that way the Dirichlet entries are already correct.
  Vector x = rhs;
  // Turn the matrix into a linear operator
  MatrixAdapter<Matrix, Vector, Vector> stiffnessOperator(stiffnessMatrix);
  // Fancy (but only) way to not have a preconditioner at all
  Richardson<Vector, Vector> preconditioner(0.8);
  // Construct the iterative solver
  RestartedGMResSolver<Vector> solver(stiffnessOperator, // Operator to invert
                                      preconditioner,
                                      // Preconditioner
                                      1e-10,
                                      // Desired residual reduction factor
                                      500,
                                      // Number of iterations between restarts,
                                      // here: no restarting
                                      500,
                                      // Maximum number of iterations
                                      2);
  // Verbosity of the solver
  // Object storing some statistics about the solving process
  InverseOperatorResult statistics;
  // Solve!
  solver.apply(x, rhs, statistics);

  ////////////////////////////////////////////////////////////////////////////
  // Make a discrete function from the FE basis and the coefficient vector
  ////////////////////////////////////////////////////////////////////////////

  using DisplacementRange = FieldVector<double, dim>;
  using PressureRange = double;
  auto velocityFunction =
      Functions::makeDiscreteGlobalBasisFunction<DisplacementRange>(
          Functions::subspaceBasis(dispPhase, _0), x);
  auto pressureFunction =
      Functions::makeDiscreteGlobalBasisFunction<PressureRange>(
          Functions::subspaceBasis(dispPhase, _1), x);

  //////////////////////////////////////////////////////////////////////////////////////////////
  // Write result to VTK file
  // We need to subsample, because the dune-grid VTKWriter cannot natively
  // display second-order functions
  //////////////////////////////////////////////////////////////////////////////////////////////

  SubsamplingVTKWriter<GridView> vtkWriter(gridView, refinementLevels(2));

  vtkWriter.addVertexData(
      velocityFunction,
      VTK::FieldInfo("displacments", VTK::FieldInfo::Type::vector, dim));
  vtkWriter.addVertexData(
      pressureFunction,
      VTK::FieldInfo("phasefield", VTK::FieldInfo::Type::scalar, 1));
  vtkWriter.write("phasefield-result");
}