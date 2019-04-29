/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include <cassert>
#include <limits>
#include <cblas.h>

#include "bitpit_IO.hpp"
#include "bitpit_LA.hpp"
#include "bitpit_operators.hpp"

#include "bitpit_private_lapacke.hpp"

#include "reconstruction.hpp"

namespace bitpit {

/*!
 * \class ReconstructionBase
 * \ingroup discretization
 *
 * \brief The ReconstructionBase class is the base calss upon which the
 * reconstruction objects are build.
 *
 * ReconstructionBase is the base calss for the reconstruction objects.
 */

/*!
 * Defualt constructor.
 */
ReconstructionBase::ReconstructionBase()
    : m_degree(0), m_dimensions(0)
{
}

/*!
 * Constructor.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] dimensions is the number of space dimensions
 */
ReconstructionBase::ReconstructionBase(uint8_t degree, uint8_t dimensions)
{
    _initialize(degree, dimensions);
}

/**
 * Exchanges the content of the reconstruction by the content the specified
 * other reconstruction.
 *
 * \param other is another reconstruction whose content is swapped with that
 * of this reconstruction
 */
void ReconstructionBase::swap(ReconstructionBase &other) noexcept
{
    _swap(other);
}

/*!
 * Initialize the reconstruction.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 */
void ReconstructionBase::initialize(uint8_t degree, uint8_t dimensions)
{
    _initialize(degree, dimensions);
}

/**
 * Internal function to exchange the content of the reconstruction by the
 * content the specified other reconstruction.
 *
 * \param other is another reconstruction whose content is swapped with that
 * of this reconstruction
 */
void ReconstructionBase::_swap(ReconstructionBase &other) noexcept
{
    std::swap(other.m_degree, m_degree);
    std::swap(other.m_dimensions, m_dimensions);
}

/*!
 * Internal function to initialize the reconstruction.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 */
void ReconstructionBase::_initialize(uint8_t degree, uint8_t dimensions)
{
    assert(degree == 0 || degree == 1 || degree == 2);
    m_degree = degree;

    assert(dimensions == 1 || dimensions == 2 || dimensions == 3);
    m_dimensions = dimensions;
}

/*!
 * Get the degree of the polynomial.
 *
 * \return The degree of the polynomial.
 */
uint8_t ReconstructionBase::getDegree() const
{
    return m_degree;
}

/*!
 * Get the number of space dimensions.
 *
 * \return The number of space dimensions.
 */
uint8_t ReconstructionBase::getDimensions() const
{
    return m_dimensions;
}

/*!
 * Evaluate the number of coefficients of the reconstruction polynomial.
 *
 * \return The number of coefficients of the reconstruction polynomial.
 */
int ReconstructionBase::countCoefficients() const
{
    return countCoefficients(m_degree);
}

/*!
 * Evaluate the number of coefficients of the reconstruction polynomial
 * up to the specified degree.
 *
 * \param degree is the degree
 * \return The number of coefficients of the reconstruction polynomial
 * up to the specified degree.
 */
int ReconstructionBase::countCoefficients(int degree) const
{
    int nCoeffs = 0;
    for (int i = 0; i <= degree; ++i) {
        nCoeffs += countDegreeCoefficients(i);
    }

    return nCoeffs;
}

/*!
 * Evaluate the number of coefficients associates with the specified degree
 * of the reconstruction polynomial.
 *
 * \param degree is the degree
 * \return The number of coefficients associates with the specified degree
 * of the reconstruction polynomial.
 */
int ReconstructionBase::countDegreeCoefficients(int degree) const
{
    switch (degree) {

    case 0:
        return 1;

    case 1:
        return m_dimensions;

    case 2:
        return (2 * m_dimensions - 1);

    default:
        return (utils::factorial(m_dimensions - 1 + degree) / utils::factorial(m_dimensions - 1) / utils::factorial(degree));

    }
}

/*!
 * Evaluates the coefficients of the equation that describes the reconstruction
 * of a point value.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction coefficients will
 * be evaluated
 * \param[out] equationCoeffs on output will contain the coefficients of the
 * point value reconstruction
 */
void ReconstructionBase::evalPointValueEquationCoeffs(int degree,
                                                      const std::array<double, 3> &origin,
                                                      const std::array<double,3> &point,
                                                      double *equationCoeffs) const
{
    // Set 0-th degree coefficients
    equationCoeffs[0] = 1.;

    // Set high degree coefficients
    if (degree >= 1) {
        int offset = 1;
        const std::array<double,3> distance = point - origin;

        // Set 1-st degree coefficients
        for (int i = 0; i < m_dimensions; ++i) {
            equationCoeffs[offset++] = distance[i];
        }

        // Set 2-nd degree coefficients
        if (degree >= 2) {
            for (int i = 0; i < m_dimensions; ++i) {
                equationCoeffs[offset++] = 0.5 * distance[i] * distance[i];
            }

            if (m_dimensions >= 2) {
                equationCoeffs[offset++] = distance[0] * distance[1];

                if (m_dimensions >= 3) {
                    equationCoeffs[offset++] = distance[0] * distance[2];
                    equationCoeffs[offset++] = distance[1] * distance[2];
                }
            }
        }

        // Check if all coefficients have been set
        assert(offset == countCoefficients());
    }
}

/*!
 * Evaluates the coefficients of the equation that describes the reconstruction
 * of the derivative in a point.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction coefficients will
 * be evaluated
 * \param[in] direction the direction of the derivative
 * \param[out] equationCoeffs on output will contain the coefficients of the
 * point derivative reconstruction
 */
void ReconstructionBase::evalPointDerivativeEquationCoeffs(int degree,
                                                           const std::array<double, 3> &origin,
                                                           const std::array<double,3> &point,
                                                           const std::array<double,3> &direction,
                                                           double *equationCoeffs) const
{
    // Set 0-th degree coefficients
    equationCoeffs[0] = 0.;

    // Set high degree coefficients
    if (degree >= 1) {
        int offset = 1;

        // Set 1-st degree coefficients
        for (int i = 0; i < m_dimensions; ++i) {
            equationCoeffs[offset++] = direction[i];
        }

        // Set 2-nd degree coefficients
        if (degree >= 2) {
            const std::array<double,3> distance = point - origin;

            for (int i = 0; i < m_dimensions; ++i) {
                equationCoeffs[offset++] = distance[i] * direction[i];
            }

            if (m_dimensions >= 2) {
                equationCoeffs[offset++] = distance[0] * direction[1] + distance[1] * direction[0];

                if (m_dimensions >= 3) {
                    equationCoeffs[offset++] = distance[0] * direction[2] + distance[2] * direction[0];
                    equationCoeffs[offset++] = distance[1] * direction[2] + distance[2] * direction[1];
                }
            }
        }

        // Check if all coefficients have been set
        assert(offset == countCoefficients());
    }
}

/*!
 * Evaluates the coefficients of the equation that describes the reconstruction
 * of a cell average.
 *
 * The method works only for ElementType::Voxel and ElementType::Pixel.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] cell is the cell
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] vertexCoords are the vertecx coordinates
 * \param[out] equationCoeffs on output will contain the coefficients of the
 * cell average reconstruction
 */
void ReconstructionBase::evalCellAverageEquationCoeffs(int degree,
                                                       const Cell &cell,
                                                       const std::array<double, 3> &origin,
                                                       const std::array<double, 3> *vertexCoords,
                                                       double *equationCoeffs) const
{
    // Check if cell type is supported
    ElementType cellType = cell.getType();

    bool cellTypeSupported = false;
    if (cellType == ElementType::PIXEL) {
        cellTypeSupported = true;
    } else if (cellType == ElementType::VOXEL) {
        cellTypeSupported = true;
    }

    if (!cellTypeSupported) {
        throw std::runtime_error("Cell type not supported.");
    }

    // Set 0-th degree coefficients
    equationCoeffs[0] = 1.;

    // Set high degree coefficients
    if (degree >= 1) {
        int offset = 1;
        const std::array<double, 3> distance = cell.evalCentroid(vertexCoords) - origin;

        // Set 1-st degree coefficients
        for (int i = 0; i < m_dimensions; ++i) {
            equationCoeffs[offset++] = distance[i];
        }

        // Set 2-nd degree coefficients
        if (degree >= 2) {
            double cellSize = cell.evalSize(vertexCoords);

            for (int i = 0; i < m_dimensions; ++i) {
                equationCoeffs[offset++] = 0.5 * (distance[i] * distance[i] + cellSize * cellSize / 12.);
            }

            if (m_dimensions >= 2) {
                equationCoeffs[offset++] = distance[0] * distance[1];

                if (m_dimensions >= 3) {
                    equationCoeffs[offset++] = distance[0] * distance[2];
                    equationCoeffs[offset++] = distance[1] * distance[2];
                }
            }
        }

        // Check if all coefficients have been set
        assert(offset == countCoefficients());
    }
}

/*!
 * Apply the limiter to the given coefficients.
 *
 * It is necessary to provide one scale factor for each polinomaial degree
 * greater than zero (a reconstruction using 2nd degree polinomia will need
 * two scale factors, one for the 1st degree coefficients and one for the
 * 2nd degree coefficients).
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] limiterScaleFactors are the scale factors of the limiter
 * \param[in,out] coeffs are the coefficients to be limited
 */
void ReconstructionBase::applyLimiter(int degree, const double *limiterScaleFactors, double *coeffs) const
{
    if (degree < 1) {
        return;
    }

    int coeffBegin = 0;
    int coeffEnd   = countDegreeCoefficients(0);
    for (int n = 1; n <= degree; ++n) {
        coeffBegin = coeffEnd;
        coeffEnd   = coeffBegin + countDegreeCoefficients(n);
        for (int k = coeffBegin; k < coeffEnd; ++k) {
            coeffs[k] *= limiterScaleFactors[n - 1];
        }
    }
}

/*!
 * \class ReconstructionAssembler
 * \ingroup discretization
 *
 * \brief The ReconstructionAssembler class allows to define reconstruction
 * polinoymial.
 */


/*!
 * Is the threshold for which a singuler value is considered zero.
 */
double ReconstructionAssembler::SVD_ZERO_THRESHOLD = 1e-9;

/*!
 * Constructor.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] dimensions is the number of space dimensions
 */
ReconstructionAssembler::ReconstructionAssembler(uint8_t degree, uint8_t dimensions)
    : ReconstructionBase(degree, dimensions)
{
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void ReconstructionAssembler::swap(ReconstructionAssembler &other) noexcept
{
    ReconstructionBase::_swap(other);

    _swap(other);
}

/**
* Internal function to exchange the content of the reconstruction by the
* content the specified other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void ReconstructionAssembler::_swap(ReconstructionAssembler &other) noexcept
{
    std::swap(other.m_constraintsOrder, m_constraintsOrder);
    std::swap(other.m_leastSquaresOrder, m_leastSquaresOrder);
    std::swap(other.m_A, m_A);
    std::swap(other.m_C, m_C);
}

/*!
 * Clear the reconstruction.
 */
void ReconstructionAssembler::clear()
{
    m_constraintsOrder.clear();
    m_leastSquaresOrder.clear();

    m_A.clear();
    m_C.clear();
}

/*!
 * Count the number of constrain-type equations added to the assembler.
 *
 * \return The number of constraint-type equations added to the assembler.
 */
int ReconstructionAssembler::countConstraints() const
{
    return m_constraintsOrder.size();
}

/*!
 * Count the number of least square-type equations added to the assembler.
 *
 * \return The number of least square-type equations added to the assembler.
 */
int ReconstructionAssembler::countLeastSquares() const
{
    return m_leastSquaresOrder.size();
}

/*!
 * Count the number of equations added to the assembler.
 *
 * \return The number of equations added to the assembler.
 */
int ReconstructionAssembler::countEquations() const
{
    int nConstraints  = countConstraints();
    int nLeastSquares = countLeastSquares();
    int nEquations    = nConstraints + nLeastSquares;

    return nEquations;
}

/*!
 * Add a point value equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction coefficients will
 * be evaluated
 */
void ReconstructionAssembler::addPointValueEquation(ReconstructionType type,
                                                    const std::array<double, 3> &origin,
                                                    const std::array<double, 3> &point)
{
    double *equationCoeffsStorage = _addEquation(type);
    evalPointValueEquationCoeffs(getDegree(), origin, point, equationCoeffsStorage);
}

/*!
 * Add a point derivative equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param origin is the point chosen as origin of the reconstruction
 * \param point is the point were the reconstruction coefficients will
 * be evaluated
 * \param direction the direction of the derivative
 */
void ReconstructionAssembler::addPointDerivativeEquation(ReconstructionType type,
                                                         const std::array<double, 3> &origin,
                                                         const std::array<double, 3> &point,
                                                         const std::array<double, 3> &direction)
{
    double *equationCoeffsStorage = _addEquation(type);
    evalPointDerivativeEquationCoeffs(getDegree(), origin, point, direction, equationCoeffsStorage);
}

/*!
 * Add a cell average equation.
 *
 * \param type is the type of reconstruction associated to the equation
 * \param cell is the cell
 * \param origin is the point chosen as origin of the reconstruction
 * \param vertexCoords are the vertecx coordinates
 */
void ReconstructionAssembler::addCellAverageEquation(ReconstructionType type,
                                                     const Cell &cell,
                                                     const std::array<double, 3> &origin,
                                                     const std::array<double, 3> *vertexCoords)
{
    double *equationCoeffsStorage = _addEquation(type);
    evalCellAverageEquationCoeffs(getDegree(), cell, origin, vertexCoords, equationCoeffsStorage);
}

/*!
 * Internal function to add an equation.
 *
 * \param type is the type of reconstruction associated to the equation
 */
double * ReconstructionAssembler::_addEquation(ReconstructionType type)
{
    // Update equation counters
    int nEquations = countEquations();
    switch (type) {

    case TYPE_CONSTRAINT:
        m_constraintsOrder.emplace_back(nEquations);
        break;

    case TYPE_LEAST_SQUARE:
        m_leastSquaresOrder.emplace_back(nEquations);
        break;

    }

    // Prepare storage for equation coefficients
    int nCoeffs = countCoefficients();

    double *equationCoeffsStorage = nullptr;
    switch (type) {

    case TYPE_CONSTRAINT:
        m_C.emplace_back(nCoeffs);
        equationCoeffsStorage = m_C.back().data();
        break;

    case TYPE_LEAST_SQUARE:
        m_A.emplace_back(nCoeffs);
        equationCoeffsStorage = m_A.back().data();
        break;

    }

    return equationCoeffsStorage;
}

/*!
 * Computes the weights to be used in degree to calculate the coefficients
 * of the polynomial. The coefficients are such that the conditions decoded
 * in equations are enforced.
 *
 * \param[out] polynomial on output will contain the reconstruciton polynomial
 */
void ReconstructionAssembler::assembly(ReconstructionPolynomial *polynomial)
{
    // Get the number of equations
    int nConstraints  = countConstraints();
    int nLeastSquares = countLeastSquares();
    int nEquations    = countEquations();

    // Get the number of polynomial coefficients
    int nCoeffs = countCoefficients();

    // Initialize polynomial
    int degree     = getDegree();
    int dimensions = getDimensions();

    polynomial->initialize(degree, dimensions, nEquations);

    // The linear-constrained are introduced in the least-squares problem
    // through Lagrange multipliers. The resulting linear system is:
    //
    // | A^t A  C^t | |x     | = |A^t b|
    // | C      0   | |lambda|   |d    |
    //
    // with A and C the least-squares and costraints equations respectively,
    // and b and d their corresponding RHSs. x are the coefficients of
    // the polynomial and lambda the lagrange multipliers.
    //
    // This system is denoted by S:
    //
    //     |  x   |   |A^t  0| |b|
    // |S| |      | = |      | | |
    //     |lambda|   |0    I| |d|
    // 
    // The matrices S and S^-1 are symmetric and only the upper portions are
    // computed
    int nUnknowns = nCoeffs + nConstraints;

    std::vector<double> S(nUnknowns * nUnknowns, 0.);
    for (int i = 0; i < nCoeffs; ++i) {
        for (int j = i; j < nCoeffs; ++j) {
            // Compute A^t A on the fly
            double ATA_ij = 0.;
            for (int k = 0; k < nLeastSquares; ++k) {
                ATA_ij += m_A[k][i] * m_A[k][j];
            }

            int l = linearalgebra::linearIndexColMajor(i, j, nUnknowns, nUnknowns);
            S[l] = ATA_ij;

            int m = linearalgebra::linearIndexColMajor(j, i, nUnknowns, nUnknowns);
            S[m] = ATA_ij;
        }

        for (int j = nCoeffs; j < nUnknowns; ++j) {
            int l = linearalgebra::linearIndexColMajor(i, j, nUnknowns, nUnknowns);
            S[l] = m_C[j - nCoeffs][i];

            int m = linearalgebra::linearIndexColMajor(j, i, nUnknowns, nUnknowns);
            S[m] = S[l];
        }
    }

    // Compute inverse S matrix
    // Since S may me be rank-deficit (eg if not enough neighbours are available)
    // the pseudo-inverse is used. This corresponds of computing the least-norm
    // solution of the problem.
    computePseudoInverse(nUnknowns, nUnknowns, SVD_ZERO_THRESHOLD, S.data());

    // Weights needed to evaluate the polynomial coefficients come from the
    // following equation:
    //
    // |  x   |        |A^t  0| |b|          |b|
    // |      | = S^-1 |      | | | = S^-1 Q | |
    // |lambda|        |0    I| |d|          |d|
    //
    // Since we are interested only in x (the polynomial coefficients) only
    // the first nCoeffs rows of the matrix S^-1 Q are computed. Those values
    // are the polynomial weights.
    //
    // Weigths are stored according the order in which the equations have been
    // added.
    double *weights = polynomial->getPolynomialWeights();
    for (int i = 0; i < nCoeffs; ++i) {
        for (int j = 0; j < nEquations; ++j) {
            double value = 0;
            for (int k = 0; k < nUnknowns; ++k) {
                int l = linearalgebra::linearIndexColMajorSymmetric(i, k, nUnknowns, nUnknowns, 'U');
                if (k < nCoeffs && j < nLeastSquares) {
                    value += S[l] * m_A[j][k];
                } else if ((k - nCoeffs) == (j - nLeastSquares)) {
                    value += S[l];
                }
            }

            int equation;
            if (j < nLeastSquares) {
                equation = m_leastSquaresOrder[j];
            } else {
                equation = m_constraintsOrder[j - nLeastSquares];
            }

            int weightLineraIndex = linearalgebra::linearIndexColMajor(i, equation, nCoeffs, nEquations);
            weights[weightLineraIndex] = value;
        }
    }
}
/*!
 * Computes the pseudo inverse of a matrix using a singular value decomposition
 *
 * See "Solving Ill-Conditioned And Singular Linear Systems: A Tutorial On
 * Regularization", by Arnold Neumaier (see https://www.mat.univie.ac.at/~neum/ms/regtutorial.pdf).
 *
 * \param[in] m number of columns
 * \param[in] n number of rows
 * \param[in] zeroThreshold is the threshold below which a singuler value is
 * considered zero
 * \param[in,out] A on input matrix in coumn-major ordering, on output its
 * pseudo-inverse
 */
void ReconstructionAssembler::computePseudoInverse(int m, int n, double zeroThreshold, double *A)
{
    int k = std::min(m,n);

    std::vector<double> sigma(k, 0.);
    std::vector<double> U(m * k, 0.);
    std::vector<double> Vt(k * n, 0.);

    // Compute SVD and check for error
    //
    // A = U * Sigma * Vt (Equation 21)
    std::vector<double> supdB(std::min(m,n) - 1, 0.);
    int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, A, m, sigma.data(),
                              U.data(), m, Vt.data(), k, supdB.data());

    if (info > 0) {
        log::cout() << "SVD failed in ReconstructionAssembler::computePseudoInverse()" <<std::endl;
        exit(1);
    }

    // Inv(A) = V * Sigma^ + *U^T (Equation 22)
    //
    // u = sigma^ + *U
    // and is stored in U
    for (int i = 0; i < k; ++i) {
       double sigma_plus = (sigma[i] > zeroThreshold) ? (1. / sigma[i]) : sigma[i];
       cblas_dscal(m, sigma_plus, &U[i*m], 1);
    }

    // Inv(A) = (Vt)^T * u^T
    cblas_dgemm(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans, CBLAS_TRANSPOSE::CblasTrans,
                n, m, k, 1., Vt.data(), k, U.data(), m, 0., A, n);
}

/*!
 * \class ReconstructionPolynomial
 * \ingroup discretization
 *
 * \brief The ReconstructionPolynomial class allows to apply a reconstruction
 * polinoymial previously assembled.
 */


const int ReconstructionPolynomial::MAX_STACK_WORKSPACE_SIZE = 10;

/*!
 * Constructor.
 */
ReconstructionPolynomial::ReconstructionPolynomial()
    : ReconstructionBase()
{
    clear();
}

/*!
 * Constructor.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param nEquations is the number of equations that defines the reconstruction
 */
ReconstructionPolynomial::ReconstructionPolynomial(uint8_t degree, uint8_t dimensions, int nEquations)
    : ReconstructionBase(degree, dimensions)
{
    _initialize(nEquations);
}

/*!
    Copy constructor

    \param other is another reconstruction whose content is copied in this
    reconstruction
*/
ReconstructionPolynomial::ReconstructionPolynomial(const ReconstructionPolynomial &other)
    : ReconstructionBase(other.getDegree(), other.getDimensions())
{
    _initialize(other.m_nEquations);

    if (m_nEquations > 0) {
        int nCoeffs  = countCoefficients();
        int nWeights = nCoeffs * m_nEquations;

        std::copy(other.m_weights.get(), other.m_weights.get() + nWeights, m_weights.get());
    }
}

/*!
    Copy-assignament operator.

    \param other is another reconstruction whose content is copied in this
    reconstruction
*/
ReconstructionPolynomial & ReconstructionPolynomial::operator=(const ReconstructionPolynomial &other)
{
    ReconstructionPolynomial tmp(other);
    swap(tmp);

    return *this;
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void ReconstructionPolynomial::swap(ReconstructionPolynomial &other) noexcept
{
    ReconstructionBase::_swap(other);

    _swap(other);
}

/*!
 * Initialize the polynomial.
 *
 * \param degree is the degree of the polynomial
 * \param dimensions is the number of space dimensions
 * \param nEquations is the number of equations that defines the reconstruction
 */
void ReconstructionPolynomial::initialize(uint8_t degree, uint8_t dimensions, int nEquations)
{
    ReconstructionBase::initialize(degree, dimensions);

    _initialize(nEquations);
}

/**
 * Internal function to exchange the content of the reconstruction by the
 * content the specified other reconstruction.
 *
 * \param other is another reconstruction whose content is swapped with that
 * of this reconstruction
 */
void ReconstructionPolynomial::_swap(ReconstructionPolynomial &other) noexcept
{
    std::swap(other.m_nEquations, m_nEquations);
    std::swap(other.m_weights, m_weights);
}

/*!
 * Internal function to initialize the polynomial.
 *
 * \param nEquations is the number of equations that defines the reconstruction
 */
void ReconstructionPolynomial::_initialize(int nEquations)
{
    m_nEquations = nEquations;

    int nCoeffs = countCoefficients();
    m_weights = std::unique_ptr<double[]>(new double[nCoeffs * nEquations]);
}

/*!
 * Clear the reconstruction.
 */
void ReconstructionPolynomial::clear()
{
    m_nEquations = 0;

    m_weights.reset();
}

/*!
 * Get the number of equations associated with the reconstruction.
 *
 * \return The number of equations associated with the reconstruction.
 */
int ReconstructionPolynomial::getEquationCount() const
{
    return m_nEquations;
}

/*!
 * Returns a constant pointer to the internal weight that are used to evaluate
 * the polynomial coefficients.
 *
 * \result A constant pointer to the internal weight that are used to evaluate
 * the polynomial coefficients.
 */
const double * ReconstructionPolynomial::getPolynomialWeights() const
{
    return m_weights.get();
}

/*!
 * Returns a pointer to the internal weight that are used to evaluate the
 * polynomial coefficients.
 *
 * \result A pointer to the internal weight that are used to evaluate the
 * polynomial coefficients.
 */
double * ReconstructionPolynomial::getPolynomialWeights()
{
    return m_weights.get();
}

/*!
 * Returns the internal weight that are used to evaluate the specified
 * polynomial coefficient.
 *
 * \param[in] index is the index of the polynomial coefficient
 * \return The internal weight that are used to evaluate the specified
 * polynomial coefficient.
 */
std::vector<double> ReconstructionPolynomial::getPolynomialWeights(uint8_t index) const
{
    int nCoeffs = countCoefficients();
    if (index >= nCoeffs) {
        throw std::runtime_error("Polynomial coefficient index out of range.");
    }

    int nEquations = getEquationCount();
    std::vector<double> ceofficientWeights(nEquations);
    for (int j = 0; j < nEquations; ++j) {
        ceofficientWeights[j] = m_weights[linearalgebra::linearIndexColMajor(index, j, nCoeffs, nEquations)];
    }

    return ceofficientWeights;
}

/*!
 * Computes the coefficients of the reconstruction polynomial for the specified
 * values.
 *
 * The values need to be passed in the same order as the equations.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] values are the known terms of the linear equations
 * \param[out] coeffs on output will contain the coefficients of the
 * reconstruction polynomial.
 */
void ReconstructionPolynomial::computePolynomialCoeffs( const double *values, double *coeffs) const
{
    computePolynomialCoeffs(getDegree(), values, coeffs);
}

/*!
 * Computes the coefficients of the reconstruction polynomial with the
 * given degree for the specified values.
 *
 * The degree of the polynomial has to be less than or equal to the degree
 * used in the assembly step.
 *
 * The values need to be passed in the same degree as the equations.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] values are the known terms of the linear equations
 * \param[out] coeffs on output will contain the coefficients of the
 * reconstruction polynomial.
 */
void ReconstructionPolynomial::computePolynomialCoeffs(int degree, const double *values, double *coeffs) const
{
    assert(degree <= getDegree());

    int nEquations = getEquationCount();
    int nCoeffs = countCoefficients(degree);
    int nTotalCoeffs = countCoefficients();

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasNoTrans,
                nCoeffs, nEquations, 1., getPolynomialWeights(), nTotalCoeffs,
                values, 1, 0, coeffs, 1);
}

/*!
 * Computes the coefficients of the reconstruction polynomial with the
 * given degree for the specified values.
 *
 * The degree of the polynomial has to be less than or equal to the degree
 * used in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param[in] values are the known terms of the linear equations
 * \param[in] limiterScaleFactors are the scale factors of the limiter, if
 * a null pointer to scale factor is passed, no limiter will be applied
 * \param[out] coeffs on output will contain the coefficients of the
 * reconstruction polynomial.
 */
void ReconstructionPolynomial::computeLimitedPolynomialCoeffs(const double *values, const double *limiterScaleFactors, double *coeffs) const
{
    computeLimitedPolynomialCoeffs(getDegree(), values, limiterScaleFactors, coeffs);
}

/*!
 * Computes the coefficients of the reconstruction polynomial for the specified
 * values.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] values are the known terms of the linear equations
 * \param[in] limiterScaleFactors are the scale factors of the limiter, if
 * a null pointer to scale factor is passed, no limiter will be applied
 * \param[out] coeffs on output will contain the coefficients of the
 * reconstruction polynomial.
 */
void ReconstructionPolynomial::computeLimitedPolynomialCoeffs(int degree, const double *values, const double *limiterScaleFactors, double *coeffs) const
{
    computePolynomialCoeffs(degree, values, coeffs);

    if (limiterScaleFactors) {
        applyLimiter(degree, limiterScaleFactors, coeffs);
    }
}

/*!
 * Computes the weights of the reconstructed value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[out] weights on output will contain the weights of the reconstruction
 */
void ReconstructionPolynomial::computePointValueWeights(const std::array<double, 3> &origin, const std::array<double,3> &point, double *valueWeights) const
{
    computeLimitedPointValueWeights(getDegree(), origin, point, nullptr, valueWeights);
}

/*!
 * Computes the weights of the reconstructed value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[out] weights on output will contain the weights of the reconstruction
 */
void ReconstructionPolynomial::computePointValueWeights(int degree, const std::array<double, 3> &origin, const std::array<double,3> &point, double *valueWeights) const
{
    computeLimitedPointValueWeights(degree, origin, point, nullptr, valueWeights);
}

/*!
 * Computes the weights of the reconstructed value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[in] limiterScaleFactors are the scale factors of the limiter
 * \param[out] weights on output will contain the weights of the reconstruction
 */
void ReconstructionPolynomial::computeLimitedPointValueWeights(const std::array<double, 3> &origin, const std::array<double,3> &point,
                                                               const double *limiterScaleFactors, double *valueWeights) const
{
    computeLimitedPointValueWeights(getDegree(), origin, point, limiterScaleFactors, valueWeights);
}

/*!
 * Computes the weights of the reconstructed value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the value at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * If a valid pointer to the limiter scale factor is passed, a limited
 * reconstruction will be used. It is necessary to provide one scale factor
 * for each polinomaial degree greater than zero (a reconstruction using
 * 2nd degree polinomia will need two scale factors, one for the 1st degree
 * coefficients and one for the 2nd degree coefficients).
 *
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[in] limiterScaleFactors are the scale factors of the limiter
 * \param[out] weights on output will contain the weights of the reconstruction
 */
void ReconstructionPolynomial::computeLimitedPointValueWeights(int degree, const std::array<double, 3> &origin, const std::array<double,3> &point,
                                                               const double *limiterScaleFactors, double *valueWeights) const
{
    assert(degree <= getDegree());

    int nEquations = getEquationCount();
    int nCoeffs = countCoefficients(degree);
    int nTotalCoeffs = countCoefficients();

    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);

    evalPointValueEquationCoeffs(degree, origin, point, csi);
    if (limiterScaleFactors) {
        applyLimiter(degree, limiterScaleFactors, csi);
    }

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans,
                nCoeffs, nEquations, 1., getPolynomialWeights(), nTotalCoeffs,
                csi, 1, 0, valueWeights, 1);
}


/*!
 * Computes the weights of the directional derivative value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[in] direction the direction of the derivative
 * \param[out] weights on output will contain the weights of the derivative
 */
void ReconstructionPolynomial::computePointDerivativeWeights(const std::array<double, 3> &origin, const std::array<double,3> &point, const std::array<double,3> &direction, double *derivativeWeights) const
{
    computePointDerivativeWeights(getDegree(), origin, point, direction, derivativeWeights);
}

/*!
 * Computes the weights of the directional derivative value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[in] direction the direction of the derivative
 * \param[out] weights on output will contain the weights of the derivative
 */
void ReconstructionPolynomial::computePointDerivativeWeights(int degree, const std::array<double, 3> &origin, const std::array<double,3> &point, const std::array<double,3> &direction, double *derivativeWeights) const
{
    assert(degree <= getDegree());

    int nEquations = getEquationCount();
    int nCoeffs = countCoefficients(degree);
    int nTotalCoeffs = countCoefficients(degree);

    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);

    evalPointDerivativeEquationCoeffs(degree, origin, point, direction, csi);

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans,
                nCoeffs, nEquations, 1., getPolynomialWeights(), nTotalCoeffs,
                csi, 1, 0, derivativeWeights, 1);
}

/*!
 * Computes the weights of the gradient value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[out] weights on output will contain the weights of the gradient
 */
void ReconstructionPolynomial::computePointGradientWeights(const std::array<double, 3> &origin, const std::array<double,3> &point, std::array<double, 3> *gradientWeights) const
{
    computePointGradientWeights(getDegree(), origin, point, gradientWeights);
}

/*!
 * Computes the weights of the gradient value in a given point.
 *
 * In other words, when multiplying the values associated with the support
 * with their corrisponding weights, the derivative at point is retrieved.
 *
 * Weigths are computed according the order in which the equations have been
 * added in the assembly step.
 *
 * The degree of the reconstruction has to be less than or equal to the
 * degree used in the assembly step.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] origin is the point chosen as origin of the reconstruction
 * \param[in] point is the point were the reconstruction will be evaluated
 * \param[out] weights on output will contain the weights of the gradient
 */
void ReconstructionPolynomial::computePointGradientWeights(int degree, const std::array<double, 3> &origin, const std::array<double,3> &point, std::array<double, 3> *gradientWeights) const
{
    assert(degree <= getDegree());

    int nEquations = getEquationCount();
    int nCoeffs = countCoefficients(degree);
    int nTotalCoeffs = countCoefficients();

    BITPIT_CREATE_WORKSPACE(csi, double, nCoeffs, MAX_STACK_WORKSPACE_SIZE);
    BITPIT_CREATE_WORKSPACE(directionalStencilWeights, double, nEquations, MAX_STACK_WORKSPACE_SIZE);

    const double *polynomialWeights = getPolynomialWeights();

    // Evalaute weights in x-direction
    evalPointDerivativeEquationCoeffs(degree, origin, point, {{1., 0., 0.}}, csi);

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans,
                nCoeffs, nEquations, 1., polynomialWeights, nTotalCoeffs,
                csi, 1, 0, directionalStencilWeights, 1);

    for (int i = 0; i < nEquations; ++i) {
        gradientWeights[i][0] = directionalStencilWeights[i];
    }

    // Evalaute weights in y-direction
    evalPointDerivativeEquationCoeffs(degree, origin, point, {{0., 1., 0.}}, csi);

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans,
                nCoeffs, nEquations, 1., polynomialWeights, nTotalCoeffs,
                csi, 1, 0, directionalStencilWeights, 1);

    for (int i = 0; i < nEquations; ++i) {
        gradientWeights[i][1] = directionalStencilWeights[i];
    }

    // Evalaute weights in z-direction
    evalPointDerivativeEquationCoeffs(degree, origin, point, {{0., 0., 1.}}, csi);

    cblas_dgemv(CBLAS_ORDER::CblasColMajor, CBLAS_TRANSPOSE::CblasTrans,
                nCoeffs, nEquations, 1., polynomialWeights, nTotalCoeffs,
                csi, 1, 0, directionalStencilWeights, 1);

    for (int i = 0; i < nEquations; ++i) {
        gradientWeights[i][2] = directionalStencilWeights[i];
    }
}

/*!
 * Displays reconstruction polynomial information to an output stream.
 *
 * \param[in] out is the output stream
 * \param[in] tolerance is the tolerance below which the weight will not be
 * displayed
 */
void ReconstructionPolynomial::display(std::ostream &out, double tolerance) const
{
    int nCoeffs = countCoefficients();
    for (int index = 0; index < nCoeffs; ++index) {
        out << index << " ";

        std::vector<double> coeffcientWeights = getPolynomialWeights(index);
        for (int j = 0; j < m_nEquations; ++j) {
            if (std::abs(coeffcientWeights[j]) < tolerance) {
                continue;
            }

            out << "(" << j << "," << coeffcientWeights[j] << ") ";
        }

        out << std::endl;
    }
}

/*!
 * \class Reconstruction
 * \ingroup discretization
 *
 * \brief The Reconstruction class allows to build and apply a polynomial
 * reconstructions.
 */

/*!
 * Default constructor.
 *
 * \param[in] degree is the degree of the polynomial
 * \param[in] dimensions is the number of space dimensions
 * \param nEquations is the number of equations that defines the reconstruction
 */
Reconstruction::Reconstruction(uint8_t degree, uint8_t dimensions)
{
    ReconstructionBase::_initialize(degree, dimensions);
}

/**
* Exchanges the content of the reconstruction by the content the specified
* other reconstruction.
*
* \param other is another reconstruction whose content is swapped with that
* of this reconstruction
*/
void Reconstruction::swap(Reconstruction &other) noexcept
{
    ReconstructionBase::_swap(other);
    ReconstructionPolynomial::_swap(other);
    ReconstructionAssembler::_swap(other);
}

/*!
 * Clear the reconstruction.
 */
void Reconstruction::clear()
{
    ReconstructionAssembler::clear();
    ReconstructionPolynomial::clear();
}

/*!
 * Computes the weights to be used in order to calculate the coefficients
 * of the polynomial. The coefficients are such that the conditions decoded
 * in equations are enforced.
 */
void Reconstruction::assembly()
{
    ReconstructionAssembler::assembly(this);
}

}
