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

#ifndef __BTPIT_RECONSTRUCTION_HPP__
#define __BTPIT_RECONSTRUCTION_HPP__

#include <array>
#include <iostream>
#include <memory>
#include <vector>

#include "bitpit_patchkernel.hpp"

namespace bitpit {

class ReconstructionBase {

public:
    uint8_t getDegree() const;
    uint8_t getDimensions() const;

    int countCoefficients() const;
    int countCoefficients(int degree) const;
    int countDegreeCoefficients(int degree) const;

protected:
    ReconstructionBase();
    ReconstructionBase(uint8_t degree, uint8_t dimensions);

    void swap(ReconstructionBase &other) noexcept;

    void initialize(uint8_t degree, uint8_t dimensions);

    void evalPointValueEquationCoeffs(int degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, double *equationCoeffs) const;
    void evalPointDerivativeEquationCoeffs(int degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *equationCoeffs) const;
    void evalCellAverageEquationCoeffs(int degree, const Cell &cell, const std::array<double, 3> &origin, const std::array<double, 3> *vertexCoords, double *equationCoeffs) const;

    void applyLimiter(int degree, const double *limiterScaleFactors, double *coeffs) const;

protected:
    void _swap(ReconstructionBase &other) noexcept;
    void _initialize(uint8_t degree, uint8_t dimensions);

private:
    uint8_t m_degree;
    uint8_t m_dimensions;

};

class ReconstructionPolynomial : public virtual ReconstructionBase {

public:
    ReconstructionPolynomial();
    ReconstructionPolynomial(uint8_t degree, uint8_t dimensions, int nEquations);

    ReconstructionPolynomial(const ReconstructionPolynomial &other);
    ReconstructionPolynomial(ReconstructionPolynomial &&other) = default;
    ReconstructionPolynomial & operator = (const ReconstructionPolynomial &other);
    ReconstructionPolynomial & operator=(ReconstructionPolynomial &&other) = default;

    void swap(ReconstructionPolynomial &other) noexcept;

    void initialize(uint8_t degree, uint8_t dimensions, int nEquations);
    void clear();

    int getEquationCount() const;

    const double * getPolynomialWeights() const;
    double * getPolynomialWeights();
    std::vector<double> getPolynomialWeights(uint8_t index) const;

    void computePolynomialCoeffs(const double *values, double *coeffs) const;
    void computePolynomialCoeffs(int degree, const double *values, double *coeffs) const;
    void computeLimitedPolynomialCoeffs(const double *values, const double *limiterScaleFactors, double *coeffs) const;
    void computeLimitedPolynomialCoeffs(int degree, const double *values, const double *limiterScaleFactors, double *coeffs) const;

    void computePointValueWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, double *valueWeights) const;
    void computePointValueWeights(int degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, double *valueWeights) const;
    void computeLimitedPointValueWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, const double *limiterScaleFactors, double *valueWeights) const;
    void computeLimitedPointValueWeights(int degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const double *limiterScaleFactors, double *valueWeights) const;

    void computePointDerivativeWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *derivativeWeights) const;
    void computePointDerivativeWeights(int degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction, double *derivativeWeights) const;

    void computePointGradientWeights(const std::array<double, 3> &origin, const std::array<double, 3> &point, std::array<double, 3> *gradientWeights) const;
    void computePointGradientWeights(int degree, const std::array<double, 3> &origin, const std::array<double, 3> &point, std::array<double, 3> *gradientWeights) const;

    void display(std::ostream &out, double tolerance = 1.e-10) const;

protected:
    void _swap(ReconstructionPolynomial &other) noexcept;
    void _initialize(int nEquations);

private:
    static const int MAX_STACK_WORKSPACE_SIZE;

    int m_nEquations;

    std::unique_ptr<double[]> m_weights;

};

class ReconstructionAssembler : public virtual ReconstructionBase {

public:
    enum ReconstructionType {
        TYPE_CONSTRAINT,
        TYPE_LEAST_SQUARE
    };

    ReconstructionAssembler(uint8_t degree, uint8_t dimensions);

    void swap(ReconstructionAssembler &other) noexcept;

    void clear();

    int countConstraints() const;
    int countLeastSquares() const;
    int countEquations() const;

    void addPointValueEquation(ReconstructionType type, const std::array<double, 3> &origin, const std::array<double, 3> &point);
    void addPointDerivativeEquation(ReconstructionType type, const std::array<double, 3> &origin, const std::array<double, 3> &point, const std::array<double, 3> &direction);
    void addCellAverageEquation(ReconstructionType type, const Cell &cell, const std::array<double, 3> &origin, const std::array<double, 3> *vertexCoords);

    void assembly(ReconstructionPolynomial *polynomial);

protected:
    ReconstructionAssembler() = default;

    void _swap(ReconstructionAssembler &other) noexcept;

private:
    static double SVD_ZERO_THRESHOLD;

    std::vector<int> m_constraintsOrder;
    std::vector<int> m_leastSquaresOrder;

    std::vector<std::vector<double>> m_A;
    std::vector<std::vector<double>> m_C;

    double * _addEquation(ReconstructionType type);

    void computePseudoInverse(int m, int n, double tolerance, double *A);

};

class Reconstruction : public ReconstructionPolynomial, public ReconstructionAssembler {

public:
    Reconstruction(uint8_t degree, uint8_t dimensions);

    void swap(Reconstruction &other) noexcept;

    void clear();

    void assembly();

};

}

#endif
