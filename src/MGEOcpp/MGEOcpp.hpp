/**
 * @file MGEOcpp.hpp
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-15
 *
 * @copyright Copyright (c) 2014, Ronan Arraes Jardim Chagas. All rights
 * reserved.
 * @license This project is released under the BSD 3-Clause License (see LICENSE
 * file).
 *
 * @brief MGEO algorithm for C++.
 *
 */


#ifndef MGEOCPP_HPP
#define MGEOCPP_HPP

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <list>
#include <iostream>
#include <iomanip>
#include <string>
#include <random>

namespace MGEOcpp {

/**
 * @brief Structure that defines the limits of the design variables.
 *
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<typename Scalar>
struct sDesignVariable
{
    /// Number of bits.
    unsigned int bits;
    /// Minimum allowed value for the design variable.
    Scalar min;
    /// Maximum allowed value for the design variable.
    Scalar max;
    /// Full scale of the variable.
    uint64_t fullScale;
    /// Index in the string.
    uint64_t index;
    /// Name of the variable.
    std::string name;
};

/**
 * @brief Structure that defines a point in the Pareto frontier.
 *
 * @tparam N Number of design variables.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<unsigned int N, unsigned int nf, typename Scalar>
struct sParetoPoint
{
    /// Design variables.
    Scalar vars[N];
    /// Objective funtions. 
    Scalar f[nf];
};

/**
 * @brief Class the implements the MGEO algorithm.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<unsigned int N, 
         unsigned int nb,
         unsigned int nf, 
         typename Scalar = double>
class MGEO
{
public:
    MGEO(double tau, int ngenMax, int runMax, int rng_seed);
    MGEO(double tau, int ngenMax, int runMax);
    ~MGEO();

    /**************************************************************************
                              Getters and Setters
    ***************************************************************************/

    /// @brief Get tau.
    double getTau() const {
        return tau_;
    }

    /// @brief Get the maximum number of independent runs.
    int getRunMax() const {
        return runMax_;
    }

    /// @brief Get the maximum number of generations.
    int getNGenMax() const {
        return ngenMax_;
    }

    /// @brief Get the string.
    uint8_t* getString() const {
        return string_;
    }

    /// @brief Set tau.
    void setTau(double tau) {
        tau_ = tau;
    }

    /// @brief Set the epsilon to compare objective functions.
    void setEps(double eps) {
        mgeoEps_ = eps;
    }

    /// @brief Set the maximum number of generations.
    void setNGenMax(int ngenMax) {
        ngenMax_ = ngenMax;
    }

    /// @brief Set the maximum number of independent runs.
    void setRunMax(int runMax) {
        runMax_ = runMax;
    }
       
    bool confDesignVars(std::initializer_list<unsigned int> bits, 
                        std::initializer_list<Scalar> min, 
                        std::initializer_list<Scalar> max,
                        std::initializer_list<std::string> varNames = {});

    bool confDesignVars(unsigned int bits,
                        Scalar min,
                        Scalar max);

    /**************************************************************************
                                   Functions
    ***************************************************************************/
    
    void initializeString();
    bool callObjectiveFunctions(std::bitset<nb> string, 
                                Scalar *vars, 
                                Scalar *f);
    bool checkDominance(sParetoPoint<N, nf, Scalar> p);
    void printParetoFrontier(std::ostream& outStream = std::cout,
                             unsigned int fieldWidth = 20,
                             unsigned int precision  = 7) const;
    bool run();
    bool sortParetoFrontier(int fobj = 1);
    void stringToScalar(std::bitset<nb> string, Scalar* vars);

    /// @brief Pointer to the objective function.
    bool (*objectiveFunctions)(Scalar* vars, Scalar* f);

    /// List to store the Pareto frontier.
    std::list<sParetoPoint<N, nf, Scalar>> paretoFrontier;

private:
    /**************************************************************************
                            Configuration parameters
    ***************************************************************************/

    /// Parameter to set the search determinism.
    double tau_;

    /// Maximum number of generations.
    int ngenMax_;

    /// Maximum number of independent runs (reinitializations).
    int runMax_;

    /**************************************************************************
                                   Constants
    ***************************************************************************/

    /**************************************************************************
                                   Variables
    ***************************************************************************/

    /// Variable to store if the design variables were already configured.
    bool designVarsConfigured_;

    /// Total number of bits used by the design variables.
    int designVarsNb_;

    /// Structure to store the configuration for each design variable.
    sDesignVariable<Scalar> designVars[N];

    /// Epsilon to compare objective functions (default = 1E-10).
    double mgeoEps_;

    /// String.
    std::bitset<nb> string_;

    /**************************************************************************
                                 Random Numbers
    ***************************************************************************/
    
    /// Random number generator.
    std::mt19937 rng_;

    /// Uniform distribution, integers, 0 or 1.
    std::uniform_int_distribution<uint8_t> rand_bit_{
        std::uniform_int_distribution<uint8_t>(0,1)
            };

    /// Uniform distribution, integers, 0 to nf-1.
    std::uniform_int_distribution<int> rand_f_{
        std::uniform_int_distribution<int>(0, nf-1)
            };

    /// Uniform distribution, double, 0 to 1.
    std::uniform_real_distribution<double> rand_{
        std::uniform_real_distribution<double>(0, 1)
            };
};

#include "MGEOcppMain.hpp"

}

#endif // MGEOCPP_HPP
 
