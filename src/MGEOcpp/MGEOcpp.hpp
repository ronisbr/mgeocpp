/**
 * @file MGEOcpp.hpp
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-15
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
    /// Minimum allowed value for the design variable.
    Scalar min;
    /// Maximum allowed value for the design variable.
    Scalar max;
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
 * @tparam nb Number of bits per design variable.
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
    MGEO(double tau, int nfobMax, int runMax, int rng_seed);
    MGEO(double tau, int nfobMax, int runMax);
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

    /// @brief Get the maximum number of evaluations of the objective functions.
    int getNfobMax() const {
        return nfobMax_;
    }

    /// @brief Get the string.
    uint8_t* getString() const {
        return string;
    }

    /// @brief Set tau.
    void setTau(double tau) {
        tau_ = tau;
    }

    /// @brief Set the maximum number of evaluations of the objective functions.
    void setNfobMax(int nfobMax) {
        nfobMax_ = nfobMax;
    }

    /// @brief Set the maximum number of independent runs.
    void setRunMax(int runMax) {
        runMax_ = runMax;
    }
       
    bool setDesignVarsLimits(std::initializer_list<Scalar> min, 
                             std::initializer_list<Scalar> max);

    /**************************************************************************
                                   Functions
    ***************************************************************************/
    
    void initializeString();
    bool callObjectiveFunctions(Scalar *vars, Scalar *f);
    bool checkDominance(sParetoPoint<N, nf, Scalar> p);
    bool run();
    void printParetoFrontier(std::ostream& outStream = std::cout) const;
    void stringToScalar(Scalar* vars);

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

    /// Maximum number of evaluations of the objective function.
    int nfobMax_;

    /// Maximum number of independent runs (reinitializations).
    int runMax_;

    /**************************************************************************
                                   Constants
    ***************************************************************************/

    /// Full scale of each variable given the number of bits.
    const uint64_t fullScale =  ((uint64_t)1 << nb) - 1;

    /**************************************************************************
                                   Variables
    ***************************************************************************/

    /// Variable to store if the design variable limits were already set.
    bool limitsSet_;

    /// Structure to store the minimum and maximum values of each variable.
    sDesignVariable<Scalar> designVars[N];

    /// String.
    std::bitset<N*nb> string;

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

    /// Uniform distribution, integers, 0 to N*nb-1.
    std::uniform_int_distribution<int> rand_string_bit_{
        std::uniform_int_distribution<int>(0, N*nb-1)
            };

    /// Uniform distribution, double, 0 to 1.
    std::uniform_real_distribution<double> rand_{
        std::uniform_real_distribution<double>(0, 1)
            };
};

#include "MGEOcppMain.hpp"

}

#endif // MGEOCPP_HPP
 
