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
        return string_;
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
    void printParetoFrontier(std::ostream& outStream = std::cout) const;
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

    /// Maximum number of evaluations of the objective function.
    int nfobMax_;

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

    /// Structure to store the minimum and maximum values of each variable.
    sDesignVariable<Scalar> designVars[N];

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

    /// Uniform distribution, integers, 0 to nb-1.
    std::uniform_int_distribution<int> rand_string_bit_{
        std::uniform_int_distribution<int>(0, nb-1)
            };

    /// Uniform distribution, double, 0 to 1.
    std::uniform_real_distribution<double> rand_{
        std::uniform_real_distribution<double>(0, 1)
            };
};

#include "MGEOcppMain.hpp"

}

#endif // MGEOCPP_HPP
 
