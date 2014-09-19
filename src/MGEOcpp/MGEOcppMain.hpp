/**
 * @file MGEOcppMain.hpp
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-15
 *
 * @copyright Copyright (c) 2014, Ronan Arraes Jardim Chagas. All rights
 * reserved.
 * @license This project is released under the BSD 3-Clause License (see LICENSE
 * file).
 *
 * @brief Main implementation file of the class MGEO.
 *
 */

/**
 * @brief Constructor of the class MGEO.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] tau Parameter to set the search determinism.
 * @param[in] ngenMax Maximum number of generations.
 * @param[in] runMax Maximum number of independent runs (reinitializations).
 * @param[in] rng_seed Seed for the random number generator.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
MGEO<N, nb, nf, Scalar>::MGEO(double tau, int ngenMax, int runMax, int rng_seed)
    : objectiveFunctions(NULL),
      tau_(tau),
      ngenMax_(ngenMax),
      runMax_(runMax),
      designVarsConfigured_(false),
      designVarsNb_(0),
      rng_(rng_seed)
{
}

/**
 * @brief Constructor of the class MGEO.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] tau Parameter to set the search determinism.
 * @param[in] ngenMax Maximum number of generations.
 * @param[in] runMax Maximum number of independent runs (reinitializations).
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 *
 * @remarks In this case, the random number generator will be initialized using
 * the std::random_device.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
MGEO<N, nb, nf, Scalar>::MGEO(double tau, int ngenMax, int runMax)
    : objectiveFunctions(NULL),
      tau_(tau),
      ngenMax_(ngenMax),
      runMax_(runMax),
      designVarsConfigured_(false),
      designVarsNb_(0)
{
    std::random_device rd;

    // Initialize the random number generator using the std::random_device.
    rng_ = std::mt19937(rd());
}

/**
 * @brief Destructor of the class MGEO
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
MGEO<N, nb, nf, Scalar>::~MGEO()
{
}

/******************************************************************************
                              Getters and Setters
 ******************************************************************************/

/**
 * @brief Configure the design variables.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] bits List containing the number of bits for each design variable.
 * @param[in] min List containing the minimum of each design variable.
 * @param[in] max List containing the maximum of each design variable.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::confDesignVars(
    std::initializer_list<unsigned int> bits,
    std::initializer_list<Scalar> min, 
    std::initializer_list<Scalar> max,
    std::initializer_list<std::string> varNames) 
{
    // Pointer to the beginning of the list that contains the number of bits.
    auto p_bits = std::begin(bits);

    // Pointer to the beginning of the list that contains the minimums.
    auto p_min = std::begin(min);

    // Pointer to the beginning of the list that contains the maximums.
    auto p_max = std::begin(max);

    // Pointer to the beginning of the list that contains the variable names.
    auto p_name = std::begin(varNames);

    // Pointer to count the total number of bits used.
    unsigned int numBits = 0;

    // Variable to check if the design variables name were defined.
    bool namesDefined = (varNames.size() != 0) ? true : false;

    // Check the list size.
    if ( (min.size() != N) || (max.size() != N) || 
         ( namesDefined && (varNames.size() != N) ) )
    {
        std::cerr << "MGEO::confDesignVars(): The number of elements on the lists in function confDesignVars() must be " 
                  << N << "." << std::endl;
        return false;
    }
    
    // Add the limits to the designVars variable.
    for(int i = 0; i < N; i++)
    {
        designVars[i].bits      = (unsigned int)*p_bits;
        designVars[i].min       = (Scalar)*p_min;
        designVars[i].max       = (Scalar)*p_max;
        designVars[i].fullScale = ((uint64_t)1 << designVars[i].bits) - 1;
        designVars[i].index     = numBits;

        if(!namesDefined)
            designVars[i].name = std::string("VARS "+std::to_string(i));
        else
            designVars[i].name = *p_name;

        if( *p_min >= *p_max )
        {
            std::cerr << "MGEO::confDesignVars(): The elements on the min list must be smaller than the elements on list max."
                      << std::endl;
            
            designVarsConfigured_ = false;
            return false;
        }

        if( *p_bits == 0 )
        {
            std::cerr << "MGEO::confDesignVars(): The number of bits for a design variable must not be zero."
                      << std::endl;

            designVarsConfigured_ = false;
            return false;
        }

        numBits += *p_bits;
        p_bits++;
        p_min++;
        p_max++;
        p_name++;
    }

    // Check if the number of bits were configured properly.
    if( numBits > nb )
    {
        std::cerr << "MGEO::confDesignVars(): The number of bits configured for the design variables is bigger than the total number of bits allocated."
                  << std::endl;

        designVarsConfigured_ = false;
        return false;
    }

    designVarsNb_ = numBits;
    designVarsConfigured_ = true;
    return true;
}

/**
 * @brief Set the limits for all variables.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] bits The number of bits for each variable.
 * @param[in] min The minimum value for all design variables.
 * @param[in] max The maximum value for all design variables.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::confDesignVars(unsigned int bits, 
                                             Scalar min, 
                                             Scalar max)
{
    // Check the number of bits.
    if ( bits*N > nb )
    {
        std::cerr << "MGEO::confDesignVars(): The number of bits configured for the design variables is bigger than the total number of bits allocated."
                  << std::endl;
        designVarsConfigured_ = false;
        return false;
    }
    
    // Check the limits.
    if ( min >= max )
    {
        std::cerr << "MGEO::confDesignVars(): The minimum value must be smaller than the maximum value."
                  << std::endl;
        designVarsConfigured_ = false;
        return false;
    }

    // Add the limits to the designVars variable.
    for(int i = 0; i < N; i++)
    {
        designVars[i].bits      = bits;
        designVars[i].min       = min;
        designVars[i].max       = max;
        designVars[i].fullScale = ((uint64_t)1 << designVars[i].bits) - 1;
        designVars[i].index     = i*bits;
        designVars[i].name      = std::string("VARS "+std::to_string(i));
    }

    designVarsNb_ = bits*N;
    designVarsConfigured_ = true;
    return true;
}

/******************************************************************************
                                   Functions
 ******************************************************************************/

/**
 * @brief Transform the string to real values and call the objective functions.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] string String.
 * @param[out] vars Real values of the design variables that are obtained using
 * the current string and the limits in designVars.
 * @param[out] f Values of the objective functions evaluated at <tt>vars</tt>.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 *
 * @retval TRUE The objective functions was successfully called.
 * @retval FALSE There was a problem calling the objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::callObjectiveFunctions(std::bitset<nb> string, 
                                                     Scalar *vars, 
                                                     Scalar *f)
{
    // Check if the output parameters are NULL.
    if ( (vars == NULL) || (f == NULL) )
    {
        std::cerr << "MGEO::callObjectiveFunctions(): The pointers of the output variables must not be NULL."
                  << std::endl;

        return false;
    }

    // Check if the objective functions were set.
    if (objectiveFunctions == NULL)
    {
        std::cerr << "MGEO::callObjectiveFunctions(): The pointer objectiveFunctions was not set."
                  << std::endl;
        
        return false;
    }

    // Transform the string of bits into scalars.
    stringToScalar(string, vars);

    // Call the objective functions.
    return objectiveFunctions(vars, f);
}

/**
 * @brief Check the dominance of a Pareto point.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] p Pareto point that will be compared with the list in
 * <tt>paretoFrontier</tt>.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 *
 * @remarks This function edits the list <tt>paretoFrontier</tt>. If <tt>p</tt>
 * is not dominated by any point in the list, then <tt>p</tt> is added to it and
 * all points that become dominated by <tt>p</tt> are deleted.
 *
 * @retval TRUE The candidate point <tt>p</tt> was added to the list.
 * @retval FALSE The candidate point <tt>p</tt> was not added to the list.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::checkDominance(sParetoPoint<N, nf, Scalar> p)
{
    bool addPoint = false;

    // Loop for the entire list of Pareto points.
    auto it = paretoFrontier.begin();

    while (it != paretoFrontier.end())
    {
        /* Variable to check if the point in the Pareto frontier list is
           dominated by the candidate. */
        bool pointDominated = true;
        bool candidateDominated = true;

        for(int i = 0; i < nf; i++)
        {
            if( p.f[i] < it->f[i] )
            {
                /* The candidate point is not dominated by the point in the
                   list. */
                candidateDominated = false;
            }
            else if( p.f[i] > it->f[i] )
            {
                /* The point in the list is not dominated by the candidate
                   point. */
                pointDominated = false;
            }
        }

        /* If the point is dominated by the candidate and dominates the
           candidate at the same time, then the candidate and the point is the
           same. Thus, mark to do not add the candidate and exit the loop. */
        if(pointDominated && candidateDominated)
        {
            addPoint = false;
            break;
        }

        /* If the point in the Pareto frontier is dominated, then exclude it
           from the list. */
        if (pointDominated)
            it = paretoFrontier.erase(it);
        else
            ++it;

        /* If the candidate is dominated by any point in the list, stop the
           search and do not add it. */
        if (candidateDominated)
        {
            addPoint = false;
            break;
        }
        // Otherwise, add the point to the list.
        else
            addPoint = true;
    }

    // Check if the new point must be added.
    if (addPoint)
        paretoFrontier.push_back(p);

    return addPoint;
}

/**
 * @brief Initialize the string with uniformly distributed bits.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
void MGEO<N, nb, nf, Scalar>::initializeString()
{        
    for(int i = 0; i<designVarsNb_; i++)
        string_[i] = rand_bit_(rng_);
}

/**
 * @brief Print the Pareto frontier.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-19
 *
 * @param[in] outStream Stream to which the list will be printed (default =
 * std::cout).
 * @param[in] fieldWidth Field width for each printed value.
 * @param[in] precision Precision of printed values.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<unsigned int N,
         unsigned int nb,
         unsigned int nf,
         typename Scalar>
void MGEO<N, nb, nf, Scalar>::printParetoFrontier(std::ostream& outStream,
                                                  unsigned int fieldWidth,
                                                  unsigned int precision) const
{
    // Header.
    for(int i = 0; i < N; i++)
        outStream << std::setw(fieldWidth) << designVars[i].name;

    for(int i = 0; i < nf; i++)
      outStream << std::setw(fieldWidth-1) << "FOBJ " << i;
    outStream << std::endl;

    outStream.precision(precision);

    // Pareto Frontier.
    for(auto it = paretoFrontier.begin(); it != paretoFrontier.end(); it++)
    {
        for(int i = 0; i < N; i++)
            outStream << std::setw(fieldWidth) << (*it).vars[i] ;
        for(int i = 0; i < nf; i++)
            outStream << std::setw(fieldWidth) << (*it).f[i];
        outStream << std::endl;
    }
}

/**
 * @brief Run the MGEO.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N,
         unsigned int nb,
         unsigned int nf,
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::run()
{
    // Rank.
    std::pair<int, Scalar> fRank[nb];

    // Variables to create the Pareto point.
    Scalar vars[N];
    Scalar f[nf];
    sParetoPoint<N, nf, Scalar> paretoPoint;

    // Number of independent runs.
    int run = 0;

    // Number of generations created per run.
    int ngenPerRun = 0;

    // Maximum number generations created per run.
    int ngenMaxPerRun = std::floor(ngenMax_/runMax_);

    // If true, the algorithm will be reinitialized.
    bool reinitialize = true;

    // Check if the project variable limits were set.
    if (!designVarsConfigured_)
    {
        std::cerr << "MGEO::run(): The design variables were not properly configured (see MGEO::confDesignVars())."
                  << std::endl;

        return false;
    }

    // Clear the Pareto frontier.
    paretoFrontier.clear();

#ifdef MGEOCPP_DEBUG
    int count = 0;
    std::cout << "MGEO: Run " << run+1 << "..." << std::endl;
#endif // MGEOCPP_DEBUG

    // Loop.
    while(run < runMax_)
    {
#ifdef MGEOCPP_DEBUG
        if( ngenPerRun - count > 0.1*ngenMaxPerRun)
        {
            std::cout << static_cast<int>(ngenPerRun*100/ngenMaxPerRun) << "%" << "...";
            std::cout.flush();
            count = ngenPerRun;
        }
#endif // MGEOCPP_DEBUG
        
        /**************************************************************************
                                 Initialization
        **************************************************************************/
        if( reinitialize == true )
        {
            // Sample a new string.
            initializeString();

            // Call the objective functions for the first time.
            if (!callObjectiveFunctions(string_, vars, f))
                return false;

            ngenPerRun++;

            // Add the results to the list of Pareto points in the first run.
            if (paretoFrontier.size() == 0)
            {
                std::copy(&vars[0], &vars[N], paretoPoint.vars);
                std::copy(&f[0], &f[nf], paretoPoint.f);
                paretoFrontier.push_back(paretoPoint);
            }
        }

        /**************************************************************************
                            Adaptability and Ranking
        **************************************************************************/

        /* Choose which objective function will be used to compute the adaptability
           and to assemble the rank. */
        int chosenFunc = rand_f_(rng_);
        
        // Copy of string to be used in the parallel loop.
        std::bitset<nb> string_p = string_;

        // Check if there was a problem with the objective functions.
        bool objectiveFunctionsProblem = false;
        
        // List of all points created after flipping the bits.
        sParetoPoint<N,nf,Scalar> candidatePoints[nb];

#pragma omp parallel for schedule(dynamic) firstprivate(string_p) private(vars, f, paretoPoint)
        for(int i = 0; i < designVarsNb_; i++)
        {
            // Toggle the j-th bit of the string.
            string_p.flip(i);
                
            // Compute the objective functions. 
            if(!callObjectiveFunctions(string_p, vars, f))
                objectiveFunctionsProblem = true;
            
            // Create the candidate Pareto point.
            std::copy(&vars[0], &vars[N], candidatePoints[i].vars);
            std::copy(&f[0], &f[nf], candidatePoints[i].f);
            
            // Add the result to the rank.
            fRank[i].first  = i;
            fRank[i].second = f[chosenFunc];

            /* Notice that the bit does not need to unflip because string_p is a
               copy of string on each thread. */
        }

        // Check if there was a problem in the objective functions.
        if(objectiveFunctionsProblem)
            return false;

        // Add the points to the Pareto frontier.
        for(int i = 0; i < designVarsNb_; i++)
            checkDominance(candidatePoints[i]);

        // Ranking.
        std::sort(fRank, fRank+designVarsNb_,
            [](std::pair<int, Scalar> const &a, std::pair<int, Scalar> const &b)
            {
                return a.second < b.second;
            });

        // Choose a bit to be changed for the new "population".
        bool bitAccepted = false;

        while (bitAccepted == false)
        {
            int b_sample = rand_string_bit_(rng_);

            /* Accept the change with probability r^(-tal), where r is the rank
               of the bit. */
            double Pk = std::pow(b_sample+1, -tau_);

            if ( rand_(rng_) <= Pk  )
            {
                // If the toggle is accepted, then exit the loop.
                string_.flip(fRank[b_sample].first);
                bitAccepted = true;
            }
        }

        // A new generation has been created.
        ngenPerRun++;

        /* If the number of generations created is more than ngenMaxPerRun, then
           reinitialize the algorithm. */
        if (ngenPerRun > ngenMaxPerRun)
        {
            ngenPerRun = 0;
            reinitialize = true;

#ifdef MGEOCPP_DEBUG
            std::cout << "100%" << std::endl << std::endl;
            std::cout << "Pareto frontier at run " << run+1 << ":" 
                      << std::endl
                      << std::endl;
            std::cout << "    Number of points: " 
                      << paretoFrontier.size() 
                      << std::endl;
            std::cout << "    Memory (MB):      " 
                      << paretoFrontier.size()*sizeof(sParetoPoint<N, nf, Scalar>)/1024.0/1024.0
                      << std::endl << std::endl << std::endl;

            count = 0;

            if( run < runMax_)
                std::cout << "MGEO: Run " << run+1 << " of " << runMax_ << "..." << std::endl;
#endif // MGEOCPP_DEBUG

            run += 1;
        }
    }

    return true;
}

/**
 * @brief Sort the points in the Pareto Frontier.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-25
 *
 * @param[in] fobj Function to sort.
 *
 * @retval TRUE The list was sorted.
 * @retval FALSE The list was not sorted.
 */
template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::sortParetoFrontier(int fobj)
{
    // Check the parameter.
    if( (fobj < 1) || (fobj > nf) )
    {
        std::cerr << "sortParetofrontier(): fobj must be between 1 and the number of objective functions."
                  << std::endl;
        return false;
    }

    paretoFrontier.sort(
        [&fobj](sParetoPoint<N, nf, Scalar> const &a, 
                sParetoPoint<N, nf, Scalar> const &b)
        {
            return a.f[fobj-1] < b.f[fobj-1];
        });

    return true;
}

/**
 * @brief Convert the string to real values.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] string String.
 * @param[out] vars The converted design variables.
 *
 * @tparam N Number of design variables.
 * @tparam nb Total number of bits allocated for the string
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
void MGEO<N, nb, nf, Scalar>::stringToScalar(std::bitset<nb> string, 
                                             Scalar* vars)
{
    // Loop for each project variable.
    for(int i = 0; i < N; i++)
    {
        // Convert the bits to integer.
        uint64_t varInt = 0;

        // Number of bits allocated for the design variable.
        unsigned int numBits = designVars[i].bits;

        // Full scale of the design variable.
        uint64_t fullScale = designVars[i].fullScale;

        // Index of the design variable in the string.
        uint64_t index = designVars[i].index;

        for(int j = 0; j < numBits; j++)
            varInt += string[index + j]*std::pow(2, j);

        // Convert to type Scalar.
        vars[i] = designVars[i].min + 
            (designVars[i].max-designVars[i].min)*varInt/fullScale;
    }
}