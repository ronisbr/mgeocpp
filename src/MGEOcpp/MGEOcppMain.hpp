/**
 * @file MGEOcppMain.hpp
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-15
 * @brief Main implementation file of the class MGEO.
 *
 */

/**
 * @brief Constructor of the class MGEO.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] tau Parameter to set the search determinism.
 * @param[in] nfobMax Maximum number of evaluations of the objective function.
 * @param[in] runMax Maximum number of independent runs (reinitializations).
 * @param[in] rng_seed Seed for the random number generator.
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
MGEO<N, nb, nf, Scalar>::MGEO(double tau, int nfobMax, int runMax, int rng_seed)
    : objectiveFunctions(NULL),
      tau_(tau),
      nfobMax_(nfobMax),
      runMax_(runMax),
      limitsSet_(false),
      rng_(rng_seed)
{
}

/**
 * @brief Constructor of the class MGEO.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] tau Parameter to set the search determinism.
 * @param[in] nfobMax Maximum number of evaluations of the objective function.
 * @param[in] runMax Maximum number of independent runs (reinitializations).
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
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
MGEO<N, nb, nf, Scalar>::MGEO(double tau, int nfobMax, int runMax)
    : objectiveFunctions(NULL),
      tau_(tau),
      nfobMax_(nfobMax),
      runMax_(runMax),
      limitsSet_(false)
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
 * @tparam nb Number of bits per design variable.
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
 * @brief Set the limits for each design variable.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[in] min List containing the minimum of each design variable.
 * @param[in] max List containing the maximum of each design variable.
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
bool MGEO<N, nb, nf, Scalar>::setDesignVarsLimits(
    std::initializer_list<Scalar> min, 
    std::initializer_list<Scalar> max) 
{
    // Pointer to the beginning of the list that contains the minimums.
    auto p_min = std::begin(min);

    // Pointer to the beginning of the list that contains the maximums.
    auto p_max = std::begin(max);
    
    // Check the list size.
    if ( (min.size() != N) || (max.size() != N) )
    {
        std::cerr << "MGEO::setDesignVars(): The number of elements on the lists in function setDesignVars() must be " 
                  << N << "." << std::endl;
        return false;
    }
    
    // Add the limits to the designVars variable.
    for(int i = 0; i < N; i++)
    {
        designVars[i].min = (Scalar)*p_min;
        designVars[i].max = (Scalar)*p_max;

        if( *p_min >= *p_max )
        {
            std::cerr << "MGEO::setDesignVars(): The elements on the min list must be smaller than the elements on list max."
                      << std::endl;
            
            return false;
        }
        
        p_min++;
        p_max++;
    }

    limitsSet_ = true;
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
 * @param[out] vars Real values of the design variables that are obtained using
 * the current string and the limits in designVars.
 * @param[out] f Values of the objective functions evaluated at <tt>vars</tt>.
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
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
bool MGEO<N, nb, nf, Scalar>::callObjectiveFunctions(Scalar *vars, Scalar *f)
{
    // Check if the output parameters are NULL.
    if ( (vars == NULL) || (f == NULL) )
    {
        std::cerr << "MGEO::callObjectiveFunctions(): The pointers of the output variables must not be NULL."
                  << std::endl;

        return false;
    }

    // Transform the string of bits into scalars.
    stringToScalar(vars);

    // Check if the objective functions were set.
    if (objectiveFunctions == NULL)
    {
        std::cerr << "MGEO::callObjectiveFunctions(): The pointer objectiveFunctions was not set."
                  << std::endl;
        
        return false;
    }

    // Call the objective functions.
    if (!objectiveFunctions(vars, f))
    {
        std::cerr << "MGEO::callObjectiveFunctions(): The objective functions returned an error."
                  << std::endl;

        return false;
    }

    return true;
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
 * @tparam nb Number of bits per design variable.
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
                /* The candidate point is not dominated by the point in the
                   list. */
                candidateDominated = false;
            else
                /* The point in the list is not dominated by the candidate
                   point. */
                pointDominated = false;
        }

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

        /* If the point in the Pareto frontier is dominated, then exclude it
           from the list. */
        if( pointDominated )
            it = paretoFrontier.erase(it);
        else
            ++it;
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
 * @tparam nb Number of bits per design variable.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */

template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
void MGEO<N, nb, nf, Scalar>::initializeString()
{        
    for(int i = 0; i<N*nb; i++)
        string[i] = rand_bit_(rng_);
}

/**
 * @brief Print the Pareto frontier.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-19
 *
 * @param[in] outStream Stream to which the list will be printed (default =
 * std::cout).
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<unsigned int N,
         unsigned int nb,
         unsigned int nf,
         typename Scalar>
void MGEO<N, nb, nf, Scalar>::printParetoFrontier(std::ostream& outStream) const
{
    const int fieldWidth = 20;

    // Header.
    for(int i = 0; i < N; i++)
        outStream << std::setw(fieldWidth-1) << "VAR " << i << ' ' ;
    for(int i = 0; i < nf; i++)
        outStream << std::setw(fieldWidth-1) << "FOBJ " << i << ' ' ;
    outStream << std::endl;

    // Pareto Frontier.
    for(auto it = paretoFrontier.begin(); it != paretoFrontier.end(); it++)
    {
        for(int i = 0; i < N; i++)
            outStream << std::setw(fieldWidth) << (*it).vars[i] << ' ' ;
        for(int i = 0; i < nf; i++)
            outStream << std::setw(fieldWidth) << (*it).f[i] << ' ' ;
        outStream << std::endl;
    }
}

/**
 * @brief Run the MGEO.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
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
    std::pair<int, Scalar> fRank[N*nb];

    // Variables to create the Pareto point.
    Scalar vars[N];
    Scalar f[nf];
    sParetoPoint<N, nf, Scalar> paretoPoint;

    // Number of independent runs.
    int run = 0;

    // Number of evaluations of the objective functions per run.
    int nfobPerRun = 0;

    // Maximum number of evaluations of the objective functions per run.
    int nfobRunMax = int(nfobMax_/runMax_);

    // If true, the algorithm will be reinitialized.
    bool reinitialize = true;

    // Check if the project variable limits were set.
    if (!limitsSet_)
    {
        std::cerr << "MGEO::run(): The project variable limits must be set before calling run()."
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
        if( nfobPerRun - count > 0.1*nfobRunMax)
        {
            std::cout << static_cast<int>(nfobPerRun*100/nfobRunMax) << "%" << "...";
            std::cout.flush();
            count = nfobPerRun;
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
            if (!callObjectiveFunctions(vars, f))
                return false;

            nfobPerRun += nf;

            // Add the results to the list of Pareto points in the first run.
            if (run == 0)
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
        
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < nb; j++)
            {
                // Toggle the j-th bit of the string.
                string.flip(i*nb+j);
                
                // Compute the objective functions. 
                if (!callObjectiveFunctions(vars, f))
                    return false;
                nfobPerRun += nf;
                
                // Create the candidate Pareto point.
                std::copy(&vars[0], &vars[N], paretoPoint.vars);
                std::copy(&f[0], &f[nf], paretoPoint.f);
                
                // Check the dominance of the Pareto point.
                checkDominance(paretoPoint);
                
                // Add the result to the rank.
                fRank[i*nb+j].first  = i*nb+j;
                fRank[i*nb+j].second = f[chosenFunc];
                
                // Reset the bit.
                string.flip(i*nb+j);
            }
        }

        // Ranking.
        std::sort(std::begin(fRank), std::end(fRank),
            [](std::pair<int, Scalar> const &a, std::pair<int, Scalar> const &b)
            {
                return a.second <= b.second;
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
                string.flip(fRank[b_sample].first);
                bitAccepted = true;
            }
        }

        /* If the objective function is evaluated more than nfobRunMax, then
           reinitialize the algorithm. */
        if (nfobPerRun > nfobRunMax)
        {
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

            std::cout << "MGEO: Run " << run+1 << "..." << std::endl;
            count = 0;
#endif // MGEOCPP_DEBUG

            nfobPerRun = 0;
            reinitialize = true;
            run += 1;
        }
    }

    return true;
}

/**
 * @brief Convert the string to real values.
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-18
 *
 * @param[out] vars The converted design variables.
 *
 * @tparam N Number of design variables.
 * @tparam nb Number of bits per design variable.
 * @tparam nf Number of objective functions.
 * @tparam Scalar Scalar type of design variables and objective functions.
 */
template<unsigned int N, 
         unsigned int nb, 
         unsigned int nf, 
         typename Scalar>
void MGEO<N, nb, nf, Scalar>::stringToScalar(Scalar* vars)
{
    // Loop for each project variable.
    for(int i = 0; i < N; i++)
    {
        // Convert the bits to integer.
        uint64_t varInt = 0;

        for(int j = 0; j < nb; j++)
            varInt += string[i*nb + j]*std::pow(2, j);

        // Convert to type Scalar.
        vars[i] = designVars[i].min + 
            (designVars[i].max-designVars[i].min)*varInt/fullScale;
    }
}