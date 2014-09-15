/**
 * @file mgeocpp_testcases.cpp
 * @author Ronan Arraes Jardim Chagas
 * @date 2014-08-15
 *
 * @copyright Copyright (c) 2014, Ronan Arraes Jardim Chagas. All rights
 * reserved.
 * @license This project is released under the BSD 3-Clause License (see LICENSE
 * file).
 *
 * @brief Test cases for the MGEOcpp.
 *
 * This file defines the test cases for the MGEOcpp.
 * The scenarios were obtained from:
 *
 *     Galski, R. L (2006). Desenvolvimento de Versões Aprimoradas, Híbridas,
 *     Paralela e Multiobjetivo do Método da Otimização Extrema Generalizada e
 *     Sua Aplicação no Projeto de Sistemas Espaciais. Ph.D. Thesis. Instituto
 *     Nacional de Pesquisas Espaciais: Brazil.
 *
 */

#include "MGEOcpp/MGEOcpp.hpp"

using namespace MGEOcpp;

bool objectiveFunctions1(double *vars, double *f)
{
    double x1 = vars[0];

    f[0] = exp(-x1) + 1.4*exp(-x1*x1);
    f[1] = exp(+x1) + 1.4*exp(-x1*x1);

    return true;
}

bool objectiveFunctions2(double *vars, double *f)
{
    double x1 = vars[0];

    f[0] = x1*x1;
    f[1] = (x1-2)*(x1-2);

    return true;
}

bool objectiveFunctions3(double *vars, double *f)
{
    double x1 = vars[0];

    if (x1 <= 1.0)
        f[0] = -x1;
    else if ( (1.0 < x1) && (x1 <= 3.0) )
        f[0] = -2.0+x1;
    else if ( (3.0 < x1) && (x1 <= 4.0) )
        f[0] = 4.0-x1;
    else
        f[0] = -4.0+x1;

    f[1] = (x1-5)*(x1-5);

    return true;
}

bool objectiveFunctions4(double *vars, double *f)
{
    double x = vars[0];
    double y = vars[1];

    f[0] = x*x+y*y;
    f[1] = (x-1.5)*(x-1.5)+(y-3.5)*(y-3.5);

    return true;
}

int main(int argc, char *argv[])
{
    const int testcase = 4;

    // For cases 1, 2, and 3.
    MGEO<1, 16, 2> mgeo1v(0.5, 15000, 50);

    // For case 4.
    MGEO<2, 32, 2> mgeo2v(0.5, 15000, 50);

    switch(testcase)
    {
    case 1:
        mgeo1v.confDesignVars({16},{-10},{+10});
        mgeo1v.objectiveFunctions = &objectiveFunctions1;
        mgeo1v.run();
        break;
    case 2:
        mgeo1v.confDesignVars({16},{-10},{+10});
        mgeo1v.objectiveFunctions = &objectiveFunctions2;
        mgeo1v.run();
        break;
    case 3:
        mgeo1v.confDesignVars({16},{-10},{+10});
        mgeo1v.objectiveFunctions = &objectiveFunctions3;
        mgeo1v.run();
        break;
    case 4:
        mgeo2v.confDesignVars({8,8},{-10,-10},{+10,+10},{"var x","var y"});
        mgeo2v.objectiveFunctions = &objectiveFunctions4;
        mgeo2v.run();
        break;
    default:
        std::cerr << "Invalid test case!" << std::endl;
        return 1;
    }

    if(testcase != 4)
    {
        mgeo1v.sortParetoFrontier();
        mgeo1v.printParetoFrontier();
    }
    else
    {
        mgeo2v.sortParetoFrontier();
        mgeo2v.printParetoFrontier();
    }

    return 0;
}
