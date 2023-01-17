/*

      This file is part of the CNUCTRAN library

      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2022 M. R. Omar

      This header file contains the definitions of all reusable constants and
      and variables in CNUCTRAN.

 */

#ifndef CNUCTRAN_H
#define CNUCTRAN_H

#include <unordered_map>
#include <concurrent_unordered_map.h>

using namespace std;
using namespace mpfr;
using namespace concurrency;

namespace cnuctran
{
    /*

        REUSABLE HIGH PRECISION CONSTANTS.
        __two__ = 2
        __one__ = 1
        __neg__ = -1
        __zer__ = 0
        __eps__ is the epsilon value. Any value falls below __eps__ is assumed to be zero.

        REUSABLE INTEGER CONSTANTS.
        __dps__ is the mpreal arithmetic precision.
        __dop__ is the decimal places of any printed mpreal numbers.
        __nop__ is an integer specifying no product.
        __vbs__ is the vervosity level; 0 (none), 1 (minimal), 2 (comprehensive)

        REUSABLE DOUBLE CONSTANTS.
        __mnr__ is the minimum removal rate allowed in the calculation.
        __mxr__ is the maximum removal rate allowed in the calculation.

    */

    mpreal __eps__ = mpreal("1e-200", digits2bits(50));
    double __mnr__ = 1e-200;
    double __mxr__ = 1e+200;
    int    __dps__ = 45;
    const int    __dop__ = 16;
    const int    __npr__ = 1;
    const int    __nop__ = -1;
    int          __vbs__ = 0;

    /*
        Type definition for sparse matrix non-zero elements container.
        The container is a nested unordered map.
    */
    typedef unordered_map<int, unordered_map<int, mpreal>> map_2d;
    typedef unordered_map<int, mpreal> map_1d;
    typedef concurrent_unordered_map<int, concurrent_unordered_map<int, mpreal>> cmap_2d;
    typedef concurrent_unordered_map<int, mpreal> cmap_1d;

    /*
        Enums for exceptions handling.
    */
    enum errex
    {
        MISSING_SUBSTEP_SIZE = 1,
        MISSING_STEP_SIZE = 2,
        MISSING_SPECIES_NAMES = 3,
        MISSING_W0 = 4,
        MISSING_RXN_RATE = 5,
        NUCLIDES_DATA_LOAD_FAILED = 6,
        XML_READING_ERROR = 7,
        UNEXPECTED_ERROR = 8,
        MISSING_W0_SOURCE = 9

    };

    const string WHITESPACE = " \n\r\t\f\v";

    static string ltrim(const string& s)
    {
        size_t start = s.find_first_not_of(WHITESPACE);
        return (start == string::npos) ? "" : s.substr(start);
    }

    static string rtrim(const string& s)
    {
        size_t end = s.find_last_not_of(WHITESPACE);
        return (end == string::npos) ? "" : s.substr(0, end+1);
    }

    static string trim(const string& s)
    {
        return rtrim(ltrim(s));
    }
}


#endif
