/*
* 
      This file is part of the CNUCTRAN library
 
      @license  MIT
      @author   M. R. Omar
      @link     https://github.com/rabieomar92/cnuctran

      This header file contains the definitions of SMATRIX, which is a component
      of CNUCTRAN. SMATRIX enables fast and accurate sparse binary exponentiation
      using the arbitrary precision floating-point library, MPFR.

 */

#ifndef SMATRIX_H
#define SMATRIX_H

#include <iostream>
#include <mpfr.h>
#include <mpir.h>
#include <mpreal.h>
#include <unordered_map>
#include <cnuctran.h>

using namespace std;
using namespace mpfr;

namespace cnuctran
{
 
    class smatrix {

    public:

        pair<int, int> shape;
        map_2d nzel;

        /*
            Constructor definitions.
        */
        smatrix(void) { return; }
        smatrix(pair<int, int> shape) { this->shape = shape; return; }
        smatrix(pair<int, int> shape, map_2d& A) { this->shape = shape; this->nzel = A; return; }

        smatrix copy()
        {
            smatrix r = smatrix(this->shape);
            r.nzel = this->nzel;
            return r;
        }

        smatrix mul(smatrix& other)
        {
            int const& sx = this->shape.first;
            int const& sy = other.shape.second;
            smatrix result = smatrix(pair<int, int>(sx, sy));


            int row;
            for (row = 0; row < sx; row++)
            {
                auto& c = result.nzel[row];
                for (const auto& [k1, v1] : this->nzel[row])
                    for (const auto& [k2, v2] : other.nzel[k1])
                        c[k2] += v1 * v2;
            }
            return result;
        }

        smatrix smul(void)
        {
            smatrix result(shape);
            int row;
            for (row = 0; row < shape.first; row++)
            {
                auto& c = result.nzel[row];
                for (const auto& [k1, v1] : nzel[row])
                    for (const auto& [k2, v2] : nzel[k1])
                        c[k2] += v1 * v2;

            }
            return result;
        }

        smatrix binpow(int k)
        {
            auto r = copy();
            for (int i = 1; i <= k; i++)
                r = r.smul();
            return r;
        }
    };
}

#endif
