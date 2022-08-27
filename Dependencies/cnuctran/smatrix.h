/*
*
      This file is part of the CNUCTRAN library

      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2022 M. R. Omar

      This header file contains  the  definitions of SMATRIX class.  SMATRIX enables fast,
      parallelized and accurate sparse binary exponentiation using the arbitrary precision
      floating-point library, MPFR.

 */

#ifndef SMATRIX_H
#define SMATRIX_H

#include <iostream>
#include <mpreal.h>
#include <unordered_map>
#include <cnuctran/cnuctran.h>
#include <ppl.h>
#include <concurrent_unordered_map.h>

using namespace mpfr;
using namespace concurrency;



namespace cnuctran
{

    class smatrix {

    public:

        std::pair<int, int> shape;
        cmap_2d nzel;

        /*
            Constructor definitions.
        */
        smatrix(void) { return; }
        smatrix(std::pair<int, int> shape) { this->shape = shape; return; }
        smatrix(std::pair<int, int> shape, cmap_2d& A)
        {
            this->shape = shape; this->nzel = A; return;
        }

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
            smatrix result = smatrix(std::pair<int, int>(sx, sy));


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

        // Serial implementation of self sparse multiplication. This routine is depricated.
        smatrix smul_depricated(void)
        {
            smatrix result(shape);
            auto& r = result.nzel;

            for (auto& [a1, a2] : nzel)
            {
                cmap_1d c;
                for (const auto& [k1, v1] : a2)
                    for (const auto& [k2, v2] : nzel[k1])
                        c[k2] += v1 * v2;
                r[a1] = c;
            }

            return result;
            
        }


        // Parallel implementation of self sparse multiplication.
        smatrix smul(void)
        {
            smatrix result(shape);
            auto& r = result.nzel;
            auto bits = digits2bits(__dps__);
            parallel_for_each(nzel.begin(), nzel.end(), [&](std::pair<int, cmap_1d> p)
                {
                    mpreal::set_default_prec(bits);
                    cmap_1d c;
                    
                    for (const auto& [k1, v1] : p.second)
                        for (const auto& [k2, v2] : nzel[k1])
                            c[k2] += v1 * v2;

                    r[p.first] = c;

                }, auto_partitioner());

            return result;

        }


        // Implementation of exponentiation by squaring.
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
