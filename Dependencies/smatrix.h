/*
*
      This file is part of the CNUCTRAN library

      @author   M. R. Omar (rabieomar@usm.my)
      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2023, Universiti Sains Malaysia

      This header file contains  the  definitions of SMATRIX class.  SMATRIX enables fast,
      parallelized and accurate sparse binary exponentiation using the arbitrary precision
      floating-point library, MPFR.

 */

#ifndef SMATRIX_H
#define SMATRIX_H

#include <cnuctran.h>
#include <ppl.h>

using namespace mpfr;
using namespace concurrency;



namespace cnuctran
{

    class smatrix {

    public:

        std::pair<int, int> shape;
        cmap_2d nzel;
        mpfr_prec_t bits = digits2bits(__dps__);
        /*
            Constructor definitions.
        */
        smatrix(void) { return; }
        smatrix(std::pair<int, int> shape) { this->shape = shape; return; }
        smatrix(std::pair<int, int> shape, cmap_2d& A)
        {
            this->shape = shape; this->nzel = A; return;
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

        // Parallel implementation of self sparse matrix-matrix multiplication.
        cmap_2d smul(void)
        {
            cmap_2d result;
            auto& r = result;
            parallel_for_each(begin(nzel), end(nzel), [&](const std::pair<int, cmap_1d>& p)
                {
                    mpreal::set_default_prec(bits);
                    cmap_1d c;
                    for (const auto& [k1, v1] : p.second)
                        for (const auto& [k2, v2] : nzel[k1])
                            c[k2] += v1 * v2;
                    r.insert(make_pair(p.first, c));
                });
            return result;
        }


        // Implementation of exponentiation by squaring.
        void binpow(int k)
        {
            
            for (int i = 1; i <= k; i++)
                nzel = smul();
        }

    };
}

#endif
