/*

      This file is part of the CNUCTRAN library

      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2022 M. R. Omar

      This header file contains the definitions of all reusable computational 
      routines of the proposed probabilistic method, which is implemented in
      CNUCTRAN.

 */

#ifndef SOLVER_H
#define SOLVER_H


#include <iostream>
#include <mpreal.h>
#include <cnuctran/smatrix.h>
#include <cnuctran/cnuctran.h>
#include <chrono>
#include <map>
#include <unordered_map>


using namespace mpfr;

namespace cnuctran
{
    
    /*
            This class was initially developed to facilitate fast, high-precision sparse
            matrix multiplications and powers. WARNING! This class does not covers all
            matrix operations, it only cover the basic operations used by CNUCTRAN, i.e.
            Multiplication and Powers.

            Of course, there is still no known library that provides high-precision sparse
            matrix operations. Therefore, I must endure writing a new specialized class
            handling sparse matrix power to preserve the accuracy.

    */

    class solver
    {
    public:

        std::vector<std::string> species_names;
        int __I__;
        std::vector<std::vector<mpreal>> lambdas;
        std::vector<std::vector<std::vector<int>>> G;
        std::vector<std::vector<mpreal>> fission_yields;

        solver(std::vector<std::string> species_names)
        {
            this->species_names = species_names;
            this->__I__ = this->species_names.size();
            for (int i = 0; i < this->__I__; i++)
            {
                this->lambdas.push_back(std::vector<mpreal>());
                std::vector<int> tmp1 = { __nop__ }; std::vector<std::vector<int>> tmp2; tmp2.push_back(tmp1);
                this->G.push_back(tmp2);
                this->fission_yields.push_back(std::vector<mpreal>());
            }
            return;
        }

        void add_removal(int species_index,
            mpreal rate,
            std::vector<int> products,
            std::vector<mpreal> fission_yields = std::vector<mpreal>({}))
        {


//..........Skips adding a removal if the removal rate is outside of the range
//          specified by the input file.
            if (rate < __mnr__ || rate > __mxr__)
                return;

            this->lambdas[species_index].push_back(rate);
            this->G[species_index].push_back(products);

            if (!fission_yields.empty() && products.size() > 1)
            {
                if (fission_yields.size() >= products.size())
                {
                    std::vector<mpreal> tmp = std::vector<mpreal>();
                    for (mpreal y : fission_yields)
                        tmp.push_back(y);
                    this->fission_yields[species_index] = tmp;
                }
                else
                {
                    std::cout << "FATAL-ERROR\t<cnuctran::solver::add_removal(...)> Insufficient fission yields given for species " <<
                        this->species_names[species_index] << " products.";
                    exit(1);
                }
            }
            else if ((fission_yields.empty() && products.size() == 1) || (fission_yields.size() == 1 && products.size() == 1))
            {
                for (int product : products)
                {
                    // This part is for preparing the transmutation matrix, which is not relevant 
                    // for this C++ version, cnuctran. For CRAM calculation, please use the PyNUCTRAN.
                }
            }
            else
            {
                std::cout << "FATAL-ERROR\t<cnuctran::solver::add_removal(...)> Invalid removal definition for isotope " <<
                    this->species_names[species_index] << std::endl;
                std::cout << "INFO\tNon-fission events MUST only have ONE daughter product." << std::endl;
                std::cout << "INFO\tWhereas fission events MUST have >1 products to track." << std::endl;
                exit(1);

            }
        }

        smatrix prepare_transfer_matrix(mpreal dt)
        {
            std::unordered_map<int, std::unordered_map<int, mpreal, modified_hash>, modified_hash> A;
            std::unordered_map<int, std::unordered_map<int, mpreal, modified_hash>, modified_hash> P;
            int i;

            std::unordered_map<int, mpreal, modified_hash> e;
            for (i = 0; i < this->__I__; i++)
            {
//..............Clears the exponentials container, e, and retrieves the total number of
//              events associated to nuclide-i.
                e.clear();
                const int n_events = this->G[i].size();


                //..............Precalculate the exponentials.
                mpreal norm = __zer__;

                for (int l = 1; l < n_events; l++)
                    e.emplace(l - 1, exp(-this->lambdas[i][l - 1] * dt));

//..............Constructs the pi-distribution.
                for (int j = 0; j < n_events; j++)
                {

                    auto& p = P[i][j];
                    p = __one__;
                    for (int l = 1; l < n_events; l++)
                    {
                        p *= l == j ? __one__ - e[l - 1] : e[l - 1];
                    }
                    norm += p;
                }

                if (norm == __zer__)
                    continue;

//..............Constructs the transfer matrix.
                auto const& gI = G[i];
                for (int j = 0; j < n_events; j++)
                {

                    auto const& a = (P[i][j] / norm);
                    auto const gJ = gI[j];
                    int n_daughters = gJ.size();
                    for (int l = 0; l < n_daughters; l++)
                    {
                        auto const& k = gJ[l];
                        if (k != __nop__)
                        {
                            n_daughters > 1? A[k][i] += a * fission_yields[i][l]:
                                 A[k][i] += a;
                        }
                    }

                    if (j == 0) A[i][i] += a;

                }
            }

            return smatrix({this->__I__, this->__I__ }, A);
        }

        std::map<std::string, mpreal> solve(std::map<std::string, mpreal> w0,
            mpreal n,
            mpreal t)
        {
            std::unordered_map<int, std::unordered_map<int, mpreal, modified_hash>, modified_hash> w0_matrix;
            for (int i = 0; i < this->__I__; i++)
                if (w0.count(this->species_names[i]) == 1)
                    w0_matrix[i][0] = w0[this->species_names[i]];
            smatrix converted_w0 = smatrix(std::pair<int, int>(this->__I__, 1), w0_matrix);

//..........Auto suggest the no. of substeps.
            int k = int(floor(log(t / pow(mpreal("10"), -n)) / log(__two__)));
            mpreal suggested_substeps = pow(__two__, k);
            mpz_t substeps; mpz_set_str(substeps, suggested_substeps.toString().c_str(), 10);

//..........Compute the transfer matrix power.
            if (__vbs__) std::cout << "INFO\ttime-step = " << t << " sec(s)."  << std::endl;
            if (__vbs__) std::cout << "INFO\tmax-rate = " << __mxr__ << " per sec.\tmin-rate = " << __mnr__ << " per sec." << std::endl;
            auto t1 = std::chrono::high_resolution_clock::now();
            smatrix T = this->prepare_transfer_matrix(t / suggested_substeps);

            auto t2 = std::chrono::high_resolution_clock::now();
            smatrix U = T.binpow(k);
            smatrix w = U.mul(converted_w0);
            auto t3 = std::chrono::high_resolution_clock::now();
            if (__vbs__) std::cout << "INFO\tDone computing concentrations.";
            if (__vbs__) std::cout << " (" <<
                std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "ms. for " << k << " mults.)" << std::endl;

            std::map<std::string, mpreal> out;
            for (int i = 0; i < this->__I__; i++)
                out[this->species_names[i]] = w.nzel[i][0];

            return out;
        }

    
       
    };
}

#endif
