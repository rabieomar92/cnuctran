/*

      This file is part of the CNUCTRAN library

      @license  MIT
      @author   M. R. Omar
      @link     https://github.com/rabieomar92/cnuctran

      This header file contains the definitions of all reusable computational 
      routines of the proposed probabilistic method, which is implemented in
      CNUCTRAN.

 */

#ifndef SOLVER_H
#define SOLVER_H


#include <iostream>
#include <mpfr.h>
#include <mpir.h>
#include <mpreal.h>
#include <smatrix.h>
#include <cnuctran.h>
#include <chrono>
#include <map>
#include <unordered_map>

using namespace std;
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
        const mpreal __two__ = mpreal("2.0");
        const mpreal __one__ = mpreal("1.0");
        const mpreal __neg__ = mpreal("-1.0");
        const mpreal __zer__ = mpreal("0.0");
        vector<string> species_names;
        int __I__;
        vector<vector<mpreal>> lambdas;
        vector<vector<vector<int>>> G;
        vector<vector<mpreal>> fission_yields;

        solver(vector<string> species_names)
        {
            this->species_names = species_names;
            this->__I__ = this->species_names.size();
            for (int i = 0; i < this->__I__; i++)
            {
                this->lambdas.push_back(vector<mpreal>());
                vector<int> tmp1 = { __nop__ }; vector<vector<int>> tmp2; tmp2.push_back(tmp1);
                this->G.push_back(tmp2);
                this->fission_yields.push_back(vector<mpreal>());
            }
            return;
        }

        void add_removal(int species_index,
            mpreal rate,
            vector<int> products,
            vector<mpreal> fission_yields = vector<mpreal>({}))
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
                    vector<mpreal> tmp = vector<mpreal>();
                    for (mpreal y : fission_yields)
                        tmp.push_back(y);
                    this->fission_yields[species_index] = tmp;
                }
                else
                {
                    cout << "FATAL-ERROR <cnuctran::solver::add_removal(...)> Insufficient fission yields given for species " <<
                        this->species_names[species_index] << " products.";
                    exit(1);
                }
            }
            else if (fission_yields.empty() && products.size() == 1)
            {
                for (int product : products)
                {
                    // This part is for preparing the transmutation matrix, which is not relevant 
                    // for this C++ version, cnuctran. For CRAM calculation, please use the PyNUCTRAN.
                }
            }
            else
            {
                cout << "FATAL-ERROR <cnuctran::solver::add_removal(...)> Invalid removal definition for isotope " <<
                    this->species_names[species_index] << endl;
                cout << "Non-fission events MUST only have ONE daughter product." << endl;
                cout << "Whereas fission events MUST have >1 products to track." << endl;
                exit(1);

            }
        }

        smatrix prepare_transfer_matrix(mpreal dt)
        {
            map_2d A;
            map_2d P;
            int i;

            unordered_map<int, mpreal, modified_hash> e;
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
                            n_daughters > 1 ? A[k][i] += a * fission_yields[i][l] :
                                A[k][i] += a;
                        }
                    }

                    if (j == 0) A[i][i] += a;

                }
            }


            return smatrix({ this->__I__, this->__I__ }, A);
        }

        map<string, mpreal> solve(map<string, mpreal> w0,
            mpreal n,
            mpreal t)
        {
            map_2d w0_matrix;
            for (int i = 0; i < this->__I__; i++)
                if (w0.count(this->species_names[i]) == 1)
                    w0_matrix[i][0] = w0[this->species_names[i]];
            smatrix converted_w0 = smatrix(pair<int, int>(this->__I__, 1), w0_matrix);

//..........Auto suggest the no. of substeps.
            int k = int(floor(log(t / pow(mpreal("10"), -n)) / log(__two__)));
            mpreal suggested_substeps = pow(__two__, k);
            mpz_t substeps; mpz_set_str(substeps, suggested_substeps.toString().c_str(), 10);

//..........Compute the transfer matrix power.
            if (__vbs__) cout << "time-step = " << t << " sec(s)."  << endl;
            if (__vbs__) cout << "max-rate = " << __mxr__ << " per sec.\tmin-rate = " << __mnr__ << " per sec." << endl;
            auto t1 = chrono::high_resolution_clock::now();
            smatrix T = this->prepare_transfer_matrix(t / suggested_substeps);
            auto t2 = chrono::high_resolution_clock::now();
            smatrix w = T.binpow(k).mul(converted_w0);
            auto t3 = chrono::high_resolution_clock::now();
            if (__vbs__) cout << "Done computing concentrations.";
            if (__vbs__) cout << chrono::duration_cast<chrono::milliseconds>(t3 - t1).count() << "ms. (" <<
                chrono::duration_cast<chrono::milliseconds>(t3 - t2).count() << "ms. for " << k << " mults.)" << endl;

            map<string, mpreal> out;
            for (int i = 0; i < this->__I__; i++)
                out[this->species_names[i]] = w.nzel[i][0];

            return out;
        }
    };
}

#endif
