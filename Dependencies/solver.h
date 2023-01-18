/*

      This file is part of the CNUCTRAN library

      @author   M. R. Omar (rabieomar@usm.my)
      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran

      Copyright (c) 2023, Universiti Sains Malaysia.

      This header file contains the definitions of all reusable computational 
      routines of the proposed probabilistic method, which is implemented in
      CNUCTRAN.

 */

#ifndef SOLVER_H
#define SOLVER_H

#include <mpreal.h>
#include <smatrix.h>
#include <cnuctran.h>
#include <map>

using namespace std;
using namespace mpfr;
using namespace concurrency;

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

        // lambdas is a 2D vector storing the transmutation constants for all nuclides.
        // Its row corresponds to the various nuclide and its column corresponds to the various removal types.
        vector<vector<mpreal>> lambdas;

        // G is a 2D vector storing the transmutation products ID. The ID maps to the nuclide stored in species_names.
        vector<vector<vector<int>>> G;
        vector<vector<mpreal>> fission_yields;

        solver(vector<string> species_names)
        {
            this->species_names = species_names;
            this->__I__ = this->species_names.size();
            for (int i = 0; i < this->__I__; i++)
            {
                this->lambdas.push_back(vector<mpreal>());
                vector<int> tmp1; 
                tmp1.push_back(__nop__);
                vector<vector<int>> tmp2; tmp2.push_back(tmp1);
                this->G.push_back(tmp2);
                this->fission_yields.push_back(vector<mpreal>());
            }
            return;
        }

        /*
            ADD_REMOVAL
            This subroutine defines the removal event of a nuclide.
        */
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
                    cout << "fatal-error <cnuctran.solver.add_removal(...)>\nInsufficient fission yields given for species " <<
                        this->species_names[species_index] << " products.";
                    exit(1);
                }
            }
            else if (fission_yields.empty() && products.size() == 1)
            {
                for (int product : products)
                {
                    // This part is for preparing the transmutation matrix, which is not relevant 
                    // for this C++ version, cnuctran.
                }
            }
            else
            {
                cout << "fatal-error <cnuctran.solver.add_removal(...)>\nInvalid removal definition for isotope " <<
                    this->species_names[species_index] << endl;
                cout << "Non-fission events MUST only have ONE daughter product." << endl;
                cout << "Whereas fission events MUST have >1 products to track." << endl;
                exit(1);

            }
        }

        /*
            PREPARE_TRANSFER_MATRIX
            This function returns the transfer matrix, P, in Eq. (17) of CNUCTRAN manual.
        */
        smatrix prepare_transfer_matrix(mpreal dt)
        {
            cmap_2d A;
            cmap_2d P;
            int i;

            unordered_map<int, mpreal> e;
            
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

//..............Constructs the pi-distribution according to Eq. (12) if CNUCTRAN manual.
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

//..............Constructs the transfer matrix according to Eq. (15) of CNUCTRAN manual.
                auto const& gI = G[i];
                for (int j = 0; j < n_events; j++)
                {

                    auto const& a = (P[i][j] / norm);
                    auto const& gJ = gI[j];
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


        /*
            SOLVE
            This function solves the final nuclides concentration according to Eq. (18) of CNUCTRAN manual.
        */
        map<string, mpreal> solve(map<string, mpreal> w0,
            mpreal n,
            mpreal t)
        {
            cmap_2d w0_matrix;
            for (int i = 0; i < this->__I__; i++)
                if (w0.count(this->species_names[i]) == 1)
                    w0_matrix[i][0] = w0[this->species_names[i]];
                    
            smatrix converted_w0 = smatrix(pair<int, int>(this->__I__, 1), w0_matrix);

            //..........Auto suggest the no. of Sparse Self Matrix Multiplication.
            int k = int(floor(log(t / pow(mpreal("10"), -n)) / log(__two__)));
            if (__vbs__) cout << "Approximation order, n = " << n << endl;

            //..........Compute the transfer matrix power.
            if (__vbs__) cout << "Time step, T = " << t << endl;
            auto t1 = chrono::high_resolution_clock::now();
            smatrix T = this->prepare_transfer_matrix(t / pow(__two__, k));

            auto t2 = chrono::high_resolution_clock::now();
            //..........Compute the matrix exponentiation and multiply with w0 to obtain w.
            T.binpow(k);
            smatrix w = T.mul(converted_w0);
            auto t3 = chrono::high_resolution_clock::now();
            if (__vbs__) cout << "Done computing concentrations. ";
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
