/*
      This file is part of the CNUCTRAN source code.
      
      @author   M. R. Omar (rabieomar@usm.my)
      @license  MIT
      @link     https://github.com/rabieomar92/cnuctran
      
      Copyright (c) 2023, Universiti Sains Malaysia.
      
      This is the MAIN program body. 
      
 */

#include <simulation.h>
#include <iostream>

using namespace cnuctran;

int main(int nargs, char** argv)
{

    nargs == 2? 
        simulation::from_input(argv[1]) :
        simulation::from_input(".\\input.xml");

    return 0;

}

