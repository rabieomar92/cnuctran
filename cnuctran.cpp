// cnuctran.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <simulation.h>
#include <iostream>

using namespace cnuctran;
using namespace std;

int main(int nargs, char** argv)
{
    nargs == 2? 
        simulation::from_input(argv[1]) :
        simulation::from_input(".\\input.xml");
    return 0;

}