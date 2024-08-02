#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>
#include "potentialFunction/potentialFunctionPrototype.hpp"


int main(int argc, char *argv[]) {

    int N=2;
    std::shared_ptr<double[]> xA=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
    std::shared_ptr<double[]> xB=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
    xA[0]=1;
    xA[1]=3;
    xB[0]=2;
    xB[1]=4;

    auto funcPtr= createPotentialFunction("V_inv_12_6","25,80,15,67,2");
    funcPtr->init();
    double val=(*funcPtr)(xA.get(),xB.get());
    std::cout<<"val="<<val<<std::endl;
    return 0;
}