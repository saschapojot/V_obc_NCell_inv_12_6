//
// Created by polya on 7/19/24.
//

#ifndef V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
#define V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
#include "../potentialFunction/potentialFunctionPrototype.hpp"
#include <boost/filesystem.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
//#include <boost/math/special_functions/next.hpp>
#include <boost/python.hpp>
#include <boost/python/object/pickle_support.hpp>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
//#include <boost/random/random_device.hpp>
//#include <boost/random/ranlux.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
//#include <boost/serialization/vector.hpp>
#include <cfenv> // for floating-point exceptions
#include <chrono>
#include <cstdlib>
#include <cxxabi.h>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <thread>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;
const auto PI=M_PI;


class mc_computation {
public:
    mc_computation(const std::string &cppInParamsFileName): e2(std::random_device{}())  {
        std::ifstream file(cppInParamsFileName);
        if (!file.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }
        std::string line;

        int paramCounter = 0;

        while (std::getline(file, line)) {
            // Check if the line is empty
            if (line.empty()) {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);

            //read T
            if (paramCounter == 0) {
                iss >> T;
                if (T <= 0) {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                }//end if
                std::cout << "T=" << T << std::endl;
                this->beta = 1 / T;
//                double stepForT1 = 0.1;

                double h_threshhold=0.004;

                this->h=h_threshhold;
//                this->h = stepForT1 * T > h_threshhold ? h_threshhold : stepForT1 * T;//stepSize;
                std::cout << "h=" << h << std::endl;
                this->M = std::pow(2.0 * PI, 0.5) * h * 1.001;
                std::cout<<"M="<<M<<std::endl;
                paramCounter++;
                continue;
            }//end reading T

            //read unit cell number
            if (paramCounter==1) {

                iss >> N;

                xAVecInit = std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
                xBVecInit = std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

                paramCounter++;
                continue;

            }

            //read coefficients
            if(paramCounter==2){
                iss>>coefsToPotFunc;
                paramCounter++;
                continue;

            }// end reading coefficients


            //read potential function name
            if(paramCounter==3){
                iss>>potFuncName;
                paramCounter++;
                continue;
            }//end reading potential function name

            //read initial values
            if(paramCounter==4){
               std::vector<double> values;
                std::string valueStr;

                while (std::getline(iss, valueStr, ',')) {
                    values.push_back(std::stod(valueStr));
                }

               for(int i=0;i<N;i++){
                   xAVecInit[i]=values[i];
               }
               for(int i=N;i<2*N;i++){
                   xBVecInit[i-N]=values[i];
               }





                paramCounter++;
                continue;




            }//end reading initial values


            //read loopToWrite
            if(paramCounter==5){
                iss>>loopToWrite;
                paramCounter++;
                continue;
            }//end reading loopToWrite

            //read newFlushNum
            if(paramCounter==6){
                iss>>newFlushNum;
                paramCounter++;
                continue;
            }//end reading newFlushNum

            //read loopLastFile
            if(paramCounter==7){
                //if loopLastFileStr is "-1", loopLastFile uses the overflowed value
                //and loopLastFile+1 will be 0
                iss>>loopLastFile;
                paramCounter++;
                continue;
            }//end reading loopLastFile

            //read TDirRoot
            if (paramCounter==8){
                iss>>TDirRoot;
                paramCounter++;
                continue;
            }//end reading TDirRoot

            //read U_dist_dataDir
            if(paramCounter==9){
                iss>>U_dist_dataDir;
                paramCounter++;
                continue;
            }//end reading U_dist_dataDir



        }//end while

//        std::cout<<"LInit="<<LInit<<std::endl;
        std::cout<<"unit cell number="<<N<<std::endl;
        std::cout<<"xAVecInit: \n";
        print_shared_ptr(xAVecInit,N);
        std::cout<<"xBVecInit: \n";
        print_shared_ptr(xBVecInit,N);
        this->potFuncPtr = createPotentialFunction(potFuncName, coefsToPotFunc);
        potFuncPtr->init();
        this->varNum = 2*N+1;//U,xA, xB
        try {
            this->U_dist_ptr= std::shared_ptr<double[]>(new double[loopToWrite * varNum],
                                                        std::default_delete<double[]>());
        }
        catch (const std::bad_alloc &e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception &e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }


        std::cout<<"loopToWrite="<<loopToWrite<<std::endl;
        std::cout<<"newFlushNum="<<newFlushNum<<std::endl;
        std::cout<<"loopLastFile+1="<<loopLastFile+1<<std::endl;
        std::cout<<"TDirRoot="<<TDirRoot<<std::endl;
        std::cout<<"U_dist_dataDir="<<U_dist_dataDir<<std::endl;

    }//end constructor




public:

        ///
    /// @param x
    /// @param leftEnd
    /// @param rightEnd
    /// @param eps
    /// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
   double generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps);

        ///
    /// @param x proposed value
    /// @param y current value
    /// @param a left end of interval
    /// @param b right end of interval
    /// @param epsilon half length
    /// @return proposal probability S(x|y)
    double S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon);

    ///
    /// @param xAVecCurr
    /// @param xBVecCurr
    /// @param xAVecNext
    /// @param xBVecNext
    void proposal_uni(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr,
                       std::shared_ptr<double[]>& xAVecNext, std::shared_ptr<double[]>& xBVecNext);

    ///
    /// @param xAVecCurr
    /// @param xBVecCurr
    /// @param UCurr
    /// @param xAVecNext
    /// @param xBVecNext
    /// @param UNext
    /// @return
    double acceptanceRatio_uni(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr
            ,const double& UCurr,
                               const std::shared_ptr<double[]>& xAVecNext,
                               const std::shared_ptr<double[]>&  xBVecNext,
                               double &UNext);

    ///
    /// @param xAVecCurr
    /// @param xBVecCurr
    /// @param xAVecNext
    /// @param xBVecNext
    void proposal(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr,
                   std::shared_ptr<double[]>& xAVecNext, std::shared_ptr<double[]>& xBVecNext);


   ///
   /// @param xAVecCurr
   /// @param xBVecCurr
   /// @param UCurr
   /// @param xAVecNext
   /// @param xBVecNext
   /// @param UNext
   /// @return
    double acceptanceRatio(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr
                           ,const double& UCurr,
                           const std::shared_ptr<double[]>& xAVecNext,
                           const std::shared_ptr<double[]>&  xBVecNext,
                           double &UNext);

    ///
/// @param y
/// @param x center
/// @param a left end
///@param b right end
/// @return known proposal function, which is normal distribution
    double Q(const double &y, const double &x, const double &a, const double &b);

    ///
    /// @param y
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return truncated Gaussian
    double f(const double &y, const double &x, const double &a, const double &b);

    ///
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return random number from truncated Gaussian
    double reject_sampling_one_data(const double &x,const double &a, const double &b);

    ///
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return integral
    double zVal(const double& x,const double &a, const double &b);

    ///
    /// @param y
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return
    double integrand(const double &y, const double& x,const double &a, const double &b);

    void execute_mc(const std::shared_ptr<double[]>& xAVec, const std::shared_ptr<double[]>& xBVec, const size_t & loopInit, const size_t & flushNum);


     void saveArrayToCSV(const std::shared_ptr<double[]>& array, const  int& arraySize, const std::string& filename, const int& numbersPerRow) ;

    void init_and_run();


    template<class T>
            void print_shared_ptr(const std::shared_ptr<T> &ptr,const int& size){
        if (!ptr) {
            std::cout << "Pointer is null." << std::endl;
            return;
        }

        for(int i=0;i<size;i++){
            if(i<size-1){
                std::cout<<ptr[i]<<",";
            }
            else{
                std::cout<<ptr[i]<<std::endl;
            }
        }

            }


    void save_array_to_pickle(double *ptr, std::size_t size, const std::string& filename);


public:
    double T;// temperature
    double beta;
    double h;// step size
    size_t loopToWrite;
    size_t newFlushNum;
    size_t loopLastFile;
    std::shared_ptr<potentialFunction> potFuncPtr;
    std::string TDirRoot;
    std::string U_dist_dataDir;
    std::shared_ptr<double[]> U_dist_ptr;
    int varNum;
//    double LInit;
    std::shared_ptr<double[]> xAVecInit;
    std::shared_ptr<double[]> xBVecInit;

//    double y0Init;
//    double z0Init;
//    double y1Init;
    std::string coefsToPotFunc;
    std::string potFuncName;
    double M;
    int N;//unit cell number
    std::ranlux24_base e2;

};


#endif //V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
