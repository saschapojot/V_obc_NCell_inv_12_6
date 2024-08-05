//
// Created by polya on 7/19/24.
//

#include "mc_read_load_compute.hpp"





void mc_computation::execute_mc(const std::shared_ptr<double[]>& xAVec, const std::shared_ptr<double[]>& xBVec, const size_t & loopInit, const size_t & flushNum){


    std::shared_ptr<double[]> xAVecCurr=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

    std::shared_ptr<double[]> xBVecCurr=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

    std::memcpy(xAVecCurr.get(),xAVec.get(),N*sizeof (double ));
    std::memcpy(xBVecCurr.get(),xBVec.get(),N*sizeof (double ));

//    std::cout<<"d0VecCurr: ";
//    print_shared_ptr(d0VecCurr,N);
//    std::cout<<"d1VecCurr: ";
//    print_shared_ptr(d1VecCurr,N-1);

    //initialize next values

    std::shared_ptr<double[]> xAVecNext=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
    std::shared_ptr<double[]> xBVecNext=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

    double UCurr;
    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    size_t loopStart = loopInit;
    for (size_t fls = 0; fls < flushNum; fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (size_t j = 0; j < loopToWrite; j++) {
            //propose a move
//            double LNext;
//            double y0Next;
//            double z0Next;
//            double y1Next;
//            double LReset;

            this->proposal(xAVecCurr,xBVecCurr,xAVecNext,xBVecNext);
            double UNext;
            UCurr=(*potFuncPtr)(xAVecCurr.get(),xBVecCurr.get());
            double r= acceptanceRatio(xAVecCurr,xBVecCurr,UCurr,
                                      xAVecNext,xBVecNext,UNext);
            double u = distUnif01(e2);
            if (u <= r) {

                std::memcpy(xAVecCurr.get(),xAVecNext.get(),N*sizeof (double ));
                std::memcpy(xBVecCurr.get(),xBVecNext.get(),N*sizeof (double ));


                UCurr = UNext;

            }//end of accept-reject
            U_dist_ptr[varNum*j+0]=UCurr;

            for(int n=1;n<1+N;n++){
                U_dist_ptr[varNum*j+n]=xAVecCurr[n-1];
            }
            for(int n=N+1;n<=2*N;n++){
                U_dist_ptr[varNum*j+n]=xBVecCurr[n-N-1];
            }

        }//end for loop
        size_t loopEnd = loopStart + loopToWrite - 1;
        std::string fileNameMiddle = "loopStart" + std::to_string(loopStart) + "loopEnd" + std::to_string(loopEnd);
        std::string out_U_distPickleFileName = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.csv";

        //save U_dist_ptr
        saveArrayToCSV(U_dist_ptr,varNum * loopToWrite,out_U_distPickleFileName,varNum);
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "loop " + std::to_string(loopStart) + " to loop " + std::to_string(loopEnd) + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

        loopStart = loopEnd + 1;
    }//end flush for loop

    std::cout<<"mc executed for "<<flushNum<<" flushes."<<std::endl;


}







void mc_computation::saveArrayToCSV(const std::shared_ptr<double[]>& array, const  int& arraySize, const std::string& filename, const int& numbersPerRow) {

    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

    outFile<<"U";
    for(int i=0;i<N;i++){
        outFile<<",x"+std::to_string(i)+"A";
    }
    for(int i=0;i<N;i++){
        outFile<<",x"+std::to_string(i)+"B";
    }
    outFile<<"\n";
    for (int i = 0; i < arraySize; ++i) {
        outFile << array[i];
        if ((i + 1) % numbersPerRow == 0) {
            outFile << '\n';
        } else {
            outFile << ',';
        }
    }

    // If the last row isn't complete, end the line
    if (arraySize % numbersPerRow != 0) {
        outFile << '\n';
    }

    outFile.close();


}

void mc_computation::init_and_run(){
    this->execute_mc(xAVecInit,xBVecInit,loopLastFile+1,newFlushNum);


}



///
/// @param y
/// @param x center
/// @return known proposal function, which is normal distribution
double mc_computation::Q(const double &y, const double &x, const double &a, const double &b){

    double val=1/(std::pow(2.0*PI,0.5)*h)
               *std::exp(-1/(2*std::pow(h,2))*std::pow(y-x,2.0));

    return val;

}


///
/// @param y
/// @param x center
/// @param a left end
/// @param b right end
/// @return truncated Gaussian
double mc_computation::f(const double &y, const double &x, const double &a, const double &b){


    if(y<=a or y>=b){
        return 0;
    }else{

        double val=std::exp(-1.0/(2.0*std::pow(h,2))*std::pow(y-x,2));
        return val;
    }

}


///
/// @param x center
/// @param a left end
/// @param b right end
/// @return random number from truncated Gaussian
double mc_computation::reject_sampling_one_data(const double &x,const double &a, const double &b){

    std::random_device rd;  // Create a random device object
    std::ranlux24_base engine(rd());  // Seed the engine with the random device

    std::normal_distribution<> normal_dist(x,h);
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    double y=normal_dist(engine);
    double u=distUnif01(engine);

    while(u>=f(y,x,a,b)/(M* Q(y,x,a,b))){
        y=normal_dist(engine);
        u=distUnif01(engine);

    }

    return y;

}

///
/// @param x center
/// @param a left end
/// @param b right end
/// @return integral
double mc_computation::zVal(const double& x,const double &a, const double &b){

    auto integrandWithParam=[x,a,b, this](const double &y){
        return this->integrand(y,x,a,b);
    };
    double tolerance = 1e-20;
    boost::math::quadrature::gauss_kronrod<double, 80> integrator;
//    double result = boost::math::quadrature::trapezoidal(integrandWithParam,a,b, tolerance);
    double result = integrator.integrate(integrandWithParam,a,b,tolerance);
    return result;

}

///
/// @param y
/// @param x center
/// @param a left end
/// @param b right end
/// @return
double mc_computation::integrand(const double &y, const double& x,const double &a, const double &b){

    return f(y,x,a,b);


}

///
/// @param xAVecCurr
/// @param xBVecCurr
/// @param xAVecNext
/// @param xBVecNext
void mc_computation::proposal(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr,
                              std::shared_ptr<double[]>& xAVecNext, std::shared_ptr<double[]>& xBVecNext){


    //proposal using truncated Gaussian
    double lm = potFuncPtr->getLm();
    for(int i=0;i<N;i++){
        xAVecNext[i]=reject_sampling_one_data(xAVecCurr[i],0,lm);
    }

    for(int i=0;i<N;i++){
        xBVecNext[i]=reject_sampling_one_data(xBVecCurr[i],0,lm);
    }


}


///
/// @param LCurr
/// @param d0VecCurr
/// @param d1VecCurr
/// @param UCurr
/// @param LNext
/// @param d0VecNext
/// @param d1VecNext
/// @param UNext
/// @return
double mc_computation::acceptanceRatio(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr
        ,const double& UCurr,
        const std::shared_ptr<double[]>& xAVecNext,
        const std::shared_ptr<double[]>&  xBVecNext,
        double &UNext){

    double lm=potFuncPtr->getLm();
    std::cout<<std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;
    UNext=(*potFuncPtr)(xAVecNext.get(),xBVecNext.get());
    std::cout<<"UNext="<<UNext<<std::endl;
    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);
    std::cout<<"R="<<R<<std::endl;





    for(int i=0;i<N;i++){
    double zxAOneCurrVal= zVal(xAVecCurr[i],0,lm);
//    std::cout<<"xAVecCurr[i]="<<xAVecCurr[i]<<std::endl;
//    std::cout<<"zxAOneCurrVal="<<zxAOneCurrVal<<std::endl;
    double zxAOneNextVal= zVal(xAVecNext[i],0,lm);
//    std::cout<<"xAVecNext[i]="<<xAVecNext[i]<<std::endl;
//    std::cout<<"zxAOneNextVal="<<zxAOneNextVal<<std::endl;
    double ratio_xAOneVal=zxAOneCurrVal/zxAOneNextVal;
//    std::cout<<"ratio_xAOneVal="<<ratio_xAOneVal<<std::endl;
    R*=ratio_xAOneVal;

    }
    std::cout<<"R="<<R<<std::endl;

    for(int i=0;i<N;i++){
        double zxBOneCurrVal= zVal(xBVecCurr[i],0,lm);
        double zxBOneNextVal= zVal(xBVecNext[i],0,lm);
        double ratio_xBOneVal=zxBOneCurrVal/zxBOneNextVal;
        R*=ratio_xBOneVal;
    }

    return std::min(1.0,R);



}
