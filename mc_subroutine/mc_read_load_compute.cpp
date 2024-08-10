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
//    std::random_device rd;
//    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    std::chrono::duration<double> totalTime_propose_uni(0);
    std::chrono::duration<double> totalTime_UCurr(0);
    std::chrono::duration<double> totalTime_acceptanceRatio_uni(0);
    std::chrono::duration<double> totalTime_acc_copy(0);
    std::chrono::duration<double> totalTime_copy(0);
    std::chrono::duration<double> totalTime_save(0);

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

//            this->proposal(xAVecCurr,xBVecCurr,xAVecNext,xBVecNext);

            //propose next
            auto startProposalUni = std::chrono::steady_clock::now();
            this->proposal_uni(xAVecCurr,xBVecCurr,xAVecNext,xBVecNext);
            auto endProposalUni = std::chrono::steady_clock::now();
            totalTime_propose_uni += endProposalUni - startProposalUni;


            //next U
            auto startUCurr = std::chrono::steady_clock::now();
            double UNext;
            UCurr=(*potFuncPtr)(xAVecCurr.get(),xBVecCurr.get());
            auto endUCurr = std::chrono::steady_clock::now();
            totalTime_UCurr+=endUCurr-startUCurr;


            //accept reject
            auto start_acc=std::chrono::steady_clock::now();
//            double r= acceptanceRatio(xAVecCurr,xBVecCurr,UCurr,
//                                      xAVecNext,xBVecNext,UNext);
            double r=this->acceptanceRatio_uni(xAVecCurr,xBVecCurr,UCurr,
                                      xAVecNext,xBVecNext,UNext);
            auto end_acc=std::chrono::steady_clock::now();
            totalTime_acceptanceRatio_uni+=end_acc-start_acc;


            //accept-copy
            auto acc_copyStart=std::chrono::steady_clock::now();
            double u = distUnif01(e2);
            if (u <= r) {

                std::memcpy(xAVecCurr.get(),xAVecNext.get(),N*sizeof (double ));
                std::memcpy(xBVecCurr.get(),xBVecNext.get(),N*sizeof (double ));


                UCurr = UNext;

            }//end of accept-reject
            auto acc_copyEnd=std::chrono::steady_clock::now();
            totalTime_acc_copy+=acc_copyEnd-acc_copyStart;


            //copy and save
            auto copyStart=std::chrono::steady_clock::now();
            U_dist_ptr[varNum*j+0]=UCurr;

            for(int n=1;n<1+N;n++){
                U_dist_ptr[varNum*j+n]=xAVecCurr[n-1];
            }
            for(int n=N+1;n<=2*N;n++){
                U_dist_ptr[varNum*j+n]=xBVecCurr[n-N-1];
            }
            auto copyEnd=std::chrono::steady_clock::now();
            totalTime_copy+=copyEnd-copyStart;

        }//end for loop

        //save to csv

        auto saveStart=std::chrono::steady_clock::now();
        size_t loopEnd = loopStart + loopToWrite - 1;
        std::string fileNameMiddle = "loopStart" + std::to_string(loopStart) + "loopEnd" + std::to_string(loopEnd);
//        std::string out_U_distPickleFileName = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.csv";
        std::string out_U_distPickleFileName_pkl = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.pkl";
                std::string out_U_distPickleFileName_csv = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.csv";

        //save U_dist_ptr
//        saveArrayToCSV(U_dist_ptr,varNum * loopToWrite,out_U_distPickleFileName,varNum);
        save_array_to_pickle(U_dist_ptr.get(),varNum * loopToWrite,out_U_distPickleFileName_pkl);
        saveLastData2Csv(U_dist_ptr,varNum * loopToWrite,out_U_distPickleFileName_csv,varNum);
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "loop " + std::to_string(loopStart) + " to loop " + std::to_string(loopEnd) + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

        loopStart = loopEnd + 1;
        auto saveEnd=std::chrono::steady_clock::now();
        totalTime_save+=saveEnd-saveStart;

    }//end flush for loop

    std::cout<<"mc executed for "<<flushNum<<" flushes."<<std::endl;

    std::cout<<"totalTime_propose_uni is "<<totalTime_propose_uni.count()<<" s\n";

    std::cout<<"totalTime_UCurr is "<<totalTime_UCurr.count()<<" s\n";

    std::cout<<"totalTime_acceptanceRatio_uni is "<<totalTime_acceptanceRatio_uni.count()<<" s\n";

    std::cout<<"totalTime_acc_copy is "<<totalTime_acc_copy.count()<<" s\n";

    std::cout<<"totalTime_copy is "<<totalTime_copy.count()<<" s\n";
    std::cout<<"totalTime_save is "<<totalTime_save.count()<<" s\n";
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
//    std::cout<<std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;
    UNext=(*potFuncPtr)(xAVecNext.get(),xBVecNext.get());
//    std::cout<<"UNext="<<UNext<<std::endl;
    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);
//    std::cout<<"R="<<R<<std::endl;





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
//    std::cout<<"R="<<R<<std::endl;

    for(int i=0;i<N;i++){
        double zxBOneCurrVal= zVal(xBVecCurr[i],0,lm);
        double zxBOneNextVal= zVal(xBVecNext[i],0,lm);
        double ratio_xBOneVal=zxBOneCurrVal/zxBOneNextVal;
        R*=ratio_xBOneVal;
    }

    return std::min(1.0,R);



}


///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
double mc_computation::generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps){


double xMinusEps=x-eps;
double xPlusEps=x+eps;

double unif_left_end=xMinusEps<leftEnd?leftEnd:xMinusEps;
double unif_right_end=xPlusEps>rightEnd?rightEnd:xPlusEps;

//    std::random_device rd;
//    std::ranlux24_base e2(rd());

double unif_left_end_double_on_the_right=std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());



    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right,unif_right_end); //[unif_left_end_double_on_the_right, unif_right_end)

    double xNext=distUnif(e2);
    return xNext;



}



///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon){

    if (a<y and y<a+epsilon){
        return 1.0/(y-a+epsilon);
    } else if( a+epsilon<=y and y<b+epsilon){
        return 1.0/(2.0*epsilon);
    }else if(b-epsilon<=y and y<b){
        return 1/(b-y+epsilon);
    } else{

        std::cerr<<"value out of range."<<std::endl;
        std::exit(10);


    }


}

///
/// @param xAVecCurr
/// @param xBVecCurr
/// @param xAVecNext
/// @param xBVecNext
void mc_computation::proposal_uni(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr,
                   std::shared_ptr<double[]>& xAVecNext, std::shared_ptr<double[]>& xBVecNext){




    ////////////////////////////////////////////////////////////
    /// serial code
    double lm=potFuncPtr->getLm();
    // Parallelize the loop for xAVecNext
//    #pragma omp parallel for
    for(int i=0;i<N;i++){

        xAVecNext[i]= generate_uni_open_interval(xAVecCurr[i],0,lm,h);
    }
    // Parallelize the loop for xBVecNext
//    #pragma omp parallel for
    for(int i=0;i<N;i++){
        xBVecNext[i]=generate_uni_open_interval(xBVecCurr[i],0,lm,h);
    }
    ///end serial code
    ////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////
    //parallel code
//    double lm = potFuncPtr->getLm();
//
//    // Define a lambda function for concurrent execution
//    auto generate_for_range = [&](int start, int end, const std::shared_ptr<double[]>& xCurr, std::shared_ptr<double[]>& xNext) {
//        for (int i = start; i < end; ++i) {
//            xNext[i] = generate_uni_open_interval(xCurr[i], 0, lm, h);
//        }
//    };
//    // Determine the number of threads to use
//    unsigned int numThreads = std::thread::hardware_concurrency();
//    int chunkSize = N / numThreads;
//    std::vector<std::thread> threads;
//    // Launch threads to process xAVecNext
//    for (unsigned int t = 0; t < numThreads; ++t) {
//        int start = t * chunkSize;
//        int end = (t == numThreads - 1) ? N : start + chunkSize;
//        threads.emplace_back(generate_for_range, start, end, std::ref(xAVecCurr), std::ref(xAVecNext));    }
//    // Join threads for xAVecNext
//    for (auto& thread : threads) {
//        thread.join();
//    }
//    // Clear the threads vector and reuse it for xBVecNext
//    threads.clear();
//
//    // Launch threads to process xBVecNext
//    for (unsigned int t = 0; t < numThreads; ++t) {
//        int start = t * chunkSize;
//        int end = (t == numThreads - 1) ? N : start + chunkSize;
//        threads.emplace_back(generate_for_range, start, end, std::ref(xBVecCurr), std::ref(xBVecNext));    }
//    // Join threads for xBVecNext
//    for (auto& thread : threads) {
//        thread.join();
//    }
    //end parallel code
    ////////////////////////////////////////////////////////////////

}



///
/// @param xAVecCurr
/// @param xBVecCurr
/// @param UCurr
/// @param xAVecNext
/// @param xBVecNext
/// @param UNext
/// @return
double mc_computation::acceptanceRatio_uni(const std::shared_ptr<double[]>& xAVecCurr ,const std::shared_ptr<double[]>&xBVecCurr
        ,const double& UCurr,
                           const std::shared_ptr<double[]>& xAVecNext,
                           const std::shared_ptr<double[]>&  xBVecNext,
                           double &UNext){


    double lm=potFuncPtr->getLm();
    UNext=(*potFuncPtr)(xAVecNext.get(),xBVecNext.get());
    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);

    for(int i=0;i<N;i++){
    double S_xACurrNext= S_uni(xAVecCurr[i],xAVecNext[i],0,lm,h);
    double S_xANextCurr= S_uni(xAVecNext[i],xAVecCurr[i],0,lm,h);
    double ratio_xAi=S_xACurrNext/S_xANextCurr;
        if (std::fetestexcept(FE_DIVBYZERO)) {
            std::cout << "Division by zero exception caught." << std::endl;
            std::exit(15);
        }

        if (std::isnan(ratio_xAi)) {
            std::cout << "The result is NaN." << std::endl;
            std::exit(15);
        }
    R*=ratio_xAi;

    }

    for(int i=0;i<N;i++){
        double S_xBCurrNext= S_uni(xBVecCurr[i],xBVecNext[i],0,lm,h);
        double S_xBNextCurr= S_uni(xBVecNext[i],xBVecCurr[i],0,lm,h);
        double ratio_xBi=S_xBCurrNext/S_xBNextCurr;
        if (std::fetestexcept(FE_DIVBYZERO)) {
            std::cout << "Division by zero exception caught." << std::endl;
            std::exit(15);
        }

        if (std::isnan(ratio_xBi)) {
            std::cout << "The result is NaN." << std::endl;
            std::exit(15);
        }
        R*=ratio_xBi;

    }
    return std::min(1.0,R);


}


void mc_computation::save_array_to_pickle(double *ptr, std::size_t size, const std::string& filename){
    using namespace boost::python;
    try {
        Py_Initialize();  // Initialize the Python interpreter
        if (!Py_IsInitialized()) {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }

        // Debug output
        std::cout << "Python interpreter initialized successfully." << std::endl;

        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Create a Python list from the C++ array
        list py_list;
        for (std::size_t i = 0; i < size; ++i) {
            py_list.append(ptr[i]);
        }

        // Serialize the list using pickle.dumps
        object serialized_array = pickle_dumps(py_list);

        // Extract the serialized data as a string
        std::string serialized_str = extract<std::string>(serialized_array);

        // Write the serialized data to a file
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open file for writing");
        }
        file.write(serialized_str.data(), serialized_str.size());
        file.close();

        // Debug output
        std::cout << "Array serialized and written to file successfully." << std::endl;
    } catch (const error_already_set&) {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    if (Py_IsInitialized()) {
        Py_Finalize();  // Finalize the Python interpreter
    }



}


void mc_computation::saveLastData2Csv(const std::shared_ptr<double[]>& array, const  int& arraySize, const std::string& filename, const int& numbersPerRow){
//saves last row to csv
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
    for(int i=arraySize-numbersPerRow;i<arraySize;i++){
        outFile<<array[i];
        if ((i + 1) % numbersPerRow == 0) {
            outFile << '\n';
        } else {
            outFile << ',';
        }


    }
//    outFile << '\n';
    outFile.close();

}
