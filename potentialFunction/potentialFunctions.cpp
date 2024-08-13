//
// Created by polya on 7/19/24.
//

#include "potentialFunctionPrototype.hpp"

class V_inv_12_6:public potentialFunction {

public:
    V_inv_12_6(const std::string &coefsStr):potentialFunction(){
        this->coefsInStr=coefsStr;
    }

public:
    void json2Coefs(const std::string &coefsStr)override{
        std::stringstream iss;
        iss<<coefsStr;
        std::string temp;
        //read a1
        if (std::getline(iss, temp, ',')){
            this->a1=std::stod(temp);
        }

        //read b1

        if (std::getline(iss, temp, ',')){
            this->b1=std::stod(temp);
        }

        //read a2
        if (std::getline(iss, temp, ',')){
            this->a2=std::stod(temp);
        }

        //read b2

        if (std::getline(iss, temp, ',')){
            this->b2=std::stod(temp);
        }

        //read N

        if (std::getline(iss, temp, ',')){
            this->N=std::stoi(temp);
        }
    }









    void init() override{
        this->json2Coefs(coefsInStr);
        this->r1=std::pow(2.0*a1/b1,1.0/6.0);
        this->r2=std::pow(2.0*a2/b2,1.0/6.0);
        this->lm=(static_cast<double >(N)*(r1+r2))*1.5;
        this->eps=((r1+r2)/2.0)/8;
//        pow_result_tmp=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

        std::cout << "a1=" << a1 << ", b1=" << b1 << ", a2=" << a2 << ", b2=" << b2 << std::endl;
        std::cout<<"r1="<<r1<<", r2="<<r2<<std::endl;
        std::cout<<"lm="<<lm<<std::endl;
//        std::cout<<"eps="<<eps<<std::endl;


    }

    ///
    /// @param L
    /// @param xAVec
    /// @param xBVec
    /// @return
    double operator()(const double *xAVec,const double *xBVec) override {

      double val=0;
      for(int j=0;j<N;j++){
          double d0j=xBVec[j]-xAVec[j];
          val+= V1(d0j);
      }

      for(int j=0;j<N-1;j++){
          double d1j=xAVec[j+1]-xBVec[j];
          val+= V2(d1j);
      }

        return val;


    }//end of () operator

    double V1(const double &r){
        double  val=a1*std::pow(r,-12.0)-b1*std::pow(r,-6.0);
        return val;
    }

    double V2(const double &r){
        double val=a2*std::pow(r,-12.0)-b2*std::pow(r,-6.0);
        return val;
    }

    double getLm() const override {
        return lm;
    }
    double get_eps() const override {
        return eps;
    }
public:
    double a1;
    double a2;
    double b1;
    double b2;
    std::string coefsInStr;
    double r1;//min position of V1
    double r2;//min position of V2
    double lm;//range of distances
    double eps;//half interval length of uniform distribution
    int N;
    std::shared_ptr<double[]>pow_result_tmp;
};

std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &coefsJsonStr) {
    if (funcName == "V_inv_12_6") {

        return std::make_shared<V_inv_12_6>(coefsJsonStr);
    }

    else {
        throw std::invalid_argument("Unknown potential function type");
    }
}