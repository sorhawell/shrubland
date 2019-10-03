//#include <algorithm>
//#include <string>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <random>
//#include <chrono>
//#include <iomanip>

#include "./shrubland_arduino/lib/dynarray/src/dynarray.h"
#include "./shrubland_arduino/lib/shrub_learner/src/shrub_learner.h"
#include "./shrubland_arduino/lib/testdata/src/testdata.h"
//#include "shrub_learner.h"

using namespace std;

CPUtimer Timer;



constexpr int ntree=50;
constexpr int nnode=50*300;
constexpr int nrow =1000;
constexpr int ncol =6;
constexpr int ycol =5;


forest_lake forlake(uint16_t(12000));
uint16_t n_train=500;
float* Dp_train = DData;
float* Dp_test = DData_test;

dynamic_array<float,uint16_t> X(Dp_train    , n_train*ncol,ncol);
dynamic_array<float,uint16_t> X_test(Dp_test, n_train*ncol,ncol);
dynamic_vector<float,uint16_t> y = X.get_vector(ycol); //
dynamic_vector<float,uint16_t> y_test = X_test.get_vector(ycol); //
dynamic_array<float,uint16_t> y_pred_arr(n_train,1);
dynamic_vector<float,uint16_t> y_pred = y_pred_arr.get_vector(0);

dynamic_array<uint16_t,uint16_t> radArr(10000,1);
dynamic_vector<uint16_t,uint16_t> radVec = radArr.get_vector(0);
    

int n_trials{0};
constexpr int nval = nrow*ncol;
float sum_time{0};
int main() {
    

    //radix test
    /* Timer.setFreq();
    Timer.reset_sumTime();
    iota(radVec.begin(),radVec.end(),0);
    radVec.print();
    for(int i=0;i<10;i++){
    my_shuff(radVec.begin(),radVec.end());
    Timer.start();
    RadixSort(radVec.begin(),radVec.get_size());
    //std::sort(radVec.begin(),radVec.end());
    Timer.stop();
    sprintln(float(Timer.elapsedTime/1000));sprintln(float(Timer.sumTime/1000));
    }
    sprint(" sort:"); sprintln(float(Timer.sumTime/1000));
    radVec.print(); */

    
    sprintln(42);
    X.print();
    X.truncate_col(); //revoke X's read/write access this last column
    X_test.truncate_col(); //revoke X's read/write access this last column

    
     //std::shared_ptr<dynamic_array<double,int>> sp3(new dynamic_array<double,int>(5,10), array_deleter<dynamic_array<double,int>>());
    //std::cout << "\n ncopy:" <<sp3.use_count();
    
    

    forlake.pre_index_features(&X);
    sprint("survived to here");
    sprintln(forlake.index_fixed.at(0,0));
    sprintln(forlake.index_fixed.at(1,1));
    //forlake.index_fixed.print();

    while(n_trials<25) {
        n_trials++;
        Timer.reset_sumTime();
        Timer.start();
        
        forlake.truncate();
        forlake.rec_grow(&X,&y);
        

        forlake.predict_all(&X_test,&y_pred_arr,0);
        print_error(y_test, y_pred);
        Timer.stop();    
        
        sum_time += Timer.elapsedTime;    
        sprint("  time is:"); sprint(float(Timer.elapsedTime/1000));
        sprint(" avg.:"); sprint(sum_time/n_trials/1000);
        sprint(" i_node:"); sprint(forlake.i_node);
        sprint(" i_tree:"); sprint(forlake.i_tree);
        sprint(" sort:"); sprintln(float(forlake.sortTimer.sumTime/1000));
        
    }
}
