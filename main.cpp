//#include <algorithm>
//#include <string>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <random>
//#include <chrono>
//#include <iomanip>
#include "Profileapi.h"
#include "./shrubland_arduino/lib/dynarray/src/dynarray.h"
#include "./shrubland_arduino/lib/shrub_learner/src/shrub_learner.h"
#include "./shrubland_arduino/lib/testdata/src/testdata.h"
//#include "shrub_learner.h"



using namespace std;

struct CPUtimer {    
    CPUtimer():
    elapsedTime(0.0), sumTime(0.0) {}

    double elapsedTime;
    double sumTime;
    LARGE_INTEGER frequency;        // ticks per second
    LARGE_INTEGER t1, t2;           // ticks

    void setFreq();
    void start();
    void reset_sumTime();
    double stop();
    
};

void   CPUtimer::setFreq() {QueryPerformanceFrequency(&frequency);}
void   CPUtimer::start()   {QueryPerformanceCounter(&t1);}
void   CPUtimer::reset_sumTime() {sumTime=0;}
double CPUtimer::stop()  {
    QueryPerformanceCounter(&t2);
    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000000.0 / frequency.QuadPart;
    sumTime += elapsedTime;
    return(elapsedTime);
}    
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
    
int sum_time{0};
float sum_error{0};
int n_trials{0};
constexpr int nval = nrow*ncol;
int main() {

    
    sprintln(42);
    X.print();
    X.truncate_col(); //revoke X's read/write access this last column
    X_test.truncate_col(); //revoke X's read/write access this last column

    
     //std::shared_ptr<dynamic_array<double,int>> sp3(new dynamic_array<double,int>(5,10), array_deleter<dynamic_array<double,int>>());
    //std::cout << "\n ncopy:" <<sp3.use_count();
    Timer.setFreq();
    int startTime;
    int endTime;
   
    while(n_trials==0) {
        n_trials++;
        Timer.reset_sumTime();
        Timer.start();
        forlake.truncate();
        forlake.rec_grow(&X,&y);


        //forlake.print_nlines(80);

        

        double mean =0;
        for(int i=0; i<n_train;i++) {mean += y.at(i);}
        mean /= n_train;
        
        double sqerr  = 0;
        double sqerr2 = 0;
        float pred=0;
        float yval=0;
        for(int i=0; i<n_train; i++) {
            pred = forlake.predict(&X_test,i);
            yval = y_test.at(i);
            sqerr += sq(yval - pred);
            sqerr2+= sq(yval - mean);

/*             sprint(yval);sprint(" ,");sprint(pred);
            sprint(" ,");sprint  (X_test.at(i,0));
            sprint(" ,");sprintln(X_test.at(i,1)); 
 */            
        }
        Timer.stop();        
        float SD = sqrt(sqerr /(n_train-1));
        float SD2= sqrt(sqerr2/(n_train-1));
        sprint("model error is: "); sprint(SD );
        sum_error += SD;
        sprint(" avg.:"); sprint(float(sum_error/n_trials),4);
        sprint("  total error is: "); sprint(SD2);

        cout << "\n time used is: " << Timer.sumTime/1000;
    }
}

    //allocate data buffer on stack at run-time (heap and static is possible too)
    
    //double data_buffer[nval];
    //placing buffer in same scope as any abstrations using it best practice
    //as they will be desctructed simultanously
    
    
    /* std::cout << std::fixed;
    std::cout << std::setprecision(3);

   
    

    std::cout<<"\n i_nodes:" << forlake.i_node;
    
    double result;
    for(int i = 0;i<101;i++) {
        result = forlake.predict(&X2,i);
        std::cout<< "\n" << result ;//<< " , " << y.get(i);

    }
     */

    //forlake.print_nlines(50);

    //Forest<50,500> my_forest;




/*
    //safe access
    for(int i=0; i<5; i++) {cout << " " << X.at(i,0);}

    Timer.reset_sumTime();
    Timer.start();
    double temp;
    for(int i=0; i<nrow; i++) { temp =+ X.at(i,4);}
    Timer.stop();
    cout << "\n safe time: " << Timer.sumTime;

    //get raw pointer at own risk
    cout << "\n x4 ";
    double* x4 = X.get_col_p(4);
    for(int i=0; i<5; i++) {cout << " " << x4[i];}

    Timer.reset_sumTime();
    Timer.start();
    double temp2;
    for(int i=0; i<nrow; i++) { temp2 =+ x4[i];}
    Timer.stop();
    cout << "\n x4 time: " << Timer.sumTime;

    
    //get raw pointer at own risk
    cout << "\n v4 ";
    dynamic_vector<double,int> v4 = X.get_vector(4);
    for(int i=0; i<5; i++) {cout << " " << v4.at(i);}


    Timer.reset_sumTime();
    Timer.start();
    double temp3;
    for(int i=0; i<nrow; i++) { temp3 =+ v4.at(i);}
    Timer.stop();
    cout << "\n safe v4.at(i) time: " << Timer.sumTime;

    Timer.reset_sumTime();
    Timer.start();
    double temp4;
    for(int i=0; i<nrow; i++) { temp4 =+ v4.next();}
    Timer.stop();
    cout << "\n v4 next time: " << Timer.sumTime;

 */
    /* int a;
    int b;

    cin >>a >>b;

    int arr[a][b];

    arr[a-1][b-1]=0;

    cout <<"l604\n";    
    shrub<double,ntree,nnode,nrow,ncol> sf; //initialize a forest with 10 trees and 500 nodes(max)
    cout <<"l606\n";    

    Timer.setFreq();
    matrix<double,nrow,ncol> X;
    matrix<double,nrow,1   > y;
    X.fillCSV("./data/train1.csv");
    //X.print();
    
    for(int i = 0 ; i<X.nrow; i++) {y.a[0][i] = X.a[ycol][i];}
    X.erase_col(ycol,ycol+1);
    
    shrub_data<double,nrow,ncol> sd(X,y);
    sd.index_features();
    sf.link_to_data(&sd);
    
    sf.p_maxnodes = 250;
    sf.p_maxdepth = 12;
    sf.p_minnode = 5;
    cout << "hello";
    Timer.start();
    sf.grow_forest();
    Timer.stop();


    

    cout << "still alive time avg=" << Timer.sumTime << "\n";

    for(int i=0; i<sf.forest.n_trees;i++) cout << sf.forest.tree_index[i] << ", ";
    
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    //for(int i=0; i<1000;i++) {sf.forest.nodes[i].print_body();}
     
    //for(int i=0; i<150;i++) cout << sf.forest.tree_index[i] << ", ";
    sf.forest.nodes[0].print_header();
    int j_tree=0;
    //for(int i=0; i<sf.forest.i_node;i++) {
    for(int i=0; i<250;i++) {
        if(i==sf.forest.tree_index[j_tree]) {
            j_tree++;
            cout << "\n";
        }
        sf.forest.nodes[i].print_body();
        cout << "\t\t" << j_tree;
    }   */

