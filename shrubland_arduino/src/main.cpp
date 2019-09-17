#include <Arduino.h>
#include <algorithm>
#include <dynarray.h>
#include <shrub_learner.h>
#include <testdata.h>


using namespace std;

constexpr int ntree=50;
constexpr int nnode=50*300;
constexpr int nrow =1000;
constexpr int ncol =6;
constexpr int ycol =5;

forest_lake forlake(10000);
uint16_t n_train=500;
float* Dp_train = DData;
float* Dp_test = DData_test;

dynamic_array<float,uint16_t> X(Dp_train    , n_train*ncol,ncol);
dynamic_array<float,uint16_t> X_test(Dp_test, n_train*ncol,ncol);
dynamic_vector<float,uint16_t> y = X.get_vector(ycol); //
dynamic_vector<float,uint16_t> y_test = X_test.get_vector(ycol); //
    

void setup() {
    Serial.begin(115200);
    Serial.println("ln29");

    X.print();
    
    X.truncate_col(); //revoke X's read/write access this last column
    X.truncate_col();
    X.truncate_col();
    X_test.truncate_col(); //revoke X's read/write access this last column
    X_test.truncate_col(); //revoke X's read/write access this last column
    X_test.truncate_col(); //revoke X's read/write access this last column
    

}


int sum_time{0};
float sum_error{0};
int n_trials{0};
void loop(){
    n_trials++;
 
    //X2 = X; // not public accessible use instead .assign method instead
    int startTime;
    int endTime;
    startTime = millis();
    forlake.truncate();
    forlake.grow(&X,&y);
    

    double mean =0;
    for(int i=0; i<n_train;i++) {mean += y  .at(i);}
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

        //Serial.print(yval);Serial.print(" ,");Serial.println(pred);
    }
    
    double SD = sqrt(sqerr /(n_train-1));
    double SD2= sqrt(sqerr2/(n_train-1));
    Serial.print("model error is: "); Serial.print(SD );
    sum_error += SD;
    Serial.print(" avg.:"); Serial.print(sum_error/n_trials,4);
    Serial.print("  total error is: "); Serial.print(SD2);
    endTime = millis();

    //forlake.print_nlines(15);    

    Serial.print("  time is:"); Serial.print(endTime-startTime);
    sum_time += endTime-startTime;
    Serial.print(" avg.:"); Serial.print(sum_time/n_trials);
    Serial.print(" forlake.i_node"); Serial.print(forlake.i_node);
    
    //Serial.println(forlake.i_node);

    Serial.println("...");
}

