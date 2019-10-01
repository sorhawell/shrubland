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

forest_lake forlake(12000);
uint16_t n_train=500;
float* Dp_train = DData;
float* Dp_test = DData_test;

dynamic_array<float,uint16_t> X(Dp_train    , n_train*ncol,ncol);
dynamic_array<float,uint16_t> X_test(Dp_test, n_train*ncol,ncol);
dynamic_vector<float,uint16_t> y = X.get_vector(ycol); //
dynamic_vector<float,uint16_t> y_test = X_test.get_vector(ycol); //
dynamic_array<float,uint16_t> y_pred_arr(n_train,1);
dynamic_vector<float,uint16_t> y_pred = y_pred_arr.get_vector(0);
    

//#define sprintf(b) (Serial.printf(b));

void setup() {


    Serial.begin(115200);
    Serial.println("ln29");
    X.print();
    
    X.truncate_col(); //revoke X's read/write access this last column
    X_test.truncate_col(); //revoke X's read/write access this last column
    
    
    sprint(42);
    sprint("hehej");

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
    forlake.rec_grow(&X,&y);

    //forlake.print_nlines(130);

    
    forlake.predict_all(&X_test,&y_pred_arr,0);
    print_error(y_test, y_pred);
    endTime = millis();

    //forlake.print_nlines(15);    

    Serial.print("  time is:"); Serial.print(endTime-startTime);
    sum_time += endTime-startTime;
    Serial.print(" avg.:"); Serial.print(sum_time/n_trials);
    Serial.print(" i_node:"); Serial.print(forlake.i_node);
    Serial.print(" i_tree:"); Serial.print(forlake.i_tree);
    
    //Serial.println(forlake.i_node);

    Serial.println("...");
}

