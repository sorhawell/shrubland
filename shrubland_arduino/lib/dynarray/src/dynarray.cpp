


#ifdef testpc
    #include "dynarray.h"
    #include <algorithm>
    
    void sprint    (int x)                 {std::cout<<x;}
    void sprint    (uint16_t x)            {std::cout<<x;}
    void sprint    (const char x[])        {std::cout<<x;}
    void sprint    (float x)               {std::cout<<x;}
    void sprint    (float x, uint16_t n)   {std::cout<<x;}
    void sprintln  (int x)                 {std::cout<<x<<'\n';}
    void sprintln  (uint16_t x)            {std::cout<<x<<'\n';}
    void sprintln  (const char x[])        {std::cout<<x<<'\n';}
    void sprintln  (float x)               {std::cout<<x<<'\n';}

    void sprintln  (float x, uint16_t n)   {std::cout<<x<<'\n';}
    void sprintln()                        {std::cout   <<'\n';}

    void error(const char* e) {throw std::runtime_error(e);}

    
    //mimic arduino sq function

    double sq(double x) {return x * x;}
    float sq(float x) {return x * x;}

    



    

#else
    #include <Arduino.h>
    #include "dynarray.h"
    void sprint(int x) {Serial.print(x);}
    void sprint(uint16_t x) {Serial.print(x);}
    void sprint(const char x[]) {Serial.print(x);}
    void sprint(float x) {Serial.print(x);}
    void sprintln(int x) {Serial.println(x);}
    void sprintln(uint16_t x) {Serial.println(x);}
    void sprintln(const char x[]) {Serial.println(x);}
    void sprintln(float x) {Serial.println(x);}
    void sprint(float x, uint16_t n) {Serial.print(x,n);}
    void sprintln(float x, uint16_t n) {Serial.println(x,n);}
    void sprintln() {Serial.println();}
     
    //simple error handle
    void error(const char* msg) {
        sprint(msg); //print msg to Serial
        while(true) {
            delay(5000); //halt execution forever
        }
    }

/*     #define sprint(x) (Serial.print(x));
    #define sprintln(x) (Serial.println(x));
    #define sprint(x,n) (Serial.print(x,n));
    #define sprintln(x,n) (Serial.println(x,n)); */

#endif




//handle different conversions
/* double stox(double* pd, std::string str) {return(double(str));}
  float stox(float* pd, std::string str) {return(float(str));}
      int stox(int* pd, std::string str) {return(int(str));}
    bool stox(bool* pd, std::string str) {return(bool(int(str)));}
uint16_t stox(uint16_t* pd, std::string str) {return(uint16_t(str));}
 */

//void error(const char* e) {throw runtime_error(e);}
/* template <class T>
dynamic_segment<T>::dynamic_segment(T* begin, T* end):
    dyn_vec(begin), dyn_end(end) {};

template <class T>
T* dynamic_segment<T>::begin() {
    return(dyn_vec);
}
template <class T>
T* dynamic_segment<T>::end() {
    return(dyn_end);
} */

template <class T, class S>
dynamic_vector<T,S>::dynamic_vector(T* begin, S size):
    dyn_vec(begin), j(begin), sz(size) {};

template <class T, class S>
T dynamic_vector<T,S>::at(S i){
    if(i>=sz) error("bound exceeded");
    return(*(dyn_vec+i));

}

template <class T, class S>
T dynamic_vector<T,S>::get(S i){return(*(dyn_vec+i));}

template <class T, class S>
T* dynamic_vector<T,S>::get_cell_ref() {return(dyn_vec);}

template <class T, class S>
T dynamic_vector<T,S>::next() {
    if(j-dyn_vec>=sz) error("bound exceeded");
    return(*(j++));
}

template <class T, class S>
void dynamic_vector<T,S>::print() {
    constexpr S num_per_line{10};
    const S print_limit = (sz<50 ? sz:50);
    sprint("\n size: ");sprintln(sz);
    for(S i =0; i<print_limit; i++) {
        sprint(get(i));sprint("\t ");
        if((i+1)%num_per_line==0 && i+1!=print_limit) sprint("\n");
    }
    if(sz!=print_limit) sprint(". . .");;
}

template <class T, class S>
T* dynamic_vector<T,S>::begin() {
    return(dyn_vec);
}

template <class T, class S>
T* dynamic_vector<T,S>::end() {
    return(dyn_vec+sz);
}

template <class T, class S>
S dynamic_vector<T,S>::get_size() {
    return(sz);
}

///dynamic_array

template <class T, class S>
dynamic_array<T,S>::dynamic_array():
    dyn_array(nullptr), ownership(false), n_row(0), n_col(0),
    sz(0), j_col(0), j_row(0), j(0) {};

template <class T, class S>
dynamic_array<T,S>::dynamic_array(T* begin, S size, S ncol):
    dyn_array(begin), ownership(false), n_row(size/ncol), n_col(ncol),
    sz(size), j_col(0), j_row(0), j(begin) {};


template <class T, class S>
dynamic_array<T,S>::dynamic_array(S size, S ncol):
    dyn_array(nullptr), ownership(true), n_row(size/ncol), n_col(ncol),
    sz(size), j_col(0), j_row(0), j(nullptr) {
    dyn_array = new T[sz];
    j = dyn_array;
};

template <class T, class S>
dynamic_array<T,S>::~dynamic_array() {
    //sprint("\n dynarray deleted");
    if(ownership) {
        //sprint("\n delete on heap");sprintln(j_col);;sprintln(sz);
        delete[] dyn_array;
    }
    
};

template <class T, class S>
dynamic_array<T,S>::dynamic_array(const dynamic_array<T,S>& t)
    : dyn_array(t.dyn_array), ownership(t.ownership), n_row(t.n_row),
    n_col(t.n_col), sz(t.sz), j_col(t.j_col), j_row(t.j_row), j(t.j) {
    if(!t.ownership) error("dyn_arr without ownership cannot copy");
}

template <class T, class S>
dynamic_array<T,S> dynamic_array<T,S>::copy() {
    dynamic_array<T,S> new_da = *this;
    ownership = false; //revoke ownership in orginal
    return(new_da);
}

template <class T, class S>
dynamic_array<T,S> dynamic_array<T,S>::copy_weak() {
    dynamic_array<T,S> new_da = *this;
    new_da.ownership = false; //revoke ownership in copy
    return(new_da);
}

/* template <class T, class S>
dynamic_vector<T,S> dynamic_vector<T,S>::copy_weak() {
    dynamic_vector<T,S> new_da = *this;
    //new_da.ownership = false; //revoke ownership in copy
    return(new_da);
}
 */
/* template <class T, class S>
dynamic_array<T,S> dynamic_array<T,S>::clone() {
    dynamic_array new_da;
    memcpy(this,&new_da,sizeof(dynamic_array<T,S>));
    new_da.dyn_array = new T[sz];
    new_da.j = new_da.dyn_array+(j-dyn_array); //j point relatively same place in block
    new_da.ownership = true;

    return(new_da);
} */

template <class T, class S>
void dynamic_array<T,S>::assign(dynamic_array<T,S>& to) {
    to = *this;
    ownership = false; //rewoke ownership in original
}

template <class T, class S>
void dynamic_array<T,S>::assign_weak(dynamic_array<T,S>& to) {
    to = *this;
    to.ownership = false; //rewoke ownership in copy
}

template <class T, class S>
dynamic_array<T,S>& dynamic_array<T,S>::operator = (const dynamic_array<T,S>& t) {
    if(!t.ownership){
        error("dyn_arr without ownership cannot assign");
    }
    if(ownership) { //release old mem
        delete[] dyn_array;
        dyn_array=nullptr;
        j=nullptr;
    }
    memcpy(this,&t,sizeof(dynamic_array<T,S>));

}

template <class T, class S>
void dynamic_array<T,S>::release_ownership() {
    ownership=false;
};

template <class T, class S>
T dynamic_array<T,S>::at(S i_row,S i_col) {
    if(i_row>=n_row||i_col>=n_col) error("bounds exceeded");
    return(*(dyn_array+i_row+(i_col*n_row)));
}

template <class T, class S>
T dynamic_array<T,S>::get(S i_row,S i_col) {
    if(i_row>=n_row||i_col>=n_col) error("bounds exceeded");
    return(*(dyn_array+i_row+(i_col*n_row)));
};

template <class T, class S>
T dynamic_array<T,S>::next() {
    return(*(j++));
}

template <class T, class S>
T dynamic_array<T,S>::prev() {
    return(*(--j));
}

template <class T, class S>
T* dynamic_array<T,S>::get_col_p(S i_col) {
    return(dyn_array+(i_col*n_row));
}

template <class T, class S>
dynamic_vector<T,S> dynamic_array<T,S>::get_vector(S i_col) {
    dynamic_vector<T,S> dyn_vec(dyn_array+i_col*n_row,n_row);
    return(dyn_vec);
}

template <class T, class S>
void dynamic_array<T,S>::print() {

    sprint("\n dim: [");
    sprint(n_row);
    sprint("]x[");
    sprint(n_col);
    sprint("] size:");
    sprintln(sz);
    
    for(int i=0 ; i<5; i++) {
        for(int j=0; j<n_col; j++) {
            sprint(this->at(i,j),5);
            sprint("  \t, ");
        }
        sprint("\n");
    } 
}

template <class T, class S>
void dynamic_array<T,S>::truncate_col() {
    if(n_col!=0) {
        n_col--;
        sz -= n_row; 
        j_col=0;
    } else {
        error("cannot truncate empty array");
    }
}

template <class T, class S>
void dynamic_array<T,S>::push_row(T value) {
    if(j_col>=n_col || j_row>=n_row) error("columns exceeded");
    dyn_array[j_row+j_col++*n_row] = value;
}

template <class T, class S>
void dynamic_array<T,S>::next_row() {
        if(j_row>=n_row) error("rows exceeded");
        j_col = 0;
        j_row++;
}

template <class T, class S>
S dynamic_array<T,S>::get_size() {
    return(n_row*n_col);
}

template <class T, class S>
S dynamic_array<T,S>::get_col() {
    return(n_col);
}

template <class T, class S>
S dynamic_array<T,S>::get_row() {
    return(n_row);
}

template<class T>
T Abs(T x) {
    if(x >= 0) {
        return  x;
    } else {
        return -x;
    }
};
template float Abs(float);
template double Abs(double);


void print_error( dynamic_vector<float,uint16_t>& y_true, dynamic_vector<float,uint16_t>& y_pred) {
    
    uint16_t n =y_true.get_size();
    double mean{0};
    for(const auto& i : y_true) {mean += i;}
    mean /= n;
    

    double sqerr  = 0;
    double sqerr2 = 0;
    float pred=0;
    float yval=0;
    

    for(int i=0; i<n; i++) {
        pred = y_true.at(i);
        yval = y_pred.at(i);
        sqerr += sq(yval - pred);
        sqerr2+= sq(yval - mean);
    }

    float SD = sqrt(sqerr /(n-1));
    float SD2= sqrt(sqerr2/(n-1));
    sprint("model error is: "); sprint(SD );
    //sum_error += SD;
    //Serial.print(" avg.:"); Serial.print(sum_error/n,4);
    sprint("  total error is: "); sprint(SD2);

}


/* 
template <class T,class S>
void dynamic_array<T,S>::fillCSV(std::string path) {
    std::ifstream data(path);
    std::string line;
    j_row=0; j_col=0;
    while(std::getline(data,line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while(std::getline(lineStream,cell,',')) {
            T value = stox(&value,cell);
            push_row(value);
        }
        next_row();
    }
};
 */

//template class dynamic_array<double, int>;
//template class dynamic_vector<double, int>;
template class dynamic_array<float, uint16_t>;
template class dynamic_vector<float, uint16_t>;
template class dynamic_array<uint16_t, uint16_t>;
template class dynamic_vector<uint16_t, uint16_t>;
template class dynamic_array <bool,uint16_t>;
template class dynamic_vector <bool,uint16_t>;