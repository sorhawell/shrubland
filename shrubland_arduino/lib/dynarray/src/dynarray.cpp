#include <Arduino.h>
#include "dynarray.h"

//simple error handle
void error(const char* msg) {
    Serial.println(msg); //print msg to Serial
    while(true) {
        delay(5000); //halt execution forever
    }
}

//handle different conversions
/* double stox(double* pd, std::string str) {return(double(str));}
  float stox(float* pd, std::string str) {return(float(str));}
      int stox(int* pd, std::string str) {return(int(str));}
    bool stox(bool* pd, std::string str) {return(bool(int(str)));}
uint16_t stox(uint16_t* pd, std::string str) {return(uint16_t(str));}
 */

//void error(const char* e) {throw runtime_error(e);}

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
   Serial.print("\n size:");
   Serial.println(sz);
   for(S i =0; i<5&&i<sz; i++) {
        Serial.println(i);
    }
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
    //Serial.println("\n dynarray deleted");
    if(ownership) {
        //Serial.println("with ownership to heap");
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
    Serial.print("\n dim: [");
    Serial.print(n_row);
    Serial.print("]x[");
    Serial.print(n_col);
    Serial.print("] size:");
    Serial.println(sz);
    
    for(int i=0 ; i<5; i++) {
        for(int j=0; j<n_col; j++) {
            Serial.print(this->at(i,j),5);
            Serial.print("  \t, ");
        }
        Serial.print("\n");
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

template class dynamic_array<double, int>;
template class dynamic_vector<double, int>;
template class dynamic_array<float, uint16_t>;
template class dynamic_vector<float, uint16_t>;
template class dynamic_array<uint16_t, uint16_t>;
template class dynamic_vector<uint16_t, uint16_t>;
template class dynamic_array <bool,uint16_t>;
template class dynamic_vector <bool,uint16_t>;

