#pragma once
#ifdef testpc
    #include <iostream>
    double sq(double x);
    float sq(float x);
    
#else
    #include <Arduino.h>
#endif
void sprint(int x);
void sprint(uint16_t x);
void sprint(float x);
void sprint(float x, uint16_t n);
void sprint(const char x[]);
void sprintln(int x);
void sprintln(uint16_t x);
void sprintln(float x);
void sprintln(float x,uint16_t n);
void sprintln(const char x[]);
void sprintln();


//#include <adaptor.h>

//print adaptors

//void error(const char* e);
void error(const char* msg);

template <class T, class S>
class dynamic_vector {
    T* dyn_vec;
    T* j;
    S sz;
    public:
    
    dynamic_vector(T* begin, S size);
    T at(S i);
    T get(S i);
    T* get_cell_ref();
    T next();
    T* begin();
    T* end();
    S get_size();

    void print();
};


template <class T, class S>
class dynamic_array{
    
    T* dyn_array;
    bool ownership;
    S n_row;
    S n_col;
    S sz;
    S j_col;
    S j_row;
    T* j;
    dynamic_array(const dynamic_array<T,S>&);
    dynamic_array<T,S>& operator = (const dynamic_array<T,S> &other);

    public:
    dynamic_array(T* begin, S size, S ncol);
    dynamic_array( S size, S ncol);
    ~dynamic_array();
    dynamic_array<T,S> copy();
    dynamic_array<T,S> copy_weak();
    void assign(dynamic_array<T,S>&);
    void assign_weak(dynamic_array<T,S>&);
    //dynamic_array<T,S> dynamic_array<T,S> clone();

    T at(S i_row,S i_col);
    T get(S i_row,S i_col);
    void push_row(T);
    void next_row();
    T next();
    T prev();
    T* get_col_p(S i_col);
    dynamic_vector<T,S> get_vector(S i_col);
    //dynamic_vector<T,S> get_row_vector(S i_row);
    void print();
    void truncate_col();
    S get_size();
    S get_col();
    S get_row();
    //void fillCSV(string path);
    //void fillCSV(std::string);
    //void fillCSV_float(string path); //could be reworked to template specialization
};//dynarray


template<class T>
T Abs(T);

void print_error( dynamic_vector<float,uint16_t>& y_true, dynamic_vector<float,uint16_t>& y_pred);