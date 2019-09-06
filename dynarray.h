#pragma once
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

double stox(double* pd, std::string str);
float  stox(float* pd, std::string str);
int    stox(int* pd, std::string str);

void error(const char* e);

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

    void print();
};

template <class T, class S>
class dynamic_array{
    
    T* dyn_array;
    bool ownership;
    T* j;
    S n_row;
    S n_col;
    S sz;
    S j_col;
    S j_row;
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
    void fillCSV(std::string);
    //void fillCSV_float(string path); //could be reworked to template specialization
};//dynarray



