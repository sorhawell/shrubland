#pragma once

#ifdef testpc
    #include "../../dynarray/src/dynarray.h"
#else
    #include <dynarray.h>
#endif

#include <algorithm>

//#include <iomanip



struct node {
    node();
    float splitval;
    uint16_t bestvar, right_child_id;
    void makeTerminal(float, uint16_t);
    void makeTerminal(uint16_t);
    void print();
    void print_body();
    void print_terminal();
    void print_header();
};

struct temp_node{    
    uint16_t n;
    float sum;
    uint16_t status;
    bool* innodep;
    node* nodep;

    float lchild_sum;
    float rchild_sum;
    uint16_t lchild_n;
    uint16_t rchild_n;

    temp_node();
    float mean_node();
};

class forest_lake {
    
   

    dynamic_array<float,uint16_t>* X;
    dynamic_vector<float,uint16_t>* y;
    dynamic_array<uint16_t,uint16_t>* index;  
    bool ownership;
    
    public:
    node* nodes;
    uint16_t n_nodes;
    uint16_t i_node;
    uint16_t n_nodes_i_tree;
    uint16_t i_tree;
    float temp_node_sum_left;
    float temp_node_sum_right;


    forest_lake(node*,uint16_t);
    forest_lake(uint16_t);
    ~forest_lake();
    node& new_node();
    node* new_node2();
    node* new_node2(float);
    node* new_node3();
    node* new_node3(float);
    node& new_tree();
    void new_tree2();
    bool more_trees();
    bool more_nodes();
    bool two_more_nodes();
    void truncate();
    void print_nlines(uint16_t);
    void print_splits(uint16_t);
    
    void bestsplit(temp_node*);
    void grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy);
    void rec_grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy);
    void grow_node(uint16_t* Sp, uint16_t* Ep, node* parent_node, uint16_t depth);
    bool recsplit(uint16_t* Sp, uint16_t* Ep, uint16_t*& Cp, node* parent_node);
    float predict(dynamic_array<float,uint16_t>* X,uint16_t i_row);
    void predict_all(dynamic_array<float,uint16_t>* X,dynamic_array<float,uint16_t>* out,uint16_t i_col);
};

