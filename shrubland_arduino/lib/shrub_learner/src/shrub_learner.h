#pragma once
#include <dynarray.h>
#include <algorithm>

//#include <iomanip



struct node {
    node();
    float splitval;
    uint16_t bestvar, right_child_id;
    void makeTerminal(float, uint16_t);
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
    
    bool ownership;
    void bestsplit(temp_node*);
    dynamic_array<float,uint16_t>* X;
    dynamic_vector<float,uint16_t>* y;
    dynamic_array<uint16_t,uint16_t>* index;  
    
    public:
    
    node* nodes;
    uint16_t n_nodes;
    uint16_t i_node;
    uint16_t n_nodes_i_tree;
    uint16_t i_tree;


    forest_lake(node*,uint16_t);
    forest_lake(uint16_t);
    ~forest_lake();
    node& new_node();
    node* new_node2();
    node& new_tree();
    void new_tree2();
    bool more_trees();
    bool more_nodes();
    bool two_more_nodes();
    void truncate();
    void print_nlines(uint16_t);
    void print_splits(uint16_t);
    
    void grow(dynamic_array<float,uint16_t>* X, dynamic_vector<float,uint16_t>* y);
    float predict(dynamic_array<float,uint16_t>* X,uint16_t i_row);
};

