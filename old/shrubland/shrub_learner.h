#pragma once
#include "dynarray.h"
#include <iostream>
//#include <iomanip



struct node {
    node();
    double splitval;
    int bestvar, right_child_id;
    void makeTerminal(double, uint16_t);
    void print();
    void print_body();
    void print_terminal();
    void print_header();
};

struct temp_node{    
    uint16_t n;
    double sum;
    uint16_t status;
    bool* innodep;
    node* nodep;

    double lchild_sum;
    double rchild_sum;
    uint16_t lchild_n;
    uint16_t rchild_n;

    temp_node();
    double mean_node();
};

class forest_lake {
    
    bool ownership;
    void bestsplit(temp_node*);
    public:
    
    node* nodes;
    int n_nodes;
    int i_node;
    int n_nodes_i_tree;
    int i_tree;
    dynamic_array<double,int>* X;
    dynamic_vector<double,int>* y;
    dynamic_array<uint16_t,uint16_t>* index;  

    forest_lake(node*,int);
    forest_lake(int);
    ~forest_lake();
    node& new_node();
    node* new_node2();
    node& new_tree();
    void new_tree2();
    bool more_trees();
    bool more_nodes();
    bool two_more_nodes();
    void truncate();
    void print_nlines(int);
    
    void grow(dynamic_array<double,int>* X, dynamic_vector<double,int>* y);
    double predict(dynamic_array<double,int>* X,int i_row);
};

