#pragma once

#ifdef testpc
    #include "../../dynarray/src/dynarray.h"
    #include "Profileapi.h"

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
#else
    #include <dynarray.h>
#endif

#include <algorithm>

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

    dynamic_array<float,uint16_t>* X = nullptr;
    dynamic_vector<float,uint16_t>* y = nullptr;
    dynamic_array<uint16_t,uint16_t>* index = nullptr;  
    
    bool ownership;
    
    public:
    dynamic_array<uint16_t,uint16_t> index_fixed;
    node* nodes = nullptr;
    uint16_t n_nodes;
    uint16_t i_node = 0;
    uint16_t n_nodes_i_tree = 0;
    uint16_t i_tree = 0;
    float temp_node_sum_left = 0.0;
    float temp_node_sum_right = 0.0;

    uint16_t p_depth = 15;
    uint16_t p_minnode = 5;
    uint16_t p_ntree = 500;
    uint16_t p_sampsize = 400;
    uint16_t p_minimal_index = 45;

    #ifdef testpc
        CPUtimer sortTimer;
    #endif
    forest_lake();
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
    bool pre_index_features(dynamic_array<float,uint16_t>* newX);

    void bestsplit(temp_node*);
    void grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy);
    void rec_grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy);
    void grow_node(uint16_t* Sp, uint16_t* Ep, node* parent_node, uint16_t depth);
    bool recsplit(uint16_t* Sp, uint16_t* Ep, uint16_t*& Cp, node* parent_node);
    float predict(dynamic_array<float,uint16_t>* X,uint16_t i_row);
    void predict_all(dynamic_array<float,uint16_t>* X,dynamic_array<float,uint16_t>* out,uint16_t i_col);
};



uint16_t* RadixSort(uint16_t * a, size_t count);
void my_shuff(uint16_t* begin, uint16_t* end);