#pragma once
#include "dynarray.h"
#include <iostream>

template<int Tntree, int Tnnode>
struct Forest {
    int n_nodes;
    int n_trees;
    int i_node;
    int n_nodes_i_tree;
    int i_tree;
    node nodes[Tnnode];
    int tree_index[Tntree];
    
    Forest();
    node& new_node();
    node& new_tree();
    bool more_trees();
    bool more_nodes();
    bool two_more_nodes();
    void truncate();
};

//old_shrub.h
struct node_stack {
    node_stack();
    double crit_parent, critmax, parent_sum;
    int parent_n, splitcurr;
    bool* innode;
    void print();
};

template <class T, int Tnrow, int Tncol>
struct shrub_data {
  //stack varibles
    matrix<T,Tnrow,Tncol>& X;
    matrix<T,Tnrow,1>& y;

    matrix<uint16_t,Tnrow,Tncol> index;
    matrix<uint16_t,Tnrow,1> ybag;
    
    shrub_data(matrix<T,Tnrow,Tncol>& newX, matrix<T,Tnrow,1>& newy);
    void index_features();
    void set_X(matrix<T,Tnrow,Tncol>& X_param);
    void set_y(matrix<T,Tnrow,1>&     y_param);
};


//shrub a class implementing a variant of random forest
template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
class shrub {
private:
public:
  
    //pointers data, indexes
    shrub_data<T,Tnrow,Tncol>* sdata; //convenient struct of data and indexes
    matrix<T,Tnrow,Tncol>* X;
    matrix<T,Tnrow,1>* y;
    matrix<uint16_t,Tnrow,Tncol>* index;
    matrix<uint16_t,Tnrow,1>* ybag;

    //permanent allocatede data for the shrub learner
    Forest<Tntree,Tnnode> forest;
    int y_size, X_ncol, p_ntree, p_mtry, p_maxdepth, p_maxnodes, p_minnode;
    bool data_ready;
    

    shrub();
    void make_innode2(
        const node& this_node, const node_stack& this_stack_node,
        node_stack& stack_l_node, node_stack& stack_r_node);
    void bestsplit(node& this_node, node_stack& this_stack_node);
    void eval_node(node& this_node, node_stack& this_stack_node, int depth);
    void grow_forest();
    void link_to_data(shrub_data<T,Tnrow,Tncol>*);
    T predict(T* newX);
};