#include "shrub_learner.h"
#include <algorithm>


typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val;           // populate somehow
MyRNG rng;                   // e.g. keep one global instance (per thread)
void initialize() {
  rng.seed(seed_val);
}
std::uniform_int_distribution<uint32_t> uint_dist;

//stack_node is tempoary node information placed on the stack
//stack_node information is dropped when node goes out of scope of recursive algorithm

//node is permenant information of a node. Nodes are saved to a forest pool (arena)
node::node() : splitval(0), bestvar(0), right_child_id(0) {};
void node::print() {
    std::cout << "value \t right_child_id \t bestvar" <<
    "nsplitval     =" << splitval          <<
    "\nright_child_id=" << right_child_id    <<
    "\nbestvar       =" << bestvar           <<'\n';
}

void node::print_header()  {std::cout << "\n node \t splitval \tright_child_id \tbestvar\t pred \t size \t tree";}
void node::print_body()    {std::cout << splitval<<"\t\t"<< right_child_id<<"\t\t"<< bestvar << "\t\t\t";}
void node::print_terminal(){std::cout << "\t\t\t\t\t"<< splitval << "\t" << bestvar << "\t";}
void forest_lake::print_nlines(int nlines) {
//for(int i=0; i<150;i++) cout << sf.forest.tree_index[i] << ", ";
    
    int j_tree=0;
    //for(int i=0; i<sf.forest.i_node;i++) {
    for(int i=0; i<nlines && i<n_nodes;i++) {
        if(nodes[i].right_child_id==0 && nodes[i].bestvar==0 ) {
            j_tree++;
            std::cout << "\n";
            nodes[0].print_header();
        }
        std::cout << '\n' << i << " \t";
        if(nodes[i].right_child_id!=0) nodes[i].print_body();
        if(nodes[i].right_child_id==0 && nodes[i].bestvar!=0) nodes[i].print_terminal();
        std::cout << j_tree;
    } 
}

void node::makeTerminal(double new_splitval, uint16_t terminal_n) {
    right_child_id = 0;
    bestvar = terminal_n;
    splitval = new_splitval;
}

temp_node::temp_node() :
    n(0),sum(0),status(0), innodep(nullptr), nodep(nullptr),
    lchild_sum(0), rchild_sum(0), lchild_n(0), rchild_n(0)  {}

double temp_node::mean_node() {
    return(sum/n);
}


//the forest pool, is memory pool for trees and nodes.
//A learner can acuire new nodes and trees from here.
//when pool is exhausted, the learner will gracefully complete its job.
forest_lake::forest_lake(node* nds,int size) :
    ownership(false), nodes(nds), n_nodes(size), i_node(0), n_nodes_i_tree(0), i_tree(0),
    X(nullptr), y(nullptr), index(nullptr) {};

forest_lake::forest_lake(int size) :
ownership(true), nodes(nullptr), n_nodes(size), i_node(0), n_nodes_i_tree(0), i_tree(0),
X(nullptr), y(nullptr), index(nullptr) {
    nodes = new node[n_nodes];
};

forest_lake::~forest_lake() {
    if(ownership) {
        delete[] nodes;
    }
}

bool forest_lake::more_nodes() {
    return(i_node<n_nodes);
}

bool forest_lake::two_more_nodes() {
    return(i_node < n_nodes-1);
}

node& forest_lake::new_node() {
    if(i_node==n_nodes) error("no more nodes");
    n_nodes_i_tree++;
    node& ref = nodes[i_node++];
    return(ref);
}

node* forest_lake::new_node2() {
    if(i_node==n_nodes) error("no more nodes");
    n_nodes_i_tree++;
    return(&nodes[i_node++]);
}

node& forest_lake::new_tree() {
    if(!two_more_nodes()) error("no more nodes");
    
    node& tree_start = new_node();
    tree_start.right_child_id=0; //all normal nodes have non-zero child ids
    tree_start.bestvar=0;        //signals this is a start node
    n_nodes_i_tree=1;
    i_tree++;
    return(new_node());
}

void forest_lake::new_tree2() {
    if(!two_more_nodes()) error("no more nodes");
    node& tree_start = new_node();
    tree_start.right_child_id=0; //all normal nodes have non-zero child ids
    tree_start.bestvar=0;        //signals this is a start node
    n_nodes_i_tree=1;
    i_tree++;
}

void forest_lake::truncate() {
    i_tree, i_node, n_nodes_i_tree = 0;
}

void forest_lake::bestsplit(temp_node* snode) {

    //internal variables
    double sum_l{0}, sum_r{0};   //rolling sum/substraction of target values
    int n_l{0}, n_r{0}, tieVal{0}; //count size of left and tight partitions
    double crit{0}, critvar{0}, critmax{0};
    int prev{0}, curr{0}; //dummy init
    bool none_found_yet{true};
    bool* innodep = snode->innodep;
   
    for(int i_var=0; i_var < X->get_col();i_var++) { //iterate each variable
        
        const auto* x = X->get_col_p(i_var); //reference to column for this variable in matrix
        
        //reset variables for rolling mse loss function
        sum_l = snode->sum;
        sum_r = 0;
        n_l = snode->n;
        n_r = 0;
     
        auto* i_var_indexed = index->get_col_p(i_var);
        none_found_yet = true;

        for(auto i_idx=0; i_idx<index->get_row();i_idx++) {
            //use index to iterate any split point sorted from low variable value to high
            curr = i_var_indexed[i_idx];
            //std::cout << "\n v: " << i_var << " c: " << curr << " p: " << prev;
             
            
            if(!innodep[curr]) {   //check if indexed obs is innode
                //std::cout << "skip because as not innode \n";
                continue;
            }
            if(none_found_yet) {//check if first found innode
                none_found_yet = false;
                prev = curr;
                //std::cout << "skip because first \n";
                continue;
            }
            if(x[curr]==x[prev]) {//check if same value
                //std::cout << "skip because same variable value \n";
                continue;
            }
            
            //compute crit
            sum_l -= y->get(curr);
            sum_r += y->get(curr);
            crit = (sum_l * sum_l / --n_l) + (sum_r * sum_r / ++n_r);// - crit_parent;            
            //std::cout << "crit: " << crit << '\n';
            
            if(n_l<=2 || n_r<=2) {//check if same value
                crit /= 3.2;
                //std::cout << "skip because same variable value \n";
                //continue;
            }


            //update crit handling, //swap crit== and crit> (the latter is more common)
            if (crit == critmax) { //handle tie, random accept new split
                //std::cout << "\n tie handling - ";
                tieVal++;
                uint32_t randVal = uint_dist(rng);
                if (randVal < UINT32_MAX / tieVal) {
                    snode->nodep->splitval = (x[curr] + x[prev]) / 2.0;
                    snode->nodep->bestvar = i_var;
                    snode->lchild_n = n_l;
                    snode->rchild_n = n_r;
                    snode->lchild_sum = sum_l;
                    snode->rchild_sum = sum_r;
                }
            } else if (crit > critmax) {  //handle better crit, accept new split
                snode->nodep->splitval = (x[curr] + x[prev]) / 2.0;
                snode->nodep->bestvar = i_var;
                tieVal = 1;
                critmax = crit;
                snode->lchild_n = n_l;
                snode->rchild_n = n_r;
                snode->lchild_sum = sum_l;
                snode->rchild_sum = sum_r;
                //std::cout << "\n new critvar: " << critvar << '\n';
            }

            //update indexes
            prev = curr;
        }
    }

}


void forest_lake::grow(dynamic_array<double,int>* newX, dynamic_vector<double,int>* newy) {
    
    X = newX;
    y = newy;

    //short-hands
    bool* in_parent = nullptr;
    bool* in_child  = nullptr;
    double* x = nullptr;
    double splitval = 0;
                    

    //run time pars
    const int p_depth = 35;
    const int p_ntree = 500;
    const uint16_t p_sampsize = 700;
    
    
    //initialize allocate temporay data
    //make index of X, same dim as X
    //uint16_t index_buffer[X->get_size()];
    //dynamic_array<uint16_t,uint16_t> index_(index_buffer,X->get_size(),X->get_col());
    dynamic_array<uint16_t,uint16_t> index_(X->get_size(),X->get_col());
    index = &index_;
    for(int i=0 ; i<X->get_col(); i++) {
        dynamic_vector<uint16_t,uint16_t> idx_vec = index_.get_vector(i);
        double* x_vec = X->get_col_p(i);
        std::iota(idx_vec.begin(),idx_vec.end(),0);
        std::sort(idx_vec.begin(),idx_vec.end(),
            [x_vec](size_t i1, size_t i2) {return x_vec[i1] < x_vec [i2];}
        );
    }

    //make inbag/innode sampling
    int innode_sz{X->get_row()*p_depth};
    //bool innode_buffer[innode_sz];
    dynamic_array <bool,uint16_t> innode(innode_sz,p_depth);
    dynamic_vector<bool,uint16_t> innode_v0 = innode.get_vector(0);
    bool* innode_p0 = innode.get_col_p(0);
    uint32_t samp_limit = 0;
    for(int i=0; i<X->get_row(); i++) {innode_p0[i]= samp_limit++ < p_sampsize;}

    bool split_bool[X->get_row()];


    //allocate configure array tempoary nodes on stack (no more than p_depth)
    temp_node tnode[p_depth];
    for(int i = 0; i<p_depth;i++) tnode[i].innodep = innode.get_col_p(i);
    int16_t i_depth=0;     
    int16_t x____test_n=0;

    while(two_more_nodes() && p_ntree>i_tree-1) {
        std::shuffle(innode_v0.begin(),innode_v0.end(),rng);

        
        
        //configure first tempoary node
        new_tree2();//tag marks a new tree starts in forest lake
        tnode[0].nodep = new_node2();
        tnode[0].status = 0;
        tnode[0].n = 0;
        tnode[0].sum = 0;
        for(int i=0; i<X->get_row();i++) {
            if(tnode[0].innodep[i]) {
                tnode[0].sum += y->get(i);
                tnode[0].n++;
            }
        }
        
        //grow tree
        while(i_depth>0 || tnode[0].status<2) {
            
            tnode[i_depth].status++; //count visits to node
            
            //conditions to break loop
            //out of nodes/trees, too deep, too few nodes
            if(tnode[i_depth].status == 1) {
                
                //a new node
                if(i_depth>=p_depth-1 || tnode[i_depth].n<7 || !two_more_nodes() ) {

                    tnode[i_depth].nodep->makeTerminal(tnode[i_depth].mean_node(),tnode[i_depth].n);
                    i_depth--;
                    continue;
                } else {
                    //start splitting

                    bestsplit(&tnode[i_depth]); //find best split

                    //allocate right child node, no temp_node for now
                    tnode[i_depth].nodep->right_child_id = i_node; //place right child id in this parent node
                    new_node2(); // allocate right child now for later

                    //allocate left child node, and configure the temp node
                    tnode[i_depth+1].nodep = new_node2(); //allocate left child
                    tnode[i_depth+1].sum = tnode[i_depth].lchild_sum;
                    tnode[i_depth+1].n   = tnode[i_depth].lchild_n;

                                   
                    //update left child temp_node innode by bestsplit
                    in_parent = tnode[i_depth  ].innodep;
                    in_child  = tnode[i_depth+1].innodep;
                    x = X->get_col_p(tnode[i_depth].nodep->bestvar);
                    splitval = tnode[i_depth].nodep->splitval;
                    x____test_n=0;
                    for(int i=0; i<X->get_row(); i++) {
                        in_child[i] = in_parent[i] && x[i] >= splitval;
                        if(in_child[i]) x____test_n++;
                    }
                    if(x____test_n != tnode[i_depth].lchild_n) {
                        error("wrong left n");
                    }

                     //go to left child
                    tnode[i_depth+1].status=0; //reset child status(n_visits) to zero
                    i_depth++;
                    continue;
                }

            } else {
                if(tnode[i_depth].status == 2) {

                    //recover right node pointer with ugly pointer arithmetics
                    tnode[i_depth+1].nodep = nodes + tnode[i_depth].nodep->right_child_id;
                    tnode[i_depth+1].sum = tnode[i_depth].rchild_sum;
                    tnode[i_depth+1].n   = tnode[i_depth].rchild_n;

                    //update right child temp_node innode by best split
                    in_parent = tnode[i_depth  ].innodep;
                    in_child  = tnode[i_depth+1].innodep;
                    x = X->get_col_p(tnode[i_depth].nodep->bestvar);
                    splitval = tnode[i_depth].nodep->splitval;
                    x____test_n=0; //debug only
                    for(int i=0; i<X->get_row(); i++) {
                     
                        in_child[i] = in_parent[i] && x[i] < splitval; //random tie handling also here?
                        if(in_child[i]) x____test_n++; //debug only
                    }
                    
                    //debug
                    if(x____test_n != tnode[i_depth].rchild_n) error("wrong right n");

                    //go to right child
                    tnode[i_depth+1].status=0;
                    i_depth++;
                    continue;

                } else {
                    if(tnode[i_depth].status == 3) {
                        //comming back from right - go up
                        i_depth--;
                        continue;
                    }else {
                        error("how did this happen!!?!?");
                    }
                }
            }
        }
    }


    //clean-up and clear pointers
    X = nullptr;
    y = nullptr;
    index = nullptr;

}


double forest_lake::predict(dynamic_array<double,int>* X,int i_row){
    
    //copy 
    double xnew[X->get_col()];
    for(int i_col=0;i_col<X->get_col();i_col++) xnew[i_col] = X->get(i_row,i_col);

    int j_tree=0; //count tree predictions found
    double prediction=0; //sum predictions
    int k_node = 0;      //iterate nodes in one tree

    for(int j_node=0; j_node < n_nodes; j_node++) { //iterate all nodes in forest lake to find trees
        if(nodes[j_node].right_child_id==0 &&  nodes[j_node].bestvar==0) { //tree found
            if(j_tree>=i_tree-1) break;   
            k_node = j_node + 1; //iterate tree
            while(true) {

                if(nodes[k_node].right_child_id!=0) {
                    
                    //this node is intermediary
                    if(xnew[nodes[k_node].bestvar] >= nodes[k_node].splitval) {
                        k_node++;
                    } else {
                        k_node = nodes[k_node].right_child_id;
                    }

                } else {
                    //this node is terminal
                    j_tree++;
                    prediction += nodes[k_node].splitval;
                    j_node = k_node; //new tree cant possibly start before k_node
                    break;
                }
            }
        }
    }
    return(prediction/j_tree);
}