#include <shrub_learner.h>

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
    Serial.println("value \t right_child_id \t bestvar");
    Serial.print("nsplitval     =");
    Serial.println(splitval);
    Serial.print("\nright_child_id=");
    Serial.println(right_child_id);
    Serial.print("\nbestvar       =");
    Serial.println(bestvar);
}

void node::print_header()  {Serial.println("\n node \t splitval \tright_child_id \tbestvar\t pred \t size \t tree");}
void node::print_body()    {
    Serial.print(splitval);Serial.print("\t\t");
    Serial.print(splitval);Serial.print("\t\t");
    Serial.print(right_child_id);Serial.print("\t\t");
    Serial.print(bestvar);Serial.print("\t\t");
}
void node::print_terminal(){
    Serial.print("\t\t\t\t\t");Serial.print(splitval);
    Serial.print("\t"); Serial.print(bestvar);Serial.print("\t");
}

void forest_lake::print_nlines(uint16_t nlines) {
//for(uint16_t i=0; i<150;i++) cout << sf.forest.tree_index[i] << ", ";
    
    uint16_t j_tree=0;
    //for(uint16_t i=0; i<sf.forest.i_node;i++) {
    for(uint16_t i=0; i<nlines && i<n_nodes;i++) {
        if(nodes[i].right_child_id==0 && nodes[i].bestvar==0 ) {
            j_tree++;
            Serial.println();
            nodes[0].print_header();
        }
        Serial.println(); Serial.print(i); Serial.print(" \t");
        if(nodes[i].right_child_id!=0) nodes[i].print_body();
        if(nodes[i].right_child_id==0 && nodes[i].bestvar!=0) nodes[i].print_terminal();
        Serial.print(j_tree);
    } 
}

void node::makeTerminal(float new_splitval, uint16_t terminal_n) {
    right_child_id = 0;
    bestvar = terminal_n;
    splitval = new_splitval;
}

temp_node::temp_node() :
    n(0),sum(0),status(0), innodep(nullptr), nodep(nullptr),
    lchild_sum(0), rchild_sum(0), lchild_n(0), rchild_n(0)  {}

float temp_node::mean_node() {
    return(sum/n);
}


//the forest pool, is memory pool for trees and nodes.
//A learner can acuire new nodes and trees from here.
//when pool is exhausted, the learner will gracefully complete its job.
forest_lake::forest_lake(node* nds,uint16_t size) :
    ownership(false), nodes(nds), n_nodes(size), i_node(0), n_nodes_i_tree(0), i_tree(0),
    X(nullptr), y(nullptr), index(nullptr) {};

forest_lake::forest_lake(uint16_t size) :
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
    i_tree=0;
    i_node=0;
    n_nodes_i_tree=0;
}

void forest_lake::bestsplit(temp_node* snode) {

    //internal variables
    uint16_t  tieVal{0};
    float crit{-1}, critmax{-1}, curr_y{0};
    
    float sum_l{0}, sum_r{0};   //rolling sum/substraction of target values
    uint16_t n_l{0}, n_r{0}; //count size of left and tight partitions
    
    uint16_t prev{0}, curr{0}, i_idx{0}; //dummy init
    const uint16_t n_obs = X->get_row();
//    uint32_t randVal{0};
//    bool none_found_yet{true};
    bool* innodep = snode->innodep;
   
    for(uint16_t i_var=0; i_var < X->get_col();i_var++) { //iterate each variable
        
        const auto* x = X->get_col_p(i_var); //reference to column for this variable in matrix
        const auto* i_var_indexed = index->get_col_p(i_var);
        
        //find fist obs which is innodes
        i_idx = 0;
        while(i_idx<n_obs) {
            curr = i_var_indexed[i_idx];
            if(innodep[curr]) break;
            i_idx++;
        }

        //place obs in right, rest in left
        sum_l = snode->sum - y->get(curr);
        sum_r = y->get(curr);
        n_l = snode->n - 1 ;
        n_r = 1;
        prev=curr;

        //check all splits
        for(i_idx=i_idx+1; i_idx<n_obs && n_l!=0;i_idx++) {
            
        
            //use index to iterate any split point sorted from low variable value to high
            curr = i_var_indexed[i_idx];

            //skip if not innode or if same value as last
            if(!innodep[curr]) {  //check if indexed obs is innode
                continue;         //look for next
            }

            //if(x[curr]!=x[prev]) {  //only calc crit if feature values are different
                
                //Serial.print(n_l);Serial.print(" ");
                
                //compute crit
                crit = (sum_l * sum_l / n_l) + (sum_r * sum_r / n_r);// - crit_parent; 
                
                //penalty for including less than 3 obs in innode
            /*  if(n_l<=3 || n_r<=3) {//check if same value
                    crit /= 12.7;
                } */

                //update crit handling, //swap crit== and crit> (the latter is more common)
                if (crit >= critmax && x[curr]!=x[prev]) {  //handle better crit, accept new split
                    if (crit == critmax) { //handle as tie
                        tieVal++;
                    } else { //handle as tie
                        tieVal = 1;
                        critmax = crit;
                    }
                    if (crit != critmax || uint_dist(rng) < UINT32_MAX / tieVal) { //update
                            snode->nodep->splitval = (x[curr] + x[prev]) / 2.0;
                            //Serial.print(crit);Serial.print(" ");Serial.print(n_l);Serial.print(" ");
                            //Serial.print(i_var);Serial.print(" ");
                            //Serial.println(snode->nodep->splitval );
                            snode->nodep->bestvar = i_var;
                            snode->lchild_n = n_l;
                            snode->rchild_n = n_r;
                            snode->lchild_sum = sum_l;
                            snode->rchild_sum = sum_r;
                    }
                }
            //}
            //update indexes etc.
            prev = curr;
            curr_y = y->get(curr);
            sum_r += curr_y;
            sum_l -= curr_y;
            n_r++; n_l--;
            
        }
    }

}


void forest_lake::grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy) {
    
    X = newX;
    y = newy;

    //short-hands
    bool* in_parent = nullptr;
    bool* in_child  = nullptr;
    float* x = nullptr;
    uint16_t* idx = nullptr;
    float splitval = 0;
                    

    //run time pars
    const uint16_t p_depth = 8;
    const uint16_t p_ntree = 50;
    const uint16_t p_sampsize = 330;
    
    
    //initialize allocate temporay data
    //make index of X, same dim as X
    //uint16_t index_buffer[X->get_size()];
    //dynamic_array<uint16_t,uint16_t> index_(index_buffer,X->get_size(),X->get_col());
    dynamic_array<uint16_t,uint16_t> index_(X->get_size(),X->get_col());
    index = &index_;
    for(uint16_t i=0 ; i<X->get_col(); i++) {
        dynamic_vector<uint16_t,uint16_t> idx_vec = index_.get_vector(i);
        float* x_vec = X->get_col_p(i);
        std::iota(idx_vec.begin(),idx_vec.end(),0);
        std::sort(idx_vec.begin(),idx_vec.end(),
            [x_vec](size_t i1, size_t i2) {return x_vec[i1] < x_vec [i2];}
        );
    }

    //make inbag/innode sampling
    uint16_t innode_sz{X->get_row()*p_depth};
    //bool innode_buffer[innode_sz];
    dynamic_array <bool,uint16_t> innode(innode_sz,p_depth);
    dynamic_vector<bool,uint16_t> innode_v0 = innode.get_vector(0);
    bool* innode_p0 = innode.get_col_p(0);
    uint32_t samp_limit = 0;
    for(uint16_t i=0; i<X->get_row(); i++) {innode_p0[i]= samp_limit++ < p_sampsize;}


    //allocate configure array tempoary nodes on stack (no more than p_depth)
    temp_node tnode[p_depth];
    for(uint16_t i = 0; i<p_depth;i++) tnode[i].innodep = innode.get_col_p(i);
    int16_t i_depth=0;     
    int16_t x____test_n=0;

    while(two_more_nodes() && p_ntree>i_tree-1) {
        std::shuffle(innode_v0.begin(),innode_v0.end(),rng);

/*         for(int i=0; i<5;i++) {
            Serial.print(innode_v0.at(i));
               
        } 
        Serial.println("\n"); */
        
        
        //configure first tempoary node
        new_tree2();//tag marks a new tree starts in forest lake
        tnode[0].nodep = new_node2();
        tnode[0].status = 0;
        tnode[0].n = 0;
        tnode[0].sum = 0;
        for(uint16_t i=0; i<X->get_row();i++) {
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
                if(i_depth>=p_depth-1 || tnode[i_depth].n<5 || !two_more_nodes() ) {

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
                    idx = index->get_col_p(tnode[i_depth].nodep->bestvar);
                    splitval = tnode[i_depth].nodep->splitval;
                    x____test_n=0;
                    //int x_test_innode=0;
                    for(uint16_t i=0; i<X->get_row(); i++) {
                        
                        in_child[i] = in_parent[i] && (x[i] >= splitval);
                        if(in_child[i]) x____test_n++;
                        //if(in_parent[i]) x_test_innode++;
                    }
                    if(x____test_n != tnode[i_depth].lchild_n) {
                      //  Serial.println(x_test_innode);
                        Serial.println(tnode[i_depth].nodep->bestvar);
                        Serial.print(" lsought:");Serial.print(tnode[i_depth].lchild_n);
                        Serial.print(" rsought:");Serial.print(tnode[i_depth].rchild_n);
                        Serial.print(" lfound: ");Serial.println(x____test_n);
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
                    idx = index->get_col_p(tnode[i_depth].nodep->bestvar);
                    splitval = tnode[i_depth].nodep->splitval;
                    x____test_n=0; //debug only
                    for(uint16_t i=0; i<X->get_row(); i++) {
                        in_child[i] = in_parent[i] && x[i] < splitval; //random tie handling also here?
                        if(in_child[i]) x____test_n++; //debug only
                    }
                    
                    //debug
                    if(x____test_n != tnode[i_depth].rchild_n) {
                        Serial.print(" lsought:");Serial.print(tnode[i_depth].lchild_n);
                        Serial.print(" rsought:");Serial.print(tnode[i_depth].rchild_n);
                        Serial.print(" rfound: ");Serial.println(x____test_n);
                        error("wrong right n");
                    }

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


float forest_lake::predict(dynamic_array<float,uint16_t>* X,uint16_t i_row){
    
    //copy 
    float xnew[X->get_col()];
    for(uint16_t i_col=0;i_col<X->get_col();i_col++) xnew[i_col] = X->get(i_row,i_col);

    uint16_t j_tree=0;   //count tree predictions found
    float prediction=0;  //sum predictions
    uint16_t k_node = 0; //iterate nodes in one tree

    for(uint16_t j_node=0; j_node < n_nodes; j_node++) { //iterate all nodes in forest lake to find trees
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