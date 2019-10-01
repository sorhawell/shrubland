#include "shrub_learner.h"

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
    sprintln("value \t right_child_id \t bestvar");
    sprint("nsplitval     =");
    sprintln(splitval);
    sprint("\nright_child_id=");
    sprintln(right_child_id);
    sprint("\nbestvar       =");
    sprintln(bestvar);
}

void node::print_header()  {sprintln("\n node \t splitval \tright_child_id \tbestvar\t pred \t size \t tree");}
void node::print_body()    {
    sprint(splitval);sprint("\t\t");
    sprint(right_child_id);sprint("\t\t");
    sprint(bestvar);sprint("\t\t\t");
}
void node::print_terminal(){
    sprint("\t\t\t\t\t");sprint(splitval);
    sprint("\t"); sprint(bestvar);sprint("\t");
}

void forest_lake::print_nlines(uint16_t nlines) {
//for(uint16_t i=0; i<150;i++) cout << sf.forest.tree_index[i] << ", ";
    
    uint16_t j_tree=0;
    //for(uint16_t i=0; i<sf.forest.i_node;i++) {
        sprintln();
        nodes[0].print_header();
    for(uint16_t i=0; i<nlines && i<n_nodes;i++) {

        if(nodes[i].right_child_id==0 && nodes[i].bestvar==0 ) {    
            j_tree++;
        }
        sprintln(); sprint(i); sprint(" \t");
        if(nodes[i].right_child_id!=0) nodes[i].print_body();
        if(nodes[i].right_child_id==0 && nodes[i].bestvar!=0) nodes[i].print_terminal();
        sprint(j_tree);
    } 
}

void forest_lake::print_splits(uint16_t nlines) {
    uint16_t j_tree=0;
    uint16_t depth=0;
    uint16_t max_depth=15;
    sprint(j_tree);
    for(uint16_t i=0; i<nlines && i<n_nodes;i++) {
        if(nodes[i].right_child_id==0 && nodes[i].bestvar==0 ) {    
            if(nodes[i+1].right_child_id==0 && nodes[i+1].bestvar==0) break;
            depth = 0;
            sprint("\n t");sprint(++j_tree);sprint(" ");
        }
        if(++depth<=max_depth) {sprint(nodes[i].splitval);sprint(',');}
    } 
};

void node::makeTerminal(float pred_value, uint16_t terminal_n) {
    right_child_id = 0;
    bestvar = terminal_n;
    splitval = pred_value;
}

void node::makeTerminal(uint16_t terminal_n) {
    right_child_id = 0;
    bestvar = terminal_n;
    splitval /= terminal_n;
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
    X(nullptr), y(nullptr), index(nullptr) ,
    ownership(false), nodes(nds), n_nodes(size), i_node(0), n_nodes_i_tree(0), i_tree(0),
    temp_node_sum_left(0),temp_node_sum_right(0) {};

forest_lake::forest_lake(uint16_t size) :
    X(nullptr), y(nullptr), index(nullptr) ,
    ownership(true), nodes(nullptr), n_nodes(size), i_node(0), n_nodes_i_tree(0), i_tree(0),
    temp_node_sum_left(0),temp_node_sum_right(0) {
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

node* forest_lake::new_node2(float sum) {
    if(i_node==n_nodes) error("no more nodes");
    n_nodes_i_tree++;
    nodes[i_node++].splitval = sum;
    return(&nodes[i_node]);
}

node* forest_lake::new_node3(float sum) {
    if(i_node==n_nodes) error("no more nodes");
    node* newnode = &nodes[i_node];
    n_nodes_i_tree++;
    nodes[i_node].splitval = sum;
    i_node++;
    return(newnode);
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

    constexpr uint16_t p_minsplit=2;
    if(p_minsplit>=snode->n) error("p_minsplit is set too high");

    //internal variables
    
    float crit{-1}, critmax{-1}, curr_y{0}, sum_r{0}, sum_l{0};
    uint16_t tieVal{0}, prev{0}, curr{0}, i_idx{0}, n_r{0}, n_l{0}; //dummy init
    const uint16_t n_obs = X->get_row();
    bool* innodep = snode->innodep;
    float* yp = y->begin();
    
    for(uint16_t i_var=0; i_var < X->get_col();i_var++) { //iterate each variable
        
        const auto* x = X->get_col_p(i_var); //reference to column for this variable in matrix
        const auto* i_var_indexed = index->get_col_p(i_var);
        i_idx = 0;
        
        //find second obs which is innodes, save information to curr and i_idx;
        sum_r = 0;
        sum_l = snode->sum;
        
        for(uint8_t i =0;i<p_minsplit;i++) {
            while(!innodep[(curr = i_var_indexed[i_idx++])]);
            curr_y = yp[curr];
            sum_r += curr_y;
            sum_l -= curr_y;
        }
        prev=curr;
        n_r = p_minsplit;
        n_l = snode->n - p_minsplit ;
        
        //check all splits
        for(; i_idx<n_obs; i_idx++) {
            curr = i_var_indexed[i_idx];
            if(!innodep[curr]) continue;  //check if indexed obs is innode
            
            //if(n_r>2 || crit==-1) {
                crit = (sum_l * sum_l / n_l) + (sum_r * sum_r / n_r);// - crit_parent;
                if (crit >= critmax && x[curr]!=x[prev]) {  //handle better crit, accept new split
                    if (crit!=critmax) tieVal=0;
                    if (crit!=critmax || uint_dist(rng) > UINT32_MAX / ++tieVal) { //update
                            snode->nodep->splitval = (x[curr] + x[prev]) / 2.0;
                            snode->nodep->bestvar = i_var;
                            snode->lchild_n = n_l;
                            snode->rchild_n = n_r;
                            snode->lchild_sum = sum_l;
                            snode->rchild_sum = sum_r;
                    }
                    critmax = crit;
                }
            //}
            //update indexes etc.
            prev = curr;
            curr_y = yp[curr];
            sum_r += curr_y;
            sum_l -= curr_y;
            n_r++;
            if(--n_l<=p_minsplit) break;
            
        }
    }

    if(snode->lchild_n == 0 ||snode->rchild_n == 0 || 
    snode->nodep->bestvar < 0 || snode->nodep->bestvar >= X->get_col()
    ) {
        sprintln("-");
        sprintln(snode->lchild_n);
        sprintln(snode->rchild_n);
        sprintln(snode->nodep->bestvar);
        sprintln(snode->sum);
        sprintln(snode->n);
        error("oups this split is invalid");
    }
}


void forest_lake::grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy) {
    
    X = newX;
    y = newy;

    //short-hands
    bool* in_parent = nullptr;
    bool* in_child  = nullptr;
    float* x = nullptr;
    //uint16_t* idx = nullptr;
    float splitval = 0;
                    

    //run time pars
    const uint16_t p_depth = 7;
    const uint16_t p_minnode = 5;
    const uint16_t p_ntree = 1;
    const uint16_t p_sampsize = 150;
    
    
    //initialize allocate temporay data
    //make index of X, same dim as X
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
    #ifdef testpc
        uint16_t innode_sz= uint16_t(p_depth*X->get_row()); //safe explicit cast as n_obs is limited to uint16_t
    #else
        uint16_t innode_sz= uint16_t(p_depth*X->get_row()); //safe explicit cast as n_obs is limited uint16_t max
    #endif

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
    
    // for debugging
    //int16_t x____test_n=0;
    //float   x____test_sum=0;
    //-----

    while(two_more_nodes() && i_tree<p_ntree) {
        std::shuffle(innode_v0.begin(),innode_v0.end(),rng);

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
            if(i_depth<0 || i_depth>p_depth) {
                error("seems we're in too deep or shallow");
            }

            tnode[i_depth].status++; //count visits to node
            
            //conditions to break loop
            //out of nodes/trees, too deep, too few nodes
            if(tnode[i_depth].status == 1) {
                
                //terminalize this node or split it?
                if(i_depth>=p_depth-1 || tnode[i_depth].n<p_minnode || !two_more_nodes() ) {
                    tnode[i_depth].nodep->makeTerminal(tnode[i_depth].mean_node(),tnode[i_depth].n);
                    if(i_depth==0) { //if main node is a terminal 
                        break; //stop growing this tree
                    } else {
                        i_depth--;  //some other not-main node was terminal
                        continue;   //go up one level
                    }
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
                    //idx = index->get_col_p(tnode[i_depth].nodep->bestvar);
                    splitval = tnode[i_depth].nodep->splitval;
                    
                    //x____test_n=0;
                    //x____test_sum=0;
                    //int x_test_innode=0;
                    for(uint16_t i=0; i<X->get_row(); i++) {
                        in_child[i] = in_parent[i] && (x[i] >= splitval);
/*                         if(in_child[i]) {
                            x____test_n++;
                            x____test_sum += y->at(i);
                        } */
                        //if(in_parent[i]) x_test_innode++;
                    }
                    
                    /* if(
                        x____test_n != tnode[i_depth].lchild_n ||
                        Abs(x____test_sum - tnode[i_depth].lchild_sum)>0.001
                    ) {
                      //  sprintln(x_test_innode);
                        sprintln(Abs(x____test_sum - tnode[i_depth].lchild_sum));
                        sprintln(tnode[i_depth].nodep->bestvar);
                        sprint(" lsought:");sprint(tnode[i_depth].lchild_n);
                        sprint(" rsought:");sprint(tnode[i_depth].rchild_n);
                        sprint(" lfound: ");sprintln(x____test_n);
                        error("wrong left n");
                    }
 */   

                     //go to left child
                    tnode[i_depth+1].status=0; //reset child status(n_visits) to zero
                    i_depth++;
                    continue;
                }

            } else {
                if(tnode[i_depth].status == 2) {

                    //temp_node tmp_l = tnode[i_depth+1];
                    //recover right node pointer with ugly pointer arithmetics
                    tnode[i_depth+1].nodep = nodes + tnode[i_depth].nodep->right_child_id;
                    tnode[i_depth+1].sum = tnode[i_depth].rchild_sum;
                    tnode[i_depth+1].n   = tnode[i_depth].rchild_n;

                    //update right child temp_node innode by best split
                    in_parent = tnode[i_depth  ].innodep;
                    in_child  = tnode[i_depth+1].innodep;
                    x = X->get_col_p(tnode[i_depth].nodep->bestvar);
                    //idx = index->get_col_p(tnode[i_depth].nodep->bestvar);
                    splitval = tnode[i_depth].nodep->splitval;
                    //x____test_n=0; //debug only
                    for(uint16_t i=0; i<X->get_row(); i++) {
                        in_child[i] = in_parent[i] && x[i] < splitval; //random tie handling also here?
                       // if(in_child[i]) x____test_n++; //debug only
                    }
                    
                    //debug
                    /* if(x____test_n != tnode[i_depth].rchild_n) {
                        sprint(" lsought:");sprint(tnode[i_depth].lchild_n);
                        sprint(" rsought:");sprint(tnode[i_depth].rchild_n);
                        sprint(" rfound: ");sprintln(x____test_n);
                        error("wrong right n");
                     }

                    
                    temp_node tmp_r = tnode[i_depth+1];

                    sprint("children mean:");
                    sprint((tmp_l.mean_node() * tmp_l.n + tmp_r.mean_node() * tmp_r.n) / (tmp_l.n+tmp_r.n) );
                    sprint("parent mean:");
                    sprintln(tnode[i_depth].mean_node()); */
                    
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
                        /* tnode[i_depth].nodep->print_header();
                        tnode[i_depth].nodep->print_body();
                        tnode[i_depth].nodep->print_terminal(); */
                        sprintln(i_depth);
                        sprintln( tnode[i_depth].status);
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



void forest_lake::rec_grow(dynamic_array<float,uint16_t>* newX, dynamic_vector<float,uint16_t>* newy) {
    
        //run time pars
    const uint16_t p_depth = 3;
    const uint16_t p_minnode = 5;
    const uint16_t p_ntree = 150;
    const uint16_t p_sampsize = 200;
    
    
    X = newX;
    y = newy;
    
    //allocate on heap, p_sampsize;
    dynamic_array<uint16_t,uint16_t> SEG(p_sampsize,1);
    dynamic_vector<uint16_t,uint16_t> vec = SEG.get_vector(0);
  
    while(two_more_nodes() && i_tree<p_ntree) {
        //grow tree
        //sprint("inode479");sprintln(i_node);
        new_tree2(); //tree tag marks a new tree starts in forest lake
        //sprint("inode482");sprintln(i_node);

        //draw random sample with replacement into seg(ment) and compute target sum of nodes
        temp_node_sum_left=0;
        for(auto& i : vec) {
            i = uint_dist(rng) % p_sampsize;
            temp_node_sum_left += y->at(i);
        }
        //sprint("inode489 ");sprintln(i_node);
        node* root_nodep = new_node3(-42);
        /* sprint("inode491 ");sprintln(i_node);
        sprint("addrdiff ");sprintln(int(root_nodep-nodes));
        sprintln(int(&nodes));
        sprintln(int(&nodes[0]));
        sprintln(int(&nodes[1]));
        sprintln(int(root_nodep)); */

       // for(auto& i : vec) {sprint(i);sprint(", ");}
        
        //grow this tree recursively
        grow_node(vec.begin(),vec.end(),root_nodep,1);
    }
}

void forest_lake::grow_node(uint16_t* Sp, uint16_t* Ep, node* parent_node, uint16_t depth) {
    //sprint("\n entering node: ");sprintln(int(parent_node-nodes));
    //sprint("\n i_node: ");sprintln(i_node);
    //sprint('\n');
    //sprint(depth);
    //either terminate this node and return...
    uint16_t* Cp{Sp+1};
    int n_parent{Ep-Sp};
    
    if(
         n_parent<=7 || depth>=9 || !two_more_nodes()  ||   //if this node should no be tried splitted
        !recsplit(Sp,Ep,Cp,parent_node)                       //or if split failed... (may fail if all feature values are the same)
    ) {
        
        float predSum{0};
        uint16_t* i=Sp;
        while(i!=Ep) {
            predSum += y->at(*i);
            i++;
        }
        parent_node->makeTerminal(predSum/(n_parent),n_parent); // "save to forest_lake_pointer"
        return;
    }
    //otherwise allocate new nodes, split and jump to child nodes
    
    
    parent_node->right_child_id = i_node;
    node* right_node = new_node3(Cp-Sp);
    node* left_node  = new_node3(Ep-Cp);
  
    grow_node(Sp,Cp,right_node,depth+1);
    grow_node(Cp,Ep,left_node ,depth+1);
   
}
    
bool forest_lake::recsplit(
    uint16_t* Sp, uint16_t* Ep, uint16_t*& Cp, node* parent_node
) {
    //sprint("\n new split:");sprintln(i_node);

    //internal variables
    uint16_t parent_n = uint16_t(Ep-Sp); //is safe as n_obs is limited
    float crit{-1}, critmax{-1.0}, prev_y{0}, sum_r{0}, sum_l{0}, best_sum{0}, best_splitval{0};
    uint16_t tieVal{0}, prev{0}, curr{0}, n_r{0}, n_l{0}, best_var{0};
    const uint16_t n_col{X->get_col()}; //dummy init
    uint16_t* i_Cp{nullptr};
    
    constexpr uint16_t p_minsplit=2;
    if(p_minsplit*2>=parent_n) {
        sprint(parent_n);
        error("p_minsplit is set too high");
    }

    float* yp = y->begin();

    float predSum{0};
    {
        //sprint(" ps");
        uint16_t* i=Sp;
        while(i!=Ep) {
            predSum += yp[*i];
            i++;
        }
    }

    //iterate each variable
    for(uint16_t i_var=0; i_var < n_col;i_var++) {

        i_Cp = Sp;
        const auto* x = X->get_col_p(i_var); //reference to column for this variable in matrix
        std::sort(Sp,Ep,[x](size_t i1, size_t i2) {return x[i1] < x [i2];}); //sort by x feature column
    
        //set sums and counters before rolling mean square error crit
        n_r = p_minsplit;
        n_l = parent_n-p_minsplit;
        sum_r=0;
        sum_l=predSum;//parent_node->splitval;
        
        for(uint16_t i =0;i<p_minsplit;i++) {
            prev    = *i_Cp;
            prev_y  = yp[prev];
            sum_r  += prev_y;
            sum_l  -= prev_y;
            ++i_Cp;
        }        
        curr = *i_Cp;
        //iterate all splits and compute rolling mean square error crit
        while(i_Cp!=Ep) {

            //skip if split on same obs
            if(curr!=prev && x[curr]!=x[prev]) {
                //compute loss/crit
                crit = (sum_l * sum_l / n_l) + (sum_r * sum_r / n_r);// - crit_parent;

                //if critmax save
                if (crit >= critmax) {  //handle better crit, accept new split
                    if (crit!=critmax) tieVal=0;
                    if (crit!=critmax || uint_dist(rng) > UINT32_MAX / ++tieVal) { //update
                        Cp  = i_Cp;
                        best_var = i_var;
                        best_sum = sum_l;
                        best_splitval = (x[curr]+x[prev])/2;
                        //sprint("bs");sprint(best_splitval);sprint(" v ");sprint(i_var);sprint(" c ");sprintln(crit);
                    }
                    critmax = crit;
                }            
            }
            //move the lowest feature value obs from left node to right node, and update counts and sums
            if(--n_l<p_minsplit) break; //
            n_r++;
            prev = curr;
            ++i_Cp;
            curr = *(i_Cp);
            prev_y = yp[prev];
            sum_r += prev_y;
            sum_l -= prev_y;
        }
    }

    //parent_node->right_child_id = i_node;
    parent_node->bestvar = best_var;
    //temp_node_sum_left = best_sum;
    //temp_node_sum_right = parent_node->splitval - best_sum;
    //sprint(" sv ");
    //sprint(int(Cp-Sp));
    
    //resort by best split
    const auto* x = X->get_col_p(best_var); //reference to column for this variable in matrix
    std::sort(Sp,Ep,[x](size_t i1, size_t i2) {return x[i1] < x [i2];}); //sort by x feature column
    parent_node->splitval = (X->at(*(Cp-1),best_var) + X->at(*Cp,best_var))/2;
    
    parent_node->splitval = best_splitval;
    
    //if split success Cp will not point as Sp
    return Cp!=Sp+1;
}


float forest_lake::predict(dynamic_array<float,uint16_t>* X,uint16_t i_row){
    
    //copy 
    float xnew[X->get_col()];
    for(uint16_t i_col=0;i_col<X->get_col();i_col++) xnew[i_col] = X->get(i_row,i_col);
    //sprint("x0=");sprint(xnew[0]);sprint("|");
    uint16_t j_tree=0;   //count tree predictions found
    float prediction=0;  //sum predictions
    uint16_t k_node = 0; //iterate nodes in one tree

    for(uint16_t j_node=0; j_node < n_nodes; j_node++) { //iterate all nodes in forest lake to find trees
        if(nodes[j_node].right_child_id==0 &&  nodes[j_node].bestvar==0) { //tree found
            if(j_tree>=i_tree) break;
          
            k_node = j_node + 1; //iterate tree

            while(true) {

                node& this_node = nodes[k_node];
                if(this_node.right_child_id!=0) {
                    
                    //this node is intermediary
                    k_node = this_node.right_child_id;
                    if(xnew[this_node.bestvar] >= this_node.splitval) {
                      k_node++;
                    }

                } else {
                    
                    //this node is terminal
                    j_tree++;
                    prediction += this_node.splitval;
                    j_node = k_node; //new tree cant possibly start before k_node
                    break;
                }
            }
        }
    }

    return(prediction/j_tree);
}

 void forest_lake::predict_all(dynamic_array<float,uint16_t>* X,dynamic_array<float,uint16_t>* out,uint16_t i_col){
    if(out->get_row()!=X->get_row()) error("out and X must have same row size");
    if(out->get_col()<=i_col) error("i_col exceeds out number of columns");
    const uint16_t n_obs = uint16_t(X->get_row());
    dynamic_vector<float,uint16_t> vec = out->get_vector(i_col);
    uint16_t j=0;
    for(auto& i : vec) i = forest_lake::predict(X,j++);
}