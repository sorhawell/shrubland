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
    Serial.print(right_child_id);Serial.print("\t\t");
    Serial.print(bestvar);Serial.print("\t\t\t");
}
void node::print_terminal(){
    Serial.print("\t\t\t\t\t");Serial.print(splitval);
    Serial.print("\t"); Serial.print(bestvar);Serial.print("\t");
}

void forest_lake::print_nlines(uint16_t nlines) {
//for(uint16_t i=0; i<150;i++) cout << sf.forest.tree_index[i] << ", ";
    
    uint16_t j_tree=0;
    //for(uint16_t i=0; i<sf.forest.i_node;i++) {
        Serial.println();
        nodes[0].print_header();
    for(uint16_t i=0; i<nlines && i<n_nodes;i++) {

        if(nodes[i].right_child_id==0 && nodes[i].bestvar==0 ) {    
            j_tree++;
        }
        Serial.println(); Serial.print(i); Serial.print(" \t");
        if(nodes[i].right_child_id!=0) nodes[i].print_body();
        if(nodes[i].right_child_id==0 && nodes[i].bestvar!=0) nodes[i].print_terminal();
        Serial.print(j_tree);
    } 
}

void forest_lake::print_splits(uint16_t nlines) {
    uint16_t j_tree=0;
    uint16_t depth=0;
    uint16_t max_depth=15;
    Serial.print(j_tree);
    for(uint16_t i=0; i<nlines && i<n_nodes;i++) {
        if(nodes[i].right_child_id==0 && nodes[i].bestvar==0 ) {    
            if(nodes[i+1].right_child_id==0 && nodes[i+1].bestvar==0) break;
            depth = 0;
            Serial.print("\n t");Serial.print(++j_tree);Serial.print(" ");
        }
        if(++depth<=max_depth) {Serial.print(nodes[i].splitval);Serial.print(',');}
    } 
};

void node::makeTerminal(float pred_value, uint16_t terminal_n) {
    right_child_id = 0;
    bestvar = terminal_n;
    splitval = pred_value;
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
        Serial.println("-");
        Serial.println(snode->lchild_n);
        Serial.println(snode->rchild_n);
        Serial.println(snode->nodep->bestvar);
        Serial.println(snode->sum);
        Serial.println(snode->n);
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
    const uint16_t p_depth = 10;
    const uint16_t p_minnode = 5;
    const uint16_t p_ntree = 100;
    const uint16_t p_sampsize = 350;
    
    
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
    uint16_t innode_sz{p_depth*X->get_row()};
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
                      //  Serial.println(x_test_innode);
                        Serial.println(Abs(x____test_sum - tnode[i_depth].lchild_sum));
                        Serial.println(tnode[i_depth].nodep->bestvar);
                        Serial.print(" lsought:");Serial.print(tnode[i_depth].lchild_n);
                        Serial.print(" rsought:");Serial.print(tnode[i_depth].rchild_n);
                        Serial.print(" lfound: ");Serial.println(x____test_n);
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
                        Serial.print(" lsought:");Serial.print(tnode[i_depth].lchild_n);
                        Serial.print(" rsought:");Serial.print(tnode[i_depth].rchild_n);
                        Serial.print(" rfound: ");Serial.println(x____test_n);
                        error("wrong right n");
                     }

                    
                    temp_node tmp_r = tnode[i_depth+1];

                    Serial.print("children mean:");
                    Serial.print((tmp_l.mean_node() * tmp_l.n + tmp_r.mean_node() * tmp_r.n) / (tmp_l.n+tmp_r.n) );
                    Serial.print("parent mean:");
                    Serial.println(tnode[i_depth].mean_node()); */
                    
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
                        Serial.println(i_depth);
                        Serial.println( tnode[i_depth].status);
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
    //Serial.print("x0=");Serial.print(xnew[0]);Serial.print("|");
    uint16_t j_tree=0;   //count tree predictions found
    float prediction=0;  //sum predictions
    uint16_t k_node = 0; //iterate nodes in one tree

    for(uint16_t j_node=0; j_node < n_nodes; j_node++) { //iterate all nodes in forest lake to find trees
        if(nodes[j_node].right_child_id==0 &&  nodes[j_node].bestvar==0) { //tree found
            if(j_tree>=i_tree) break;
          
            k_node = j_node + 1; //iterate tree
            //Serial.println("k_node init:");Serial.println(k_node);
            while(true) {
                /* nodes[k_node].print_body();
                Serial.println(" |"); */
                if(nodes[k_node].right_child_id!=0) {
                    
                    //this node is intermediary
                    if(xnew[nodes[k_node].bestvar] >= nodes[k_node].splitval) {
                      k_node = nodes[k_node].right_child_id + 1 ;
                    } else {
                      k_node = nodes[k_node].right_child_id;
                    }
                    
                    
                //Serial.print("[");Serial.print(k_node);Serial.print("]");
                } else {
                    
                    //this node is terminal
                    j_tree++;
                    prediction += nodes[k_node].splitval;
                   

                    //Serial.print(", (");Serial.print(nodes[k_node].bestvar);Serial.print(")");
                    //Serial.print(nodes[k_node].splitval);
                    j_node = k_node; //new tree cant possibly start before k_node
                    break;
                }
            }
        }
    }
    //Serial.print(" = ");Serial.println(prediction/j_tree);
    return(prediction/j_tree);
}