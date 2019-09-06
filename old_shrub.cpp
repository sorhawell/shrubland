#include "old_shrub.h"

template<int Tntree, int Tnnode>
Forest<Tntree, Tnnode>::Forest() : n_nodes(Tnnode), n_trees(Tntree), i_node(0), n_nodes_i_tree(0), i_tree(0) {
        for(int i; i<Tntree;i++) {tree_index[i]=0;}
};

template<int Tntree, int Tnnode>
node& Forest<Tntree,Tnnode>::new_node() {
    if(i_node==n_nodes) error("no more nodes");
    n_nodes_i_tree++;
    return(nodes[i_node++]);
}
template<int Tntree, int Tnnode>
node& Forest<Tntree,Tnnode>::new_tree() {
    if(i_tree==Tntree) error("no more trees");
    n_nodes_i_tree=1;
    tree_index[i_tree++] = i_node;
    return(new_node());
}
template<int Tntree, int Tnnode>
bool Forest<Tntree,Tnnode>::more_trees() {
    return(i_tree<Tntree);
}
template<int Tntree, int Tnnode>
bool Forest<Tntree,Tnnode>::more_nodes() {
    return(i_node<Tnnode);
}
template<int Tntree, int Tnnode>
bool Forest<Tntree,Tnnode>::two_more_nodes() {
    return(i_node < Tnnode-1);
}
template<int Tntree, int Tnnode>
void Forest<Tntree,Tnnode>::truncate() {
    i_tree, i_node = 0;
    for(int i; i<Tntree;i++) {tree_index[i]=0;}
}




node_stack::node_stack() :
    crit_parent(0), critmax(0), parent_sum(0), 
    parent_n(0), innode(nullptr), splitcurr(-1) {
};

void node_stack::print() {
    std::cout << '\n' <<
      "crit_parent =" << crit_parent <<
    "\ncritmax     =" << critmax     <<
    "\nparent_sum  =" << parent_sum  <<
    "\nparent_n    =" << parent_n    <<
    "\ninnode point=" << &innode     << '\n';
}
 

template <class T, int Tnrow, int Tncol>
struct shrub_data {
  //stack varibles
    matrix<T,Tnrow,Tncol>& X;
    matrix<T,Tnrow,1>& y;

    matrix<uint16_t,Tnrow,Tncol> index;
    matrix<uint16_t,Tnrow,1> ybag;
    
    shrub_data(matrix<T,Tnrow,Tncol>& newX, matrix<T,Tnrow,1>& newy) : 
    X(newX),y(newy) {
        
    };
//    X(nullptr) ,y(nullptr) {};

    void index_features();
    void set_X(matrix<T,Tnrow,Tncol>& X_param);
    void set_y(matrix<T,Tnrow,1>&     y_param);
};

template <class T, int Tnrow, int Tncol>
shrub_data<T,Tnrow,Tncol>::shrub_data(matrix<T,Tnrow,Tncol>& newX, matrix<T,Tnrow,1>& newy) : 
    X(newX),y(newy) {        
};

template <class T, int Tnrow, int Tncol>
void shrub_data<T,Tnrow,Tncol>::set_X(matrix<T,Tnrow,Tncol>& X_param) {
    X = X_param;
}

template <class T, int Tnrow, int Tncol>
void shrub_data<T,Tnrow,Tncol>::set_y(matrix<T,Tnrow,1>& y_param) {
    y = y_param;
}

template <class T, int Tnrow, int Tncol>
void shrub_data<T,Tnrow,Tncol>::index_features() {
    
    index.resize(X.nrow,X.ncol);
    index.truncate();
    //index.fill(0);
    //for(auto& i : index.a) iota(i.begin(),i.end(),0);

    for(int i=0 ; i<X.ncol; i++) {
        for(int j=0; j<index.nrow;j++) index.a[i][j] = j;
        const T* x = &(X.a[i][0]);
        uint16_t* pidx = &(index.a[i][0]);
        sort(pidx, (pidx+index.nrow),
            [x](size_t i1, size_t i2) {return x[i1] < x[i2];}
        );
    }
}

template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
shrub<T,Tntree,Tnnode,Tnrow,Tncol>::shrub()
    : sdata(nullptr), X(nullptr), y(nullptr), index(nullptr), ybag(nullptr),
    y_size(0), p_ntree(Tntree), p_mtry(0), p_maxdepth(10), p_maxnodes(50), p_minnode(5), data_ready(false) {
};

template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
void shrub<T,Tntree,Tnnode,Tnrow,Tncol>::link_to_data(shrub_data<T,Tnrow,Tncol>* new_sdata) {
    sdata = new_sdata;

    //load short hand pointers
    X = &sdata->X;
    y = &sdata->y;
    index = &sdata->index;
    ybag = &sdata->ybag;
    y_size = sdata->y.nrow;
    X_ncol = sdata->X.ncol;

    //check data 
    if(y->nrow != X->nrow) error("incompatible size of y and X");
    if(X->nrow != index->nrow) error("incompatible size of X and index");

    data_ready = true;
}


template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
void shrub<T,Tntree,Tnnode,Tnrow,Tncol>::bestsplit(node& this_node, node_stack& this_stack_node) {

    //internal variables
    double sum_l{0}, sum_r{0};   //rolling sum/substraction of target values
    int n_l{0}, n_r{0}, tieVal{0}; //count size of left and tight partitions
    double crit{0}, critvar{0};
    int prev{0}, curr{0}; //dummy init
    bool none_found_yet{true};
    
    
    const T* y = &sdata->y.a[0][0];

    for(int i_var=0; i_var < sdata->X.ncol;i_var++) { //iterate each variable
        
        const auto& X = sdata->X.a[i_var]; //reference to column for this variable in matrix
        

        //reset variables for rolling mse loss function
        sum_l = this_stack_node.parent_sum;
        sum_r = 0;
        n_l = this_stack_node.parent_n;
        n_r = 0;
     
        auto& i_var_indexed = sdata->index.a[i_var];
        none_found_yet = true;

        for(auto i_idx=0; i_idx<sdata->index.nrow;i_idx++) {
            //use index to iterate any split point sorted from low variable value to high
            curr = i_var_indexed[i_idx];
            //std::cout << "\n v: " << i_var << " c: " << curr << " p: " << prev;
             
            
            if(!  this_stack_node.innode[curr]) {   //check if indexed obs is innode
                //std::cout << "skip because as not innode \n";
                continue;
            }
            if(none_found_yet) {//check if first found innode
                none_found_yet = false;
                prev = curr;
                //std::cout << "skip because first \n";
                continue;
            }
            if(X[curr]==X[prev]) {//check if same value
                //std::cout << "skip because same variable value \n";
                continue;
            }

            //compute crit
            sum_l -= y[curr];
            sum_r += y[curr];
            crit = (sum_l * sum_l / --n_l) + (sum_r * sum_r / ++n_r);// - crit_parent;            
            //std::cout << "crit: " << crit << '\n';
            
            //update crit handling
            if (crit == this_stack_node.critmax) { //handle tie, random accept new split
                //std::cout << "\n tie handling - ";
                tieVal++;
                uint32_t randVal = uint_dist(rng);
                if (randVal < UINT32_MAX / tieVal) {
                    this_node.splitval = (X[curr] + X[prev]) / 2.0;
                    this_node.bestvar = i_var;
                    //this_node.splitcurr = curr;
                    //std::cout << " tie win \n";
                }
            } else if (crit > this_stack_node.critmax) {  //handle better crit, accept new split
                this_node.splitval = (X[curr] + X[prev]) / 2.0;
                this_node.bestvar = i_var;
                //this_node.splitcurr = curr;
                tieVal = 1;
                this_stack_node.critmax = crit;
                //std::cout << "\n new critvar: " << critvar << '\n';
            }

            //update indexes
            prev = curr;
        }
    }

}

//nonn index version
template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
void shrub<T,Tntree,Tnnode,Tnrow,Tncol>::make_innode2(
    const node& this_node, const node_stack& this_stack_node,
    node_stack& stack_l_node, node_stack& stack_r_node
) {
    
    //i'th obs is in child l_innode if i'th obs is this node (this_node.innode[i]=true) and not passed splitval
    bool smaller_than_split{false};
    T* x = sdata->X.a[this_node.bestvar];
    T* y = sdata->y.a[0];

    bool to_left{false}, to_right{false}, in_here{false};

    for(int i=0; i<y_size;i++) {
        in_here = this_stack_node.innode[i]; //here
        smaller_than_split = x[i] < this_node.splitval;  //opt: insert into own loop to find split indice
        to_left   = in_here &&  smaller_than_split;
        to_right  = in_here && !smaller_than_split;
        stack_l_node.innode[i] = to_left; 
        stack_r_node.innode[i] = to_right;
        if(to_left) {
            stack_l_node.parent_n++;
            stack_l_node.parent_sum += y[i];
        }
        if(to_right) {
            stack_r_node.parent_n++;
            stack_r_node.parent_sum += y[i];
        }

        //cout << '\n' << this_node.innode[i] << "-" << l_node.innode[i] << r_node.innode[i] <<
        //"\t" << l_node.parent_n << '-' << r_node.parent_n;
    }

}

template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
void shrub<T,Tntree,Tnnode,Tnrow,Tncol>::eval_node(node& this_node, node_stack& this_stack_node, int depth) {
    
    //cout << forest.n_nodes_i_tree << " ";
    //if(forest.n_nodes_i_tree>=p_maxnodes) cout << "\n >p_maxnodes" << forest.n_nodes_i_tree << " forest_i_node" << forest.i_node;

    //reasons to stop terminalize this node/every open node
    if(
        
        forest.n_nodes_i_tree>=p_maxnodes ||  // at/over limit of allowed nodes in each tree (terminate all nodes)
        this_stack_node.parent_n<p_minnode ||         // less than 5 obs in node (terminate this node)
        depth >= p_maxdepth ||                // depth of node limits the allowed (terminate this node)
        !forest.two_more_nodes()              // allocated memory in forest pool exhausted (terminate all nodes)
    ) {
        uint16_t flag=0;
        if(forest.n_nodes_i_tree>=p_maxnodes) flag=1;
        if(this_stack_node.parent_n<p_minnode) flag=2;       // less than 5 obs in node (terminate this node)
        if(depth >= p_maxdepth) flag=3;
        this_node.bestvar = this_stack_node.parent_n;



        this_node.makeTerminal(this_stack_node);
        return;   
    }

    ///if not returned yet here, this node will be splitted further
    //find for this node...
    shrub::bestsplit(this_node,this_stack_node);

    //reserve permanent nodes from Forest pool
    node& l_node = forest.new_node();
    node& r_node = forest.new_node();
    this_node.right_child_id = forest.i_node;
    
    //initialize transient nodes on stack
    node_stack stack_l_node;
    node_stack stack_r_node;
    bool l_innode[y_size];  //innodes index placed on stack, an hence deleted when finally leaving this node/scope
    bool r_innode[y_size];
    stack_l_node.innode = l_innode; // stack_l_node.innode pointers point to stack innode
    stack_r_node.innode = r_innode;
    
    //configure child nodes
    shrub::make_innode2(
        this_node,this_stack_node,
        stack_l_node,stack_r_node
    );

    //go evaluate child nodes
    eval_node(l_node, stack_l_node, depth+1);
    eval_node(r_node, stack_r_node, depth+1);

    //returned from every descending child node, return from this node also.
    return;
};

template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
void shrub<T,Tntree,Tnnode,Tnrow,Tncol>::grow_forest() {

    forest.truncate(); // reset/overwrite Forest pool
    const T* y = &sdata->y.a[0][0];
    uint16_t* ybag = &sdata->ybag.a[0][0];

    //grow trees until no more trees or nodes in Forest pool
    while(forest.more_trees() && forest.two_more_nodes()) {
    
    //reserve permenant node from Forest pool in shrubs object
    node& nd = forest.new_tree();
    

    //initialize transient node on stack (stack_nodes are deleted with the stack)
    node_stack stack_nd;   //object with all transient information
    bool inno[y_size];      //bool array, describing what obs are in this node
    stack_nd.innode = inno; //attach array to node_stack object

    //bootstrapping and configuration of first stack_nodes
    bool tempbool;


    
    for(int i=0; i<y_size; i++) {
        uint32_t randVal = uint_dist(rng);
        //cout << "_" << uint16_t(randVal%y_size);
        ybag[i] = uint16_t(randVal%y_size);
        
        stack_nd.parent_sum += y[ybag[i]];
        stack_nd.parent_n++;
        stack_nd.innode[i] = true;
    }
    
    //construct tree from first node
    eval_node(nd,stack_nd,1);

    }
};
