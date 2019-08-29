#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <chrono>
#include <iomanip>
#include "Profileapi.h"

using namespace std;


void error(string e) {throw runtime_error(e);}


template <class T, class S>
class dynamic_vector {
    T* dyn_vec;
    T* j;
    S sz;
    

    public:
    dynamic_vector(T* begin, S size):
    sz(size), dyn_vec(begin), j(begin){};
    
    T at(S i){
        if(i>=sz) error("bound exceeded");
        return(*(dyn_vec+i));
    }

    T get(S i){
        return(*(dyn_vec+i));
    }

    T* get_cell_ref() {
        return(dyn_vec);
    }

    T next() {
        if(j-dyn_vec>=sz) error("bound exceeded");
        *(j++);
    }

};

template <class T, class S>
class dynamic_array{
    
    T* dyn_array;
    T* j;
    S n_row;
    S n_col;
    S sz;
    S j_col;
    S j_row;
    

    public:
    dynamic_array(T* begin, S size, S ncol):
    dyn_array(begin), n_row(size/ncol), n_col(ncol), sz(size), j_col(0), j_row(0), j(begin) {};

    T at(S i_row,S i_col) {
        if(i_row>=n_row||i_col>=n_col) error("bounds exceeded");
        return(*(dyn_array+i_row+(i_col*n_row)));
    }

    T get(S i_row,S i_col) {
        if(i_row>=n_row||i_col>=n_col) error("bounds exceeded");
        return(*(dyn_array+i_row+(i_col*n_row)));
    };

    T next() {*(j++);}
    T prev() {*(--j);}

    T* get_col_p(S i_col) {
        return(dyn_array+(i_col*n_row));
    }

    dynamic_vector<T,S> get_vector(S i_col) {
        dynamic_vector<T,S> dyn_vec(dyn_array+i_col*n_row,n_row);
        return(dyn_vec);
    }
    void print() {
        cout << "\nnrow:" << n_row << "ncol: " << n_col << "size: " << sz << '\n';
    }
};



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

void   CPUtimer::setFreq() {QueryPerformanceFrequency(&frequency);}
void   CPUtimer::start()   {QueryPerformanceCounter(&t1);}
void   CPUtimer::reset_sumTime() {sumTime=0;}
double CPUtimer::stop()  {
    QueryPerformanceCounter(&t2);
    elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000000.0 / frequency.QuadPart;
    sumTime += elapsedTime;
    return(elapsedTime);
}    
CPUtimer Timer;

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
uint32_t seed_val;           // populate somehow
MyRNG rng;                   // e.g. keep one global instance (per thread)
void initialize() {
  rng.seed(seed_val);
}
std::uniform_int_distribution<uint32_t> uint_dist;

//define a matrix struct for storing data and accessing during training and prediction
template <class T,int Tnrow, int Tncol> struct matrix{
    T a[Tncol][Tnrow];
    int nrow, ncol, size, i_row, i_col;
    bool overflow;

    matrix():
    nrow(Tnrow), ncol(Tncol), size(Tnrow*Tncol),
    i_row(0), i_col(0), overflow(false) {};

    matrix(int pnrow, int pncol):
    nrow(pnrow), ncol(pncol), size(pnrow*pnrow),
    i_row(0), i_col(0), overflow(false) {};

    bool nextRow();
    bool pushRow(T* x);
    bool pushRow(T  x);
    bool push   (T* x);
    bool push   (T  x);
    void fillCSV(string path);
    void fillCSV_float(string path); //could be reworked to template specialization
    void truncate();
    void print();
    void resize(int rows, int cols);
    void fill(T value);
    void erase_col(int from,int to);

};

template <class T,int Tnrow, int Tncol>
void matrix<T,Tnrow,Tncol>::truncate() {
    i_col=0;
    i_row=0;
    overflow =false;
}

template <class T,int Tnrow, int Tncol>
void matrix<T,Tnrow,Tncol>::resize(int rows, int cols) {
    if(rows>Tnrow) error("number of rows exceeds allocated space of matrix");
    if(cols>Tncol) error("number of cols exceeds allocated space of matrix");
    ncol = cols;
    nrow = rows;
}

template <class T,int Tnrow, int Tncol>
void matrix<T,Tnrow,Tncol>::fill(T value) {
    for(int i_col=0;i_col<ncol;i_col++) {
        for(int i_row=0;i_row<nrow;i_row++) {
            a[i_col][i_row] = value;
        }
    }
}

template <class T,int Tnrow, int Tncol>
void matrix<T,Tnrow,Tncol>::erase_col(int from, int to) {
    int n_col_erased = to - from;
    for(int i_col_copy=to;i_col_copy<ncol;i_col_copy++) {
        for(int i_row_copy=0;i_row_copy<nrow;i_row_copy++) {
            a[from][i_row_copy] = a[i_col_copy][i_row_copy];
        }
        from++;
    }
    ncol = ncol - n_col_erased;
}

template <class T,int Tnrow, int Tncol>
bool matrix<T,Tnrow,Tncol>::nextRow() {
    i_col=0;
    i_row++;
    if(i_row>=nrow ||i_row<0) overflow=true; else overflow = false;
    return(!overflow);
}

template <class T,int Tnrow, int Tncol>
bool matrix<T,Tnrow,Tncol>::pushRow(T* x) {
    if(i_col>=ncol || i_col<0 ) overflow=true;
    if(!overflow) a[i_col++][i_row] = *x;
    return(!overflow);
}

template <class T,int Tnrow, int Tncol>
bool matrix<T,Tnrow,Tncol>::pushRow(T x) {
    if(i_col>=ncol || i_col<0 ) overflow=true;
    if(!overflow) a[i_col++][i_row] = x;
    return(!overflow);
}

template <class T,int Tnrow, int Tncol>
bool matrix<T,Tnrow,Tncol>::push(T* x) {
    if(!pushRow(x)) {
        nextRow();
        if(overflow) return(false);
        if(!pushRow(x)) error("matrix::push failed to push");
    }
    return(!overflow);
}

template <class T,int Tnrow, int Tncol>
bool matrix<T,Tnrow,Tncol>::push(T  x) {
    if(!pushRow(x)) {
        nextRow();
        if(overflow) return(false);
        if(!pushRow(x)) error("matrix::push failed to push");
    }
    return(!overflow);
}

template <class T,int Tnrow, int Tncol>
void matrix<T,Tnrow,Tncol>::print() {
    cout << "\ndim [" << nrow << "]x[" << ncol << "]\n";
    for(int i=0 ; i<nrow; i++) {
        for(int j=0; j<ncol; j++) {
            cout << a[j][i] << "\t, ";
        }
        cout << '\n';
    } 
    
}

//handle different conversions
double stox(double* pd, string str) {return(stod(str));}
double stox(float* pd, string str) {return(stof(str));}
double stox(int* pd, string str) {return(stoi(str));}

template <class T,int Tnrow, int Tncol>
void matrix<T,Tnrow,Tncol>::fillCSV(string path) {
    std::ifstream data(path);
    std::string line;
    while(std::getline(data,line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while(std::getline(lineStream,cell,',')) {
            //cout << cell << " - [" << i_row << ' ' << i_col << "] \n";
            T value = stox(&value,cell);
            if(!pushRow(value)) error("overflow:"+to_string(i_row)+","+to_string(i_col));
        }
        nextRow();
    }
};

//data structures for shrubs
struct node_stack {
    node_stack() :
    crit_parent(0), critmax(0), parent_sum(0), 
    parent_n(0), innode(nullptr), splitcurr(-1) {
    };
    double crit_parent, critmax, parent_sum;
    int parent_n, splitcurr;
    bool* innode;
    void print();

};
void node_stack::print() {
    cout << '\n' <<
      "crit_parent =" << crit_parent <<
    "\ncritmax     =" << critmax     <<
    "\nparent_sum  =" << parent_sum  <<
    "\nparent_n    =" << parent_n    <<
    "\ninnode point=" << &innode     << '\n';
}

struct node {
    node() : splitval(0), bestvar(0), right_child_id(0) {};
    double splitval;
    int bestvar, right_child_id;
    void makeTerminal(node_stack&);
    void print();
    void print_body();
    void print_header();
    
};
void node::print() {
    cout << "value \t right_child_id \t bestvar" <<
    "nsplitval     =" << splitval          <<
    "\nright_child_id=" << right_child_id    <<
    "\nbestvar       =" << bestvar           <<'\n';
}
void node::print_header() {cout << "\nsplitval \tright_child_id \tbestvar \ttree";}
void node::print_body()   {cout << '\n'<<splitval<<"\t\t"<< right_child_id<<"\t\t"<< bestvar;}
void node::makeTerminal(node_stack& ns) {
    right_child_id = 0;
    splitval = ns.parent_sum / ns.parent_n;
}




template<int Tntree, int Tnnode>
struct Forest {
    int n_nodes;
    int n_trees;
    int i_node;
    int n_nodes_i_tree;
    int i_tree;
    node nodes[Tnnode];
    int tree_index[Tntree];
    
    Forest() : n_nodes(Tnnode), n_trees(Tntree), i_node(0), n_nodes_i_tree(0), i_tree(0) {
        for(int i; i<Tntree;i++) {tree_index[i]=0;}
    };
    node& new_node();
    node& new_tree();
    bool more_trees();
    bool more_nodes();
    bool two_more_nodes();
    void truncate();
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
    
    shrub() 
    : sdata(nullptr), X(nullptr), y(nullptr), index(nullptr), ybag(nullptr),
    y_size(0), p_ntree(Tntree), p_mtry(0), p_maxdepth(10), p_maxnodes(50), p_minnode(5), data_ready(false) {
    };
    
    
    
    void make_innode2(
        const node& this_node, const node_stack& this_stack_node,
        node_stack& stack_l_node, node_stack& stack_r_node);
    void bestsplit(node& this_node, node_stack& this_stack_node);
    void eval_node(node& this_node, node_stack& this_stack_node, int depth);
    void grow_forest();
    void link_to_data(shrub_data<T,Tnrow,Tncol>*);
    T predict(T* newX);
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

/* template <class T, int Tntree, int Tnnode, int Tnrow, int Tncol>
T shrub<T,Tntree,Tnnode,Tnrow,Tncol>::predict(T* newX) {
    Xtest = newX;
    yresult = 0;
    node* pn = nullptr;
    
    pn->   
    for(int i=0;i<forest.i_tree;i++) {
        pn = &forest.nodes[forest.tree_index[i]];
        while(pn->right_child_id!=0) {
            if(Xtest[pn->bestvar]<=pn->splitval) {
                pn = pn+1; // go to left node
            } else {
                pn = &forest.nodes[pn->right_child_id]; //go to right node
            }
        }
        yresult += pn->splitval;
    }
    yresult /= forest.i_trees:
    return(yresult);
};
 */



constexpr int ntree=50;
constexpr int nnode=50*300;
constexpr int nrow =10000;
constexpr int ncol =6;
constexpr int ycol =5;

int main() {

    Timer.setFreq();

    int nval = nrow*ncol;
    double Xbuffer[nval];
    for(int i=0;i<nval;i++) Xbuffer[i]=i;
    dynamic_array<double,int> X(Xbuffer,nval,ncol);
    X.print();
    
    //safe access
    for(int i=0; i<5; i++) {cout << " " << X.at(i,0);}

    Timer.reset_sumTime();
    Timer.start();
    double temp;
    for(int i=0; i<nrow; i++) { temp =+ X.at(i,4);}
    Timer.stop();
    cout << "\n safe time: " << Timer.sumTime;

    //get raw pointer at own risk
    cout << "\n x4 ";
    double* x4 = X.get_col_p(4);
    for(int i=0; i<5; i++) {cout << " " << x4[i];}

    Timer.reset_sumTime();
    Timer.start();
    double temp2;
    for(int i=0; i<nrow; i++) { temp2 =+ x4[i];}
    Timer.stop();
    cout << "\n x4 time: " << Timer.sumTime;

    
    //get raw pointer at own risk
    cout << "\n v4 ";
    dynamic_vector<double,int> v4 = X.get_vector(4);
    for(int i=0; i<5; i++) {cout << " " << v4.at(i);}


    Timer.reset_sumTime();
    Timer.start();
    double temp3;
    for(int i=0; i<nrow; i++) { temp3 =+ v4.at(i);}
    Timer.stop();
    cout << "\n safe v4.at(i) time: " << Timer.sumTime;

    Timer.reset_sumTime();
    Timer.start();
    double temp4;
    for(int i=0; i<nrow; i++) { temp4 =+ v4.next();}
    Timer.stop();
    cout << "\n v4 next time: " << Timer.sumTime;


    /* int a;
    int b;

    cin >>a >>b;

    int arr[a][b];

    arr[a-1][b-1]=0;

    cout <<"l604\n";    
    shrub<double,ntree,nnode,nrow,ncol> sf; //initialize a forest with 10 trees and 500 nodes(max)
    cout <<"l606\n";    

    Timer.setFreq();
    matrix<double,nrow,ncol> X;
    matrix<double,nrow,1   > y;
    X.fillCSV("./data/train1.csv");
    //X.print();
    
    for(int i = 0 ; i<X.nrow; i++) {y.a[0][i] = X.a[ycol][i];}
    X.erase_col(ycol,ycol+1);
    
    shrub_data<double,nrow,ncol> sd(X,y);
    sd.index_features();
    sf.link_to_data(&sd);
    
    sf.p_maxnodes = 250;
    sf.p_maxdepth = 12;
    sf.p_minnode = 5;
    cout << "hello";
    Timer.start();
    sf.grow_forest();
    Timer.stop();


    

    cout << "still alive time avg=" << Timer.sumTime << "\n";

    for(int i=0; i<sf.forest.n_trees;i++) cout << sf.forest.tree_index[i] << ", ";
    
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    //for(int i=0; i<1000;i++) {sf.forest.nodes[i].print_body();}
     
    //for(int i=0; i<150;i++) cout << sf.forest.tree_index[i] << ", ";
    sf.forest.nodes[0].print_header();
    int j_tree=0;
    //for(int i=0; i<sf.forest.i_node;i++) {
    for(int i=0; i<250;i++) {
        if(i==sf.forest.tree_index[j_tree]) {
            j_tree++;
            cout << "\n";
        }
        sf.forest.nodes[i].print_body();
        cout << "\t\t" << j_tree;
    }   */
}
