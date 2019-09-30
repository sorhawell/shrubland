
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
            if(!pushRow(value)) error("overflow:");
        }
        nextRow();
    }
};

