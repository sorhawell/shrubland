#include <iostream> 
using namespace std; 

// Class representing a reference counter class 
class Counter { 
public: 
    // Constructor 
    Counter(); 
    Counter(const Counter&) = delete; 
    Counter& operator=(const Counter&) = delete; 
  
    // Destructor 
    ~Counter() {}; 
  
    void reset(); 
    unsigned int get();
    
    // Overload post/pre increment 
    void operator++();
    void operator++(int); 
    
    // Overload post/pre decrement 
    void operator--();  
    void operator--(int); 
    
    // Overloading << operator 
    friend ostream& operator<<(ostream& os, const Counter& counter);
  
private: 
    unsigned int m_counter{}; 
}; 
  
// Class representing a shared pointer 
template <typename T> 
class Shared_ptr { 
public: 
    // Constructor 
    explicit Shared_ptr(T* ptr = nullptr);
  
    // Copy constructor 
    Shared_ptr(Shared_ptr<T>& sp);
  
    // Reference count 
    unsigned int use_count();
    // Get the pointer 
    T* get();
  
    // Destructor 
    ~Shared_ptr();
  
    friend ostream& operator<<(ostream& os, Shared_ptr<T>& sp) ;
  
private: 
    // Reference counter 
    Counter* m_counter; 
    // Shared pointer 
    T* m_ptr; 
}; 