// https://www.geeksforgeeks.org/how-to-implement-user-defined-shared-pointers-in-c/
#include "shared_ptr.h" 
using namespace std; 
  
// Class representing a reference counter class 
  
Counter::Counter() 
    : m_counter(0){};
void Counter::reset() { 
        m_counter = 0; 
} 
unsigned int Counter::get() { 
    return m_counter; 
} 
void Counter::operator++() { 
    m_counter++; 
} 
void Counter::operator++(int) { 
    m_counter++; 
} 
  
// Overload post/pre decrement 
void Counter::operator--() { 
    m_counter--; 
} 

void Counter::operator--(int) { 
    m_counter--; 
} 
  
// Overloading << operator 
ostream& operator<<(ostream& os, const Counter& counter) { 
    os << "Counter Value : " << counter.m_counter << endl; 
} 
  

// Class representing a shared pointer 
template <typename T>
Shared_ptr<T>::Shared_ptr(T* ptr = nullptr) { 
    m_ptr = ptr; 
    m_counter = new Counter(); 
    if (ptr) { 
        (*m_counter)++; 
    } 
} 
  
// Copy constructor
template <typename T>
Shared_ptr<T>::Shared_ptr(Shared_ptr<T>& sp) { 
    m_ptr = sp.m_ptr; 
    m_counter = sp.m_counter; 
    (*m_counter)++; 
} 

// Reference count
template <typename T>
unsigned int Shared_ptr<T>::use_count() { 
        return m_counter->get(); 
}     

template <typename T>
T* Shared_ptr<T>::get() { 
    return m_ptr; 
} 
  
template <typename T>
Shared_ptr<T>::~Shared_ptr() { 
    (*m_counter)--; 
    if (m_counter->get() == 0) { 
        delete m_counter; 
        delete m_ptr; 
    } 
} 

template <typename T>
ostream& operator<<(ostream& os, Shared_ptr<T>& sp) { 
    os << "Address pointed : "
        << sp.get() << endl; 
    cout << *(sp.m_counter) << endl; 
} 
  
