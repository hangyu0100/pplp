#include <cassert>
#include <cstdlib>
#include <atomic>
#include <stdexcept>

#ifndef VERIMAG_SIMPLE_CONCURRENT_VECTOR_MAX_SIZE
#define VERIMAG_SIMPLE_CONCURRENT_VECTOR_MAX_SIZE 10000
#endif

template <typename T> class simple_concurrent_vector {
 private:
  size_t max_size;
  std::atomic<size_t> filled;
  T* data;
  
 public:
  size_t size() const { return filled; }
      
  static const size_t DEFAULT_SIZE = VERIMAG_SIMPLE_CONCURRENT_VECTOR_MAX_SIZE;

  simple_concurrent_vector(size_t size0 = DEFAULT_SIZE):
  max_size(size0), filled(0) {
    data = (T*) malloc(sizeof(T) * size0);
    if (data==0) throw std::bad_alloc();
  }

  ~simple_concurrent_vector() {
    int fil = filled;
    for(int i=0; i<fil; i++) {
      data[i].~T();
    }
    free(data);
  }

  const T& operator[](size_t i) const {
    assert(i < filled);
    return data[i]; 
  }

  T& operator[](size_t i) {
    assert(i < filled);
    return data[i]; 
  }

  const T& at(size_t i) const {
    if (i >= filled) throw std::out_of_range("simple_concurrent_vector::at");
    return data[i];
  }

  T& at(size_t i) {
    if (i >= filled) throw std::out_of_range("simple_concurrent_vector::at");
    return data[i];
  }

  class iterator : public std::iterator < std::random_access_iterator_tag, T> {
    friend class simple_concurrent_vector<T>;
    
    simple_concurrent_vector<T>& vector;
    size_t index;

    iterator(simple_concurrent_vector<T>& vector0, size_t index0):
      vector(vector0), index(index0) {
    }

  public:
    T& operator*() const {
      assert(index < vector.filled);
      return vector.data[index];
    }

    T* operator->() const {
      assert(index < vector.filled);
      return vector.data+index;
    }

    void operator++(int) {
      index++;
    }

    void operator--(int) {
      index--;
    }

    void operator++() {
      index++;
    }

    void operator--() {
      index--;
    }

    bool operator!=(const iterator& b) const {
      return index != b.index;
    }

    bool operator==(const iterator& b) const {
      return index == b.index;
    }

    ptrdiff_t operator-(const iterator& b) const {
      return index - b.index;
    }
  };
  
  class const_iterator : public std::iterator < std::random_access_iterator_tag, T> {
    friend class simple_concurrent_vector<T>;
    
    simple_concurrent_vector<T>& vector;
    size_t index;

    const_iterator(simple_concurrent_vector<T>& vector0, size_t index0):
      vector(vector0), index(index0) {
    }

  public:
    const T& operator*() const {
      assert(index < vector.filled);
      return vector.data[index];
    }

    const T* operator->() const {
      assert(index < vector.filled);
      return vector.data+index;
    }

    void operator++(int) {
      index++;
    }

    void operator--(int) {
      if (index == 0) throw std::out_of_range();
      index--;
    }

    void operator++() {
      index++;
    }

    void operator--() {
      if (index == 0) throw std::out_of_range();
      index--;
    }

    bool operator!=(const const_iterator& b) const {
      return index != b.index;
    }

    bool operator==(const const_iterator& b) const {
      return index == b.index;
    }

    ptrdiff_t operator-(const const_iterator& b) const {
      return index - b.index;
    }
  };

  friend class iterator;
  friend class const_iterator;

  iterator begin() {
    return iterator(*this, 0);
  }

  iterator end() {
    return iterator(*this, filled);
  }

  const_iterator begin() const {
    return iterator(*this, 0);
  }

  const_iterator end() const {
    return iterator(*this, filled);
  }

  iterator push_back(T& v) {
    size_t pos = filled++;
    if (filled > max_size) throw std::bad_alloc();
    new(data+pos) T(v);
    return iterator(*this, pos);
  }

  iterator push_back(T&& v) {
    size_t pos = filled++;
    if (filled > max_size) throw std::bad_alloc();
    new(data+pos) T(std::move(v));
    return iterator(*this, pos);
  }

  template<typename... Args> iterator emplace_back(Args&&... args) {
    size_t pos = filled++;
    if (filled > max_size) throw std::bad_alloc();
    new(data+pos) T(std::forward<Args>(args)...);
    return iterator(*this, pos);
  }
};
