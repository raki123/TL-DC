#pragma once

// copied from https://github.com/Laakeri/sharpsat-td and modified a bit

#include <cstdint>
#include <cstdlib>

#define BITS 64

class Bitset {
 public:
  uint64_t* data_;
  size_t chunks_;
  Bitset() {
    chunks_ = 0;
    data_ = nullptr;
  }
  explicit Bitset(size_t size) {
    chunks_ = (size + BITS - 1) / BITS;
    data_ = (uint64_t*)std::malloc(chunks_*sizeof(uint64_t));
    for (size_t i=0;i<chunks_;i++){
      data_[i] = 0;
    }
  }
  ~Bitset() {
    std::free(data_);
  }
  Bitset(const Bitset& other) {
    chunks_ = other.chunks_;
    data_ = (uint64_t*)std::malloc(chunks_*sizeof(uint64_t));
    for (size_t i=0;i<chunks_;i++){
      data_[i] = other.data_[i];
    }
  }
  Bitset& operator=(const Bitset& other) {
    if (this != &other) {
      if (chunks_ != other.chunks_) {
        std::free(data_);
        chunks_ = other.chunks_;
        data_ = (uint64_t*)std::malloc(chunks_*sizeof(uint64_t));
      }
      for (size_t i=0;i<chunks_;i++){
        data_[i] = other.data_[i];
      }
    }
    return *this;
  }
  Bitset(Bitset&& other) {
    data_ = other.data_;
    chunks_ = other.chunks_;
    other.data_ = nullptr;
    other.chunks_ = 0;
  }
  Bitset& operator=(Bitset&& other) {
    if (this != &other) {
      std::free(data_);
      data_ = other.data_;
      chunks_ = other.chunks_;
      other.data_ = nullptr;
      other.chunks_ = 0;
    }
    return *this;
  }
  bool operator==(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i] != other.data_[i]) return false;
    }
    return true;
  }
  bool operator!=(const Bitset& other) const {
    for (size_t i=0;i<chunks_;i++){
      if (data_[i] != other.data_[i]) return true;
    }
    return false;
  }
  void Set(size_t i, bool v) {
    if (v) {
      data_[i/BITS] |= ((uint64_t)1 << (uint64_t)(i%BITS));
    } else {
      data_[i/BITS] &= (~((uint64_t)1 << (uint64_t)(i%BITS)));
    }
  }
  void SetTrue(size_t i) {
    data_[i/BITS] |= ((uint64_t)1 << (uint64_t)(i%BITS));
  }
  void SetFalse(size_t i) {
    data_[i/BITS] &= (~((uint64_t)1 << (uint64_t)(i%BITS)));
  }
  bool Get(size_t i) const {
    return data_[i/BITS] & ((uint64_t)1 << (uint64_t)(i%BITS));
  }
  void Clear() {
    for (size_t i=0;i<chunks_;i++){
      data_[i] = 0;
    }
  }
  
};