// TLDC "Too long; Didn't Count" A length limited path counter.
// Copyright (C) 2023-2024 Rafael Kiesel, Markus Hecher

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#pragma once

#include <assert.h>

#include "types.h"

namespace fpc {

template <typename count_t> struct Limited_count {
  Limited_count(Edge_length max_length, Edge_length offset = 0,
                std::vector<count_t> initial_count = std::vector<count_t>(1, 1))
      : offset_(offset), max_length_(max_length), counts_(initial_count) {}
  count_t total_count(Edge_length max) const {
    count_t rv = 0;
    for (Edge_length length = 0;
         length < counts_.size() && offset_ + length <= max; length++) {
      rv += counts_[length];
    }
    return rv;
  }
  Edge_length offset() const { return offset_; }
  void increment_offset() { ++offset_; }
  void decrement_offset() { --offset_; }
  Limited_count<count_t> &operator+=(Limited_count<count_t> &other) {
    assert(offset_ <= max_length_);
    assert(other.offset_ <= max_length_);

    if (offset_ > other.offset_) {
      // we reached this state with a smaller offset
      // add the paths we currently cached to the new result
      if (other.counts_.size() + other.offset_ < counts_.size() + offset_) {
        other.counts_.resize(
            std::min(counts_.size() + offset_ - other.offset_,
                     std::size_t(max_length_ - other.offset_ + 1)));
      }
      for (Edge_length length = 0;
           length < counts_.size() && length + offset_ <= max_length_;
           length++) {
        other.counts_[offset_ - other.offset_ + length] += counts_[length];
      }
      offset_ = other.offset_;
      counts_ = std::move(other.counts_);
    } else {
      if (counts_.size() + offset_ < other.counts_.size() + other.offset_) {
        counts_.resize(std::min(other.counts_.size() + other.offset_ - offset_,
                                std::size_t(max_length_ - offset_ + 1)));
      }
      for (Edge_length length = 0; length < other.counts_.size() &&
                                   length + other.offset_ <= max_length_;
           length++) {
        counts_[other.offset_ - offset_ + length] += other.counts_[length];
      }
    }
    return *this;
  }
  Limited_count<count_t> &operator*=(Limited_count<count_t> const &other) {
    assert(offset_ + other.offset_ <= max_length_);
    auto orig_size = counts_.size();
    auto new_size =
        std::min(counts_.size() + other.counts_.size() - 1,
                 std::size_t(max_length_ - (offset_ + other.offset_) + 1));
    counts_.resize(new_size);
    for (size_t sum_length =
             std::min(orig_size - 1 + other.counts_.size() - 1, new_size - 1);
         ; --sum_length) {
      size_t first = std::min(orig_size - 1, sum_length);
      size_t second = sum_length - first;
      count_t sum = 0;
      for (; second < other.counts_.size(); --first, ++second) {
        sum += counts_[first] * other.counts_[second];
        if (first == 0) {
          break;
        }
      }
      counts_[sum_length] = sum;
      if (sum_length == 0) {
        break;
      }
    }
    offset_ += other.offset_;
    return *this;
  }

  void add_one(Edge_length length) {
    assert(offset_ == 0);
    if (counts_.size() <= length) {
      counts_.resize(length + 1);
    }
    counts_[length] += 1;
  }

  Edge_length size() const { return offset_ + counts_.size() - 1; }

  count_t const &operator[](Edge_length at) const {
    assert(size() >= at);
    assert(at >= offset_);
    return counts_[at - offset_];
  }

  static Limited_count<count_t> zero(Edge_length max_length) {
    Limited_count<count_t> rv(max_length);
    rv.counts_ = {0};
    return rv;
  }

  Edge_length offset_ = 0;
  Edge_length max_length_ = 0;
  std::vector<count_t> counts_;
};

template <typename count_t> struct Unlimited_count {
  Unlimited_count(Edge_length, count_t count = 1)
      : offset_(0), count_(std::move(count)) {}
  count_t const &total_count(Edge_length) const { return count_; }
  Edge_length offset() const { return offset_; }
  void increment_offset() { ++offset_; }
  void decrement_offset() { --offset_; }
  Unlimited_count<count_t> &operator+=(Unlimited_count<count_t> &other) {
    count_ += other.count_;
    offset_ = std::max(offset_, other.offset_);
    return *this;
  }
  Unlimited_count<count_t> &operator*=(Unlimited_count<count_t> const &other) {
    count_ *= other.count_;
    offset_ += other.offset_;
    return *this;
  }

private:
  Edge_length offset_ = 0;
  count_t count_;
};

template <template <typename...> class, template <typename...> class>
struct is_same_template : std::false_type {};

template <template <typename...> class T>
struct is_same_template<T, T> : std::true_type {};

template <typename C, template <typename> typename S>
S<C> convert_to(Limited_count<mpz_class> const &original,
                Edge_length max_length) {
  if constexpr (is_same_template<S, Limited_count>::value) {
    std::stringstream ss;
    std::vector<C> values;
    for (auto const &val : original.counts_) {
      ss << val;
      C new_val;
      ss >> new_val;
      values.push_back(new_val);
      ss.str("");
      ss.clear();
    }
    return Limited_count<C>(max_length, original.offset(), values);
  } else {
    std::stringstream ss;
    C final_val = 0;
    for (auto const &val : original.counts_) {
      ss << val;
      C new_val;
      ss >> new_val;
      final_val += new_val;
      ss.str("");
      ss.clear();
    }
    return Unlimited_count<C>(0, final_val);
  }
}
} // namespace fpc