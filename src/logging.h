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
#include <iostream>

#ifdef VERBOSE
#define LOG std::cerr
#endif

#ifdef NVERBOSE
struct No_log {
  template <typename T> No_log const &operator<<(T const &) const {
    return *this;
  }
  No_log const &operator<<(std::ostream &(*)(std::ostream &)) const {
    return *this;
  }
};
#define LOG No_log()
#endif
