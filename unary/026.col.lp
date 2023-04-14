reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17), not reach(X, 18), not reach(X, 19), not reach(X, 20), not reach(X, 21), not reach(X, 22), not reach(X, 23), not reach(X, 24), not reach(X, 25), not reach(X, 26), not reach(X, 27), not reach(X, 28), not reach(X, 29), not reach(X, 30), not reach(X, 31), not reach(X, 32), not reach(X, 33), not reach(X, 34), not reach(X, 35), not reach(X, 36), not reach(X, 37), not reach(X, 38), not reach(X, 39), not reach(X, 40), not reach(X, 41), not reach(X, 42), not reach(X, 43), not reach(X, 44), not reach(X, 45).
start(3).
goal(28).
:- reach(X, L), reach(X, L'), L != L'.
reach(2,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1), not reach(44, L+1).
reach(4,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1), not reach(44, L+1).
reach(5,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1), not reach(44, L+1).
reach(6,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(8, L+1), not reach(44, L+1).
reach(7,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(8, L+1), not reach(44, L+1).
reach(8,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(1, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1).
reach(1,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(10, L+1).
reach(4,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(10, L+1).
reach(5,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(10, L+1).
reach(6,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(10, L+1).
reach(7,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1), not reach(10, L+1).
reach(9,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1), not reach(10, L+1).
reach(8,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(10, L+1).
reach(10,L + 1) :- reach(2, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(1,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(2, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(2,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(4,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(2, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(5,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(2, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(6,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(2, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(7,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(2, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1).
reach(9,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(2, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(3, L), L <= 40, L>= 0, not reach(1, L+1), not reach(2, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1).
reach(1,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(2,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(5,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(6,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(7,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1).
reach(9,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(4, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1).
reach(1,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(2,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(4,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(6,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(7,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(9,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(8,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(11, L+1), not reach(16, L+1), not reach(41, L+1).
reach(11,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(16, L+1), not reach(41, L+1).
reach(16,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(5, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1), not reach(11, L+1), not reach(16, L+1).
reach(1,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(2,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(4,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(5,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(7, L+1), not reach(9, L+1), not reach(8, L+1).
reach(7,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(9, L+1), not reach(8, L+1).
reach(9,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(6, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1).
reach(1,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1), not reach(12, L+1).
reach(2,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1), not reach(12, L+1).
reach(4,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1), not reach(12, L+1).
reach(5,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1), not reach(12, L+1).
reach(6,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(9, L+1), not reach(8, L+1), not reach(12, L+1).
reach(9,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(8, L+1), not reach(12, L+1).
reach(8,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(12, L+1).
reach(12,L + 1) :- reach(7, L), L <= 41, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1), not reach(8, L+1).
reach(1,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1).
reach(2,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(1, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1).
reach(4,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1).
reach(5,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(9, L+1).
reach(6,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(9, L+1).
reach(7,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(9, L+1).
reach(9,L + 1) :- reach(8, L), L <= 40, L>= 1, not reach(1, L+1), not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1).
reach(1,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(12, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(12,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(37,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(38,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(37, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(40,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(42,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(43, L+1), not reach(45, L+1).
reach(43,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(44, L), L <= 42, L>= 2, not reach(1, L+1), not reach(12, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1).
reach(2,L + 1) :- reach(9, L), L <= 40, L>= 1, not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1).
reach(4,L + 1) :- reach(9, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1).
reach(5,L + 1) :- reach(9, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(6, L+1), not reach(7, L+1), not reach(8, L+1).
reach(6,L + 1) :- reach(9, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(7, L+1), not reach(8, L+1).
reach(7,L + 1) :- reach(9, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(9, L), L <= 40, L>= 1, not reach(2, L+1), not reach(3, L+1), not reach(4, L+1), not reach(5, L+1), not reach(6, L+1), not reach(7, L+1).
reach(2,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(11, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(27, L+1).
reach(11,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(27, L+1).
reach(13,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(27, L+1).
reach(12,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(27, L+1).
reach(14,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(13, L+1), not reach(12, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(27, L+1).
reach(16,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(17, L+1), not reach(18, L+1), not reach(27, L+1).
reach(17,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(18, L+1), not reach(27, L+1).
reach(18,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(10, L), L <= 42, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(5,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(10, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(10,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(13,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(12,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(14,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(13, L+1), not reach(12, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(15,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(16,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1).
reach(17,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(18, L+1).
reach(18,L + 1) :- reach(11, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(13, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1).
reach(5,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(18, L+1), not reach(20, L+1).
reach(10,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(18, L+1), not reach(20, L+1).
reach(11,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(18, L+1), not reach(20, L+1).
reach(12,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(15, L+1), not reach(18, L+1), not reach(20, L+1).
reach(13,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(15, L+1), not reach(18, L+1), not reach(20, L+1).
reach(15,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(18, L+1), not reach(20, L+1).
reach(18,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(20, L+1).
reach(20,L + 1) :- reach(16, L), L <= 42, L>= 2, not reach(5, L+1), not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(18, L+1).
reach(5,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(37,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(38,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(37, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(40,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(42,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(43, L+1), not reach(45, L+1).
reach(43,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(41, L), L <= 42, L>= 2, not reach(5, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1).
reach(7,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(10,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(11,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(13,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(14,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(15,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(16,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(17,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(18, L+1), not reach(42, L+1), not reach(44, L+1).
reach(18,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(42, L+1), not reach(44, L+1).
reach(42,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(12, L), L <= 42, L>= 2, not reach(7, L+1), not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(42, L+1).
reach(10,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(36, L+1).
reach(11,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(36, L+1).
reach(12,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(36, L+1).
reach(14,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(36, L+1).
reach(15,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(36, L+1).
reach(16,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1), not reach(36, L+1).
reach(17,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(18, L+1), not reach(36, L+1).
reach(18,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(13, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1).
reach(10,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1), not reach(19, L+1).
reach(11,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1), not reach(19, L+1).
reach(12,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1), not reach(19, L+1).
reach(13,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1), not reach(19, L+1).
reach(15,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(17, L+1), not reach(18, L+1), not reach(19, L+1).
reach(17,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(18, L+1), not reach(19, L+1).
reach(18,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(17, L+1), not reach(19, L+1).
reach(19,L + 1) :- reach(14, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(17, L+1), not reach(18, L+1).
reach(10,L + 1) :- reach(17, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1).
reach(11,L + 1) :- reach(17, L), L <= 42, L>= 3, not reach(10, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1).
reach(12,L + 1) :- reach(17, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1).
reach(13,L + 1) :- reach(17, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1).
reach(14,L + 1) :- reach(17, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1).
reach(15,L + 1) :- reach(17, L), L <= 42, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1).
reach(10,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(32, L+1).
reach(11,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(32, L+1).
reach(12,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(32, L+1).
reach(13,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1), not reach(32, L+1).
reach(14,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(15, L+1), not reach(16, L+1), not reach(32, L+1).
reach(15,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(32, L+1).
reach(16,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(32, L+1).
reach(32,L + 1) :- reach(18, L), L <= 43, L>= 3, not reach(10, L+1), not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(15, L+1), not reach(16, L+1).
reach(10,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1).
reach(19,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1).
reach(20,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1).
reach(21,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1).
reach(22,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1).
reach(23,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1).
reach(24,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1).
reach(25,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(27, L), L <= 43, L>= 3, not reach(10, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1).
reach(11,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(20, L+1), not reach(26, L+1).
reach(12,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(20, L+1), not reach(26, L+1).
reach(13,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(20, L+1), not reach(26, L+1).
reach(14,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(20, L+1), not reach(26, L+1).
reach(16,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(17, L+1), not reach(18, L+1), not reach(20, L+1), not reach(26, L+1).
reach(17,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(18, L+1), not reach(20, L+1), not reach(26, L+1).
reach(18,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(20, L+1), not reach(26, L+1).
reach(20,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(15, L), L <= 42, L>= 3, not reach(11, L+1), not reach(12, L+1), not reach(13, L+1), not reach(14, L+1), not reach(16, L+1), not reach(17, L+1), not reach(18, L+1), not reach(20, L+1).
reach(12,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(38, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(38,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(12, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(12, L+1), not reach(38, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(41,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(12, L+1), not reach(38, L+1), not reach(39, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(43,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(12, L+1), not reach(38, L+1), not reach(39, L+1), not reach(41, L+1), not reach(44, L+1), not reach(45, L+1).
reach(44,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(12, L+1), not reach(38, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(42, L), L <= 42, L>= 3, not reach(12, L+1), not reach(38, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1).
reach(13,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1).
reach(28,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1).
reach(29,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1).
reach(30,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1).
reach(31,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1).
reach(32,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(35, L+1).
reach(33,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(35, L+1).
reach(35,L + 1) :- reach(36, L), L <= 44, L>= 4, not reach(13, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1).
reach(14,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(20, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(30, L+1).
reach(21,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(22, L+1), not reach(24, L+1), not reach(20, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(30, L+1).
reach(22,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(24, L+1), not reach(20, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(30, L+1).
reach(24,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(22, L+1), not reach(20, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(30, L+1).
reach(20,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(30, L+1).
reach(25,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(20, L+1), not reach(26, L+1), not reach(27, L+1), not reach(30, L+1).
reach(26,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(20, L+1), not reach(25, L+1), not reach(27, L+1), not reach(30, L+1).
reach(27,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(20, L+1), not reach(25, L+1), not reach(26, L+1), not reach(30, L+1).
reach(30,L + 1) :- reach(19, L), L <= 43, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(20, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(15,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(16,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(19,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(21,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(22,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(23,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(24,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(25,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(26,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1), not reach(29, L+1), not reach(35, L+1).
reach(27,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(29, L+1), not reach(35, L+1).
reach(29,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(35, L+1).
reach(35,L + 1) :- reach(20, L), L <= 43, L>= 3, not reach(15, L+1), not reach(16, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(29, L+1).
reach(15,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1).
reach(19,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1).
reach(20,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1).
reach(21,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1).
reach(22,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1).
reach(24,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(25, L+1), not reach(27, L+1).
reach(25,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(26, L), L <= 43, L>= 4, not reach(15, L+1), not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1).
reach(18,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(28,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(30,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(28, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(31,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(28, L+1), not reach(30, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(33,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(34,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(35, L+1), not reach(36, L+1).
reach(35,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(32, L), L <= 44, L>= 4, not reach(18, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1).
reach(19,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(20,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(22,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(23,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(24,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(25,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(26,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1), not reach(33, L+1).
reach(27,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(33, L+1).
reach(33,L + 1) :- reach(21, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(19,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(28, L+1).
reach(20,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(28, L+1).
reach(21,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(28, L+1).
reach(23,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(28, L+1).
reach(24,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1), not reach(28, L+1).
reach(25,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1), not reach(28, L+1).
reach(26,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1), not reach(28, L+1).
reach(27,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(28, L+1).
reach(28,L + 1) :- reach(22, L), L <= 44, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(19,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(20,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(21,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(22,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(23,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(25, L+1), not reach(26, L+1), not reach(27, L+1).
reach(25,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(26, L+1), not reach(27, L+1).
reach(26,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(24, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1).
reach(19,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1).
reach(20,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1).
reach(21,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1).
reach(22,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1).
reach(23,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(26, L+1), not reach(27, L+1).
reach(24,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(26, L+1), not reach(27, L+1).
reach(26,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(25, L), L <= 43, L>= 4, not reach(19, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(23, L+1), not reach(24, L+1), not reach(26, L+1).
reach(19,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(28,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(29,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(31,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(29, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(32,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(33,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(34,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1), not reach(36, L+1).
reach(35,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(30, L), L <= 44, L>= 5, not reach(19, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1).
reach(20,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1), not reach(28, L+1).
reach(21,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(20, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1), not reach(28, L+1).
reach(22,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1), not reach(28, L+1).
reach(24,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(25, L+1), not reach(27, L+1), not reach(28, L+1).
reach(25,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(27, L+1), not reach(28, L+1).
reach(27,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(28, L+1).
reach(28,L + 1) :- reach(23, L), L <= 44, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(24, L+1), not reach(25, L+1), not reach(27, L+1).
reach(20,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1), not reach(39, L+1).
reach(28,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1), not reach(39, L+1).
reach(30,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1), not reach(39, L+1).
reach(31,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(30, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1), not reach(39, L+1).
reach(33,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1), not reach(39, L+1).
reach(34,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(35, L+1), not reach(36, L+1), not reach(39, L+1).
reach(35,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(36, L+1), not reach(39, L+1).
reach(36,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(29, L), L <= 44, L>= 4, not reach(20, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(20,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(43, L+1).
reach(29,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(43, L+1).
reach(30,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(43, L+1).
reach(31,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(43, L+1).
reach(32,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(34, L+1), not reach(36, L+1), not reach(43, L+1).
reach(34,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(36, L+1), not reach(43, L+1).
reach(36,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(35, L), L <= 43, L>= 4, not reach(20, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1).
reach(21,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(37, L+1).
reach(28,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(37, L+1).
reach(29,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(37, L+1).
reach(30,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(37, L+1).
reach(31,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1), not reach(37, L+1).
reach(32,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(34, L+1), not reach(36, L+1), not reach(37, L+1).
reach(34,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(36, L+1), not reach(37, L+1).
reach(36,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(33, L), L <= 44, L>= 4, not reach(21, L+1), not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(34, L+1), not reach(36, L+1).
reach(28,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(29,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(30,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(32,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(33,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(34, L+1), not reach(35, L+1), not reach(36, L+1).
reach(34,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1), not reach(36, L+1).
reach(35,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(31, L), L <= 44, L>= 5, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(34, L+1), not reach(35, L+1).
reach(28,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1), not reach(38, L+1).
reach(29,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1), not reach(38, L+1).
reach(30,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1), not reach(38, L+1).
reach(31,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1), not reach(38, L+1).
reach(32,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(33, L+1), not reach(35, L+1), not reach(38, L+1).
reach(33,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(35, L+1), not reach(38, L+1).
reach(35,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(34, L), L <= 44, L>= 4, not reach(28, L+1), not reach(29, L+1), not reach(30, L+1), not reach(31, L+1), not reach(32, L+1), not reach(33, L+1), not reach(35, L+1).
reach(29,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(37,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(38,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(40,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(38, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(41,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(42,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(43,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1), not reach(45, L+1).
reach(44,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(39, L), L <= 43, L>= 3, not reach(29, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1).
reach(33,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(33, L+1), not reach(40, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(40,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(33, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(41,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(33, L+1), not reach(39, L+1), not reach(40, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(43,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(33, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(44, L+1), not reach(45, L+1).
reach(44,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(33, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(43, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(37, L), L <= 43, L>= 3, not reach(33, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1).
reach(34,L + 1) :- reach(38, L), L <= 43, L>= 3, not reach(39, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1).
reach(39,L + 1) :- reach(38, L), L <= 43, L>= 3, not reach(34, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1).
reach(41,L + 1) :- reach(38, L), L <= 43, L>= 3, not reach(34, L+1), not reach(39, L+1), not reach(42, L+1), not reach(43, L+1), not reach(44, L+1).
reach(42,L + 1) :- reach(38, L), L <= 43, L>= 3, not reach(34, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1).
reach(43,L + 1) :- reach(38, L), L <= 43, L>= 3, not reach(34, L+1), not reach(39, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(38, L), L <= 43, L>= 3, not reach(34, L+1), not reach(39, L+1), not reach(41, L+1), not reach(42, L+1), not reach(43, L+1).
reach(35,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(37,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(38,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(37, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(39,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(40,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(41,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(44, L+1).
reach(42,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(43, L), L <= 42, L>= 3, not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1).
reach(37,L + 1) :- reach(40, L), L <= 42, L>= 3, not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(40, L), L <= 42, L>= 3, not reach(37, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(41,L + 1) :- reach(40, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(43, L+1), not reach(44, L+1), not reach(45, L+1).
reach(43,L + 1) :- reach(40, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(41, L+1), not reach(44, L+1), not reach(45, L+1).
reach(44,L + 1) :- reach(40, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(40, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(41, L+1), not reach(43, L+1), not reach(44, L+1).
reach(37,L + 1) :- reach(45, L), L <= 42, L>= 3, not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(39,L + 1) :- reach(45, L), L <= 42, L>= 3, not reach(37, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(40,L + 1) :- reach(45, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(41, L+1), not reach(42, L+1), not reach(44, L+1).
reach(41,L + 1) :- reach(45, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(40, L+1), not reach(42, L+1), not reach(44, L+1).
reach(42,L + 1) :- reach(45, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(45, L), L <= 42, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(40, L+1), not reach(41, L+1), not reach(42, L+1).

