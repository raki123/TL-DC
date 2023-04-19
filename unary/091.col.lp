reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14).
start(18).
goal(61).
:- reach(X, L), reach(X, L'), L != L'.
reach(13,L + 1) :- reach(1, L), L <= 11, L>= 4, not reach(14, L+1), not reach(79, L+1), not reach(37, L+1).
reach(14,L + 1) :- reach(1, L), L <= 11, L>= 4, not reach(13, L+1), not reach(79, L+1), not reach(37, L+1).
reach(79,L + 1) :- reach(1, L), L <= 11, L>= 4, not reach(13, L+1), not reach(14, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(1, L), L <= 11, L>= 4, not reach(13, L+1), not reach(14, L+1), not reach(79, L+1).
reach(1,L + 1) :- reach(13, L), L <= 12, L>= 5, not reach(2, L+1), not reach(66, L+1).
reach(2,L + 1) :- reach(13, L), L <= 12, L>= 5, not reach(1, L+1), not reach(66, L+1).
reach(66,L + 1) :- reach(13, L), L <= 12, L>= 5, not reach(1, L+1), not reach(2, L+1).
reach(1,L + 1) :- reach(14, L), L <= 11, L>= 5, not reach(52, L+1), not reach(82, L+1), not reach(29, L+1).
reach(52,L + 1) :- reach(14, L), L <= 11, L>= 5, not reach(1, L+1), not reach(82, L+1), not reach(29, L+1).
reach(82,L + 1) :- reach(14, L), L <= 11, L>= 5, not reach(1, L+1), not reach(52, L+1), not reach(29, L+1).
reach(29,L + 1) :- reach(14, L), L <= 11, L>= 5, not reach(1, L+1), not reach(52, L+1), not reach(82, L+1).
reach(1,L + 1) :- reach(79, L), L <= 11, L>= 3, not reach(39, L+1), not reach(73, L+1), not reach(86, L+1).
reach(39,L + 1) :- reach(79, L), L <= 11, L>= 3, not reach(1, L+1), not reach(73, L+1), not reach(86, L+1).
reach(73,L + 1) :- reach(79, L), L <= 11, L>= 3, not reach(1, L+1), not reach(39, L+1), not reach(86, L+1).
reach(86,L + 1) :- reach(79, L), L <= 11, L>= 3, not reach(1, L+1), not reach(39, L+1), not reach(73, L+1).
reach(13,L + 1) :- reach(2, L), L <= 11, L>= 5, not reach(21, L+1), not reach(59, L+1).
reach(21,L + 1) :- reach(2, L), L <= 11, L>= 5, not reach(13, L+1), not reach(59, L+1).
reach(59,L + 1) :- reach(2, L), L <= 11, L>= 5, not reach(13, L+1), not reach(21, L+1).
reach(2,L + 1) :- reach(21, L), L <= 11, L>= 4, not reach(82, L+1), not reach(41, L+1), not reach(23, L+1).
reach(82,L + 1) :- reach(21, L), L <= 11, L>= 4, not reach(2, L+1), not reach(41, L+1), not reach(23, L+1).
reach(41,L + 1) :- reach(21, L), L <= 11, L>= 4, not reach(2, L+1), not reach(82, L+1), not reach(23, L+1).
reach(23,L + 1) :- reach(21, L), L <= 11, L>= 4, not reach(2, L+1), not reach(82, L+1), not reach(41, L+1).
reach(2,L + 1) :- reach(59, L), L <= 12, L>= 5, not reach(27, L+1), not reach(33, L+1), not reach(71, L+1).
reach(27,L + 1) :- reach(59, L), L <= 12, L>= 5, not reach(2, L+1), not reach(33, L+1), not reach(71, L+1).
reach(33,L + 1) :- reach(59, L), L <= 12, L>= 5, not reach(2, L+1), not reach(27, L+1), not reach(71, L+1).
reach(71,L + 1) :- reach(59, L), L <= 12, L>= 5, not reach(2, L+1), not reach(27, L+1), not reach(33, L+1).
reach(39,L + 1) :- reach(3, L), L <= 12, L>= 4, not reach(40, L+1), not reach(73, L+1), not reach(19, L+1).
reach(40,L + 1) :- reach(3, L), L <= 12, L>= 4, not reach(39, L+1), not reach(73, L+1), not reach(19, L+1).
reach(73,L + 1) :- reach(3, L), L <= 12, L>= 4, not reach(39, L+1), not reach(40, L+1), not reach(19, L+1).
reach(19,L + 1) :- reach(3, L), L <= 12, L>= 4, not reach(39, L+1), not reach(40, L+1), not reach(73, L+1).
reach(3,L + 1) :- reach(39, L), L <= 11, L>= 4, not reach(79, L+1), not reach(42, L+1).
reach(79,L + 1) :- reach(39, L), L <= 11, L>= 4, not reach(3, L+1), not reach(42, L+1).
reach(42,L + 1) :- reach(39, L), L <= 11, L>= 4, not reach(3, L+1), not reach(79, L+1).
reach(3,L + 1) :- reach(40, L), L <= 13, L>= 3, not reach(25, L+1), not reach(82, L+1), not reach(61, L+1), not reach(86, L+1).
reach(25,L + 1) :- reach(40, L), L <= 13, L>= 3, not reach(3, L+1), not reach(82, L+1), not reach(61, L+1), not reach(86, L+1).
reach(82,L + 1) :- reach(40, L), L <= 13, L>= 3, not reach(3, L+1), not reach(25, L+1), not reach(61, L+1), not reach(86, L+1).
reach(61,L + 1) :- reach(40, L), L <= 13, L>= 3, not reach(3, L+1), not reach(25, L+1), not reach(82, L+1), not reach(86, L+1).
reach(86,L + 1) :- reach(40, L), L <= 13, L>= 3, not reach(3, L+1), not reach(25, L+1), not reach(82, L+1), not reach(61, L+1).
reach(3,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(52, L+1), not reach(79, L+1), not reach(38, L+1).
reach(23,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(3, L+1), not reach(25, L+1), not reach(26, L+1), not reach(52, L+1), not reach(79, L+1), not reach(38, L+1).
reach(25,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(3, L+1), not reach(23, L+1), not reach(26, L+1), not reach(52, L+1), not reach(79, L+1), not reach(38, L+1).
reach(26,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(3, L+1), not reach(23, L+1), not reach(25, L+1), not reach(52, L+1), not reach(79, L+1), not reach(38, L+1).
reach(52,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(3, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(79, L+1), not reach(38, L+1).
reach(79,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(3, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(52, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(73, L), L <= 11, L>= 4, not reach(3, L+1), not reach(23, L+1), not reach(25, L+1), not reach(26, L+1), not reach(52, L+1), not reach(79, L+1).
reach(26,L + 1) :- reach(4, L), L <= 10, L>= 6, not reach(52, L+1), not reach(19, L+1).
reach(52,L + 1) :- reach(4, L), L <= 10, L>= 6, not reach(26, L+1), not reach(19, L+1).
reach(19,L + 1) :- reach(4, L), L <= 10, L>= 6, not reach(26, L+1), not reach(52, L+1).
reach(4,L + 1) :- reach(26, L), L <= 10, L>= 5, not reach(65, L+1), not reach(73, L+1), not reach(50, L+1).
reach(65,L + 1) :- reach(26, L), L <= 10, L>= 5, not reach(4, L+1), not reach(73, L+1), not reach(50, L+1).
reach(73,L + 1) :- reach(26, L), L <= 10, L>= 5, not reach(4, L+1), not reach(65, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(26, L), L <= 10, L>= 5, not reach(4, L+1), not reach(65, L+1), not reach(73, L+1).
reach(4,L + 1) :- reach(52, L), L <= 10, L>= 5, not reach(14, L+1), not reach(73, L+1), not reach(50, L+1).
reach(14,L + 1) :- reach(52, L), L <= 10, L>= 5, not reach(4, L+1), not reach(73, L+1), not reach(50, L+1).
reach(73,L + 1) :- reach(52, L), L <= 10, L>= 5, not reach(4, L+1), not reach(14, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(52, L), L <= 10, L>= 5, not reach(4, L+1), not reach(14, L+1), not reach(73, L+1).
reach(4,L + 1) :- reach(19, L), L <= 11, L>= 5, not reach(6, L+1), not reach(3, L+1).
reach(6,L + 1) :- reach(19, L), L <= 11, L>= 5, not reach(4, L+1), not reach(3, L+1).
reach(3,L + 1) :- reach(19, L), L <= 11, L>= 5, not reach(4, L+1), not reach(6, L+1).
reach(36,L + 1) :- reach(5, L), L <= 11, L>= 1, not reach(16, L+1), not reach(86, L+1), not reach(18, L+1).
reach(16,L + 1) :- reach(5, L), L <= 11, L>= 1, not reach(36, L+1), not reach(86, L+1), not reach(18, L+1).
reach(86,L + 1) :- reach(5, L), L <= 11, L>= 1, not reach(36, L+1), not reach(16, L+1), not reach(18, L+1).
reach(5,L + 1) :- reach(36, L), L <= 10, L>= 2, not reach(43, L+1), not reach(41, L+1), not reach(34, L+1).
reach(43,L + 1) :- reach(36, L), L <= 10, L>= 2, not reach(5, L+1), not reach(41, L+1), not reach(34, L+1).
reach(41,L + 1) :- reach(36, L), L <= 10, L>= 2, not reach(5, L+1), not reach(43, L+1), not reach(34, L+1).
reach(34,L + 1) :- reach(36, L), L <= 10, L>= 2, not reach(5, L+1), not reach(43, L+1), not reach(41, L+1).
reach(5,L + 1) :- reach(16, L), L <= 10, L>= 2, not reach(7, L+1), not reach(8, L+1), not reach(46, L+1).
reach(7,L + 1) :- reach(16, L), L <= 10, L>= 2, not reach(5, L+1), not reach(8, L+1), not reach(46, L+1).
reach(8,L + 1) :- reach(16, L), L <= 10, L>= 2, not reach(5, L+1), not reach(7, L+1), not reach(46, L+1).
reach(46,L + 1) :- reach(16, L), L <= 10, L>= 2, not reach(5, L+1), not reach(7, L+1), not reach(8, L+1).
reach(5,L + 1) :- reach(86, L), L <= 12, L>= 2, not reach(79, L+1), not reach(40, L+1).
reach(79,L + 1) :- reach(86, L), L <= 12, L>= 2, not reach(5, L+1), not reach(40, L+1).
reach(40,L + 1) :- reach(86, L), L <= 12, L>= 2, not reach(5, L+1), not reach(79, L+1).
reach(19,L + 1) :- reach(6, L), L <= 10, L>= 4, not reach(46, L+1), not reach(71, L+1).
reach(46,L + 1) :- reach(6, L), L <= 10, L>= 4, not reach(19, L+1), not reach(71, L+1).
reach(71,L + 1) :- reach(6, L), L <= 10, L>= 4, not reach(19, L+1), not reach(46, L+1).
reach(6,L + 1) :- reach(46, L), L <= 10, L>= 3, not reach(8, L+1), not reach(16, L+1), not reach(43, L+1), not reach(34, L+1).
reach(8,L + 1) :- reach(46, L), L <= 10, L>= 3, not reach(6, L+1), not reach(16, L+1), not reach(43, L+1), not reach(34, L+1).
reach(16,L + 1) :- reach(46, L), L <= 10, L>= 3, not reach(6, L+1), not reach(8, L+1), not reach(43, L+1), not reach(34, L+1).
reach(43,L + 1) :- reach(46, L), L <= 10, L>= 3, not reach(6, L+1), not reach(8, L+1), not reach(16, L+1), not reach(34, L+1).
reach(34,L + 1) :- reach(46, L), L <= 10, L>= 3, not reach(6, L+1), not reach(8, L+1), not reach(16, L+1), not reach(43, L+1).
reach(16,L + 1) :- reach(7, L), L <= 11, L>= 3, not reach(25, L+1), not reach(37, L+1).
reach(25,L + 1) :- reach(7, L), L <= 11, L>= 3, not reach(16, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(7, L), L <= 11, L>= 3, not reach(16, L+1), not reach(25, L+1).
reach(7,L + 1) :- reach(25, L), L <= 12, L>= 4, not reach(37, L+1), not reach(40, L+1), not reach(73, L+1).
reach(37,L + 1) :- reach(25, L), L <= 12, L>= 4, not reach(7, L+1), not reach(40, L+1), not reach(73, L+1).
reach(40,L + 1) :- reach(25, L), L <= 12, L>= 4, not reach(7, L+1), not reach(37, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(25, L), L <= 12, L>= 4, not reach(7, L+1), not reach(37, L+1), not reach(40, L+1).
reach(16,L + 1) :- reach(8, L), L <= 10, L>= 3, not reach(46, L+1), not reach(47, L+1).
reach(46,L + 1) :- reach(8, L), L <= 10, L>= 3, not reach(16, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(8, L), L <= 10, L>= 3, not reach(16, L+1), not reach(46, L+1).
reach(73,L + 1) :- reach(23, L), L <= 11, L>= 5, not reach(80, L+1), not reach(21, L+1).
reach(80,L + 1) :- reach(23, L), L <= 11, L>= 5, not reach(73, L+1), not reach(21, L+1).
reach(21,L + 1) :- reach(23, L), L <= 11, L>= 5, not reach(73, L+1), not reach(80, L+1).
reach(71,L + 1) :- reach(10, L), L <= 12, L>= 4, not reach(33, L+1), not reach(34, L+1).
reach(33,L + 1) :- reach(10, L), L <= 12, L>= 4, not reach(71, L+1), not reach(34, L+1).
reach(34,L + 1) :- reach(10, L), L <= 12, L>= 4, not reach(71, L+1), not reach(33, L+1).
reach(10,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(29, L+1), not reach(34, L+1), not reach(41, L+1), not reach(59, L+1), not reach(6, L+1), not reach(27, L+1).
reach(29,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(10, L+1), not reach(34, L+1), not reach(41, L+1), not reach(59, L+1), not reach(6, L+1), not reach(27, L+1).
reach(34,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(10, L+1), not reach(29, L+1), not reach(41, L+1), not reach(59, L+1), not reach(6, L+1), not reach(27, L+1).
reach(41,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(10, L+1), not reach(29, L+1), not reach(34, L+1), not reach(59, L+1), not reach(6, L+1), not reach(27, L+1).
reach(59,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(10, L+1), not reach(29, L+1), not reach(34, L+1), not reach(41, L+1), not reach(6, L+1), not reach(27, L+1).
reach(6,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(10, L+1), not reach(29, L+1), not reach(34, L+1), not reach(41, L+1), not reach(59, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(71, L), L <= 11, L>= 4, not reach(10, L+1), not reach(29, L+1), not reach(34, L+1), not reach(41, L+1), not reach(59, L+1), not reach(6, L+1).
reach(27,L + 1) :- reach(41, L), L <= 11, L>= 3, not reach(36, L+1), not reach(38, L+1), not reach(71, L+1), not reach(42, L+1), not reach(21, L+1).
reach(36,L + 1) :- reach(41, L), L <= 11, L>= 3, not reach(27, L+1), not reach(38, L+1), not reach(71, L+1), not reach(42, L+1), not reach(21, L+1).
reach(38,L + 1) :- reach(41, L), L <= 11, L>= 3, not reach(27, L+1), not reach(36, L+1), not reach(71, L+1), not reach(42, L+1), not reach(21, L+1).
reach(71,L + 1) :- reach(41, L), L <= 11, L>= 3, not reach(27, L+1), not reach(36, L+1), not reach(38, L+1), not reach(42, L+1), not reach(21, L+1).
reach(42,L + 1) :- reach(41, L), L <= 11, L>= 3, not reach(27, L+1), not reach(36, L+1), not reach(38, L+1), not reach(71, L+1), not reach(21, L+1).
reach(21,L + 1) :- reach(41, L), L <= 11, L>= 3, not reach(27, L+1), not reach(36, L+1), not reach(38, L+1), not reach(71, L+1), not reach(42, L+1).
reach(65,L + 1) :- reach(12, L), L <= 11, L>= 5, not reach(42, L+1), not reach(80, L+1).
reach(42,L + 1) :- reach(12, L), L <= 11, L>= 5, not reach(65, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(12, L), L <= 11, L>= 5, not reach(65, L+1), not reach(42, L+1).
reach(12,L + 1) :- reach(65, L), L <= 11, L>= 5, not reach(26, L+1), not reach(82, L+1).
reach(26,L + 1) :- reach(65, L), L <= 11, L>= 5, not reach(12, L+1), not reach(82, L+1).
reach(82,L + 1) :- reach(65, L), L <= 11, L>= 5, not reach(12, L+1), not reach(26, L+1).
reach(12,L + 1) :- reach(42, L), L <= 11, L>= 4, not reach(27, L+1), not reach(41, L+1), not reach(43, L+1), not reach(39, L+1).
reach(27,L + 1) :- reach(42, L), L <= 11, L>= 4, not reach(12, L+1), not reach(41, L+1), not reach(43, L+1), not reach(39, L+1).
reach(41,L + 1) :- reach(42, L), L <= 11, L>= 4, not reach(12, L+1), not reach(27, L+1), not reach(43, L+1), not reach(39, L+1).
reach(43,L + 1) :- reach(42, L), L <= 11, L>= 4, not reach(12, L+1), not reach(27, L+1), not reach(41, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(42, L), L <= 11, L>= 4, not reach(12, L+1), not reach(27, L+1), not reach(41, L+1), not reach(43, L+1).
reach(12,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(23, L+1), not reach(33, L+1), not reach(34, L+1), not reach(37, L+1), not reach(38, L+1), not reach(47, L+1).
reach(23,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(12, L+1), not reach(33, L+1), not reach(34, L+1), not reach(37, L+1), not reach(38, L+1), not reach(47, L+1).
reach(33,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(12, L+1), not reach(23, L+1), not reach(34, L+1), not reach(37, L+1), not reach(38, L+1), not reach(47, L+1).
reach(34,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(12, L+1), not reach(23, L+1), not reach(33, L+1), not reach(37, L+1), not reach(38, L+1), not reach(47, L+1).
reach(37,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(12, L+1), not reach(23, L+1), not reach(33, L+1), not reach(34, L+1), not reach(38, L+1), not reach(47, L+1).
reach(38,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(12, L+1), not reach(23, L+1), not reach(33, L+1), not reach(34, L+1), not reach(37, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(80, L), L <= 12, L>= 4, not reach(12, L+1), not reach(23, L+1), not reach(33, L+1), not reach(34, L+1), not reach(37, L+1), not reach(38, L+1).
reach(13,L + 1) :- reach(66, L), L <= 13, L>= 5, not reach(31, L+1), not reach(61, L+1), not reach(27, L+1).
reach(31,L + 1) :- reach(66, L), L <= 13, L>= 5, not reach(13, L+1), not reach(61, L+1), not reach(27, L+1).
reach(61,L + 1) :- reach(66, L), L <= 13, L>= 5, not reach(13, L+1), not reach(31, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(66, L), L <= 13, L>= 5, not reach(13, L+1), not reach(31, L+1), not reach(61, L+1).
reach(14,L + 1) :- reach(82, L), L <= 12, L>= 4, not reach(21, L+1), not reach(40, L+1), not reach(65, L+1), not reach(38, L+1).
reach(21,L + 1) :- reach(82, L), L <= 12, L>= 4, not reach(14, L+1), not reach(40, L+1), not reach(65, L+1), not reach(38, L+1).
reach(40,L + 1) :- reach(82, L), L <= 12, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(65, L+1), not reach(38, L+1).
reach(65,L + 1) :- reach(82, L), L <= 12, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(40, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(82, L), L <= 12, L>= 4, not reach(14, L+1), not reach(21, L+1), not reach(40, L+1), not reach(65, L+1).
reach(25,L + 1) :- reach(37, L), L <= 12, L>= 4, not reach(80, L+1), not reach(1, L+1), not reach(7, L+1), not reach(33, L+1).
reach(80,L + 1) :- reach(37, L), L <= 12, L>= 4, not reach(25, L+1), not reach(1, L+1), not reach(7, L+1), not reach(33, L+1).
reach(1,L + 1) :- reach(37, L), L <= 12, L>= 4, not reach(25, L+1), not reach(80, L+1), not reach(7, L+1), not reach(33, L+1).
reach(7,L + 1) :- reach(37, L), L <= 12, L>= 4, not reach(25, L+1), not reach(80, L+1), not reach(1, L+1), not reach(33, L+1).
reach(33,L + 1) :- reach(37, L), L <= 12, L>= 4, not reach(25, L+1), not reach(80, L+1), not reach(1, L+1), not reach(7, L+1).
reach(5,L + 1) :- reach(18, L), L <= 10, L>= 0.
reach(29,L + 1) :- reach(27, L), L <= 12, L>= 4, not reach(41, L+1), not reach(59, L+1), not reach(42, L+1), not reach(66, L+1), not reach(71, L+1).
reach(41,L + 1) :- reach(27, L), L <= 12, L>= 4, not reach(29, L+1), not reach(59, L+1), not reach(42, L+1), not reach(66, L+1), not reach(71, L+1).
reach(59,L + 1) :- reach(27, L), L <= 12, L>= 4, not reach(29, L+1), not reach(41, L+1), not reach(42, L+1), not reach(66, L+1), not reach(71, L+1).
reach(42,L + 1) :- reach(27, L), L <= 12, L>= 4, not reach(29, L+1), not reach(41, L+1), not reach(59, L+1), not reach(66, L+1), not reach(71, L+1).
reach(66,L + 1) :- reach(27, L), L <= 12, L>= 4, not reach(29, L+1), not reach(41, L+1), not reach(59, L+1), not reach(42, L+1), not reach(71, L+1).
reach(71,L + 1) :- reach(27, L), L <= 12, L>= 4, not reach(29, L+1), not reach(41, L+1), not reach(59, L+1), not reach(42, L+1), not reach(66, L+1).
reach(27,L + 1) :- reach(29, L), L <= 11, L>= 4, not reach(71, L+1), not reach(43, L+1), not reach(14, L+1).
reach(71,L + 1) :- reach(29, L), L <= 11, L>= 4, not reach(27, L+1), not reach(43, L+1), not reach(14, L+1).
reach(43,L + 1) :- reach(29, L), L <= 11, L>= 4, not reach(27, L+1), not reach(71, L+1), not reach(14, L+1).
reach(14,L + 1) :- reach(29, L), L <= 11, L>= 4, not reach(27, L+1), not reach(71, L+1), not reach(43, L+1).
reach(31,L + 1) :- reach(28, L), L <= 11, L>= 5, not reach(47, L+1), not reach(50, L+1).
reach(47,L + 1) :- reach(28, L), L <= 11, L>= 5, not reach(31, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(28, L), L <= 11, L>= 5, not reach(31, L+1), not reach(47, L+1).
reach(28,L + 1) :- reach(31, L), L <= 12, L>= 5, not reach(66, L+1), not reach(38, L+1).
reach(66,L + 1) :- reach(31, L), L <= 12, L>= 5, not reach(28, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(31, L), L <= 12, L>= 5, not reach(28, L+1), not reach(66, L+1).
reach(28,L + 1) :- reach(47, L), L <= 11, L>= 4, not reach(80, L+1), not reach(8, L+1).
reach(80,L + 1) :- reach(47, L), L <= 11, L>= 4, not reach(28, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(47, L), L <= 11, L>= 4, not reach(28, L+1), not reach(80, L+1).
reach(29,L + 1) :- reach(43, L), L <= 10, L>= 3, not reach(36, L+1), not reach(42, L+1), not reach(46, L+1).
reach(36,L + 1) :- reach(43, L), L <= 10, L>= 3, not reach(29, L+1), not reach(42, L+1), not reach(46, L+1).
reach(42,L + 1) :- reach(43, L), L <= 10, L>= 3, not reach(29, L+1), not reach(36, L+1), not reach(46, L+1).
reach(46,L + 1) :- reach(43, L), L <= 10, L>= 3, not reach(29, L+1), not reach(36, L+1), not reach(42, L+1).
reach(38,L + 1) :- reach(50, L), L <= 10, L>= 5, not reach(26, L+1), not reach(28, L+1), not reach(52, L+1).
reach(26,L + 1) :- reach(50, L), L <= 10, L>= 5, not reach(38, L+1), not reach(28, L+1), not reach(52, L+1).
reach(28,L + 1) :- reach(50, L), L <= 10, L>= 5, not reach(38, L+1), not reach(26, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(50, L), L <= 10, L>= 5, not reach(38, L+1), not reach(26, L+1), not reach(28, L+1).
reach(31,L + 1) :- reach(38, L), L <= 11, L>= 4, not reach(41, L+1), not reach(50, L+1), not reach(80, L+1), not reach(82, L+1), not reach(73, L+1).
reach(41,L + 1) :- reach(38, L), L <= 11, L>= 4, not reach(31, L+1), not reach(50, L+1), not reach(80, L+1), not reach(82, L+1), not reach(73, L+1).
reach(50,L + 1) :- reach(38, L), L <= 11, L>= 4, not reach(31, L+1), not reach(41, L+1), not reach(80, L+1), not reach(82, L+1), not reach(73, L+1).
reach(80,L + 1) :- reach(38, L), L <= 11, L>= 4, not reach(31, L+1), not reach(41, L+1), not reach(50, L+1), not reach(82, L+1), not reach(73, L+1).
reach(82,L + 1) :- reach(38, L), L <= 11, L>= 4, not reach(31, L+1), not reach(41, L+1), not reach(50, L+1), not reach(80, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(38, L), L <= 11, L>= 4, not reach(31, L+1), not reach(41, L+1), not reach(50, L+1), not reach(80, L+1), not reach(82, L+1).
reach(59,L + 1) :- reach(33, L), L <= 13, L>= 5, not reach(61, L+1), not reach(80, L+1), not reach(10, L+1), not reach(37, L+1).
reach(61,L + 1) :- reach(33, L), L <= 13, L>= 5, not reach(59, L+1), not reach(80, L+1), not reach(10, L+1), not reach(37, L+1).
reach(80,L + 1) :- reach(33, L), L <= 13, L>= 5, not reach(59, L+1), not reach(61, L+1), not reach(10, L+1), not reach(37, L+1).
reach(10,L + 1) :- reach(33, L), L <= 13, L>= 5, not reach(59, L+1), not reach(61, L+1), not reach(80, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(33, L), L <= 13, L>= 5, not reach(59, L+1), not reach(61, L+1), not reach(80, L+1), not reach(10, L+1).
reach(71,L + 1) :- reach(34, L), L <= 11, L>= 3, not reach(80, L+1), not reach(10, L+1), not reach(46, L+1), not reach(36, L+1).
reach(80,L + 1) :- reach(34, L), L <= 11, L>= 3, not reach(71, L+1), not reach(10, L+1), not reach(46, L+1), not reach(36, L+1).
reach(10,L + 1) :- reach(34, L), L <= 11, L>= 3, not reach(71, L+1), not reach(80, L+1), not reach(46, L+1), not reach(36, L+1).
reach(46,L + 1) :- reach(34, L), L <= 11, L>= 3, not reach(71, L+1), not reach(80, L+1), not reach(10, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(34, L), L <= 11, L>= 3, not reach(71, L+1), not reach(80, L+1), not reach(10, L+1), not reach(46, L+1).

