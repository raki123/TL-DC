reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13).
start(21).
goal(68).
:- reach(X, L), reach(X, L'), L != L'.
reach(4,L + 1) :- reach(1, L), L <= 10, L>= 4, not reach(2, L+1), not reach(67, L+1).
reach(2,L + 1) :- reach(1, L), L <= 10, L>= 4, not reach(4, L+1), not reach(67, L+1).
reach(67,L + 1) :- reach(1, L), L <= 10, L>= 4, not reach(4, L+1), not reach(2, L+1).
reach(1,L + 1) :- reach(4, L), L <= 11, L>= 5, not reach(38, L+1), not reach(39, L+1), not reach(49, L+1), not reach(32, L+1), not reach(36, L+1).
reach(38,L + 1) :- reach(4, L), L <= 11, L>= 5, not reach(1, L+1), not reach(39, L+1), not reach(49, L+1), not reach(32, L+1), not reach(36, L+1).
reach(39,L + 1) :- reach(4, L), L <= 11, L>= 5, not reach(1, L+1), not reach(38, L+1), not reach(49, L+1), not reach(32, L+1), not reach(36, L+1).
reach(49,L + 1) :- reach(4, L), L <= 11, L>= 5, not reach(1, L+1), not reach(38, L+1), not reach(39, L+1), not reach(32, L+1), not reach(36, L+1).
reach(32,L + 1) :- reach(4, L), L <= 11, L>= 5, not reach(1, L+1), not reach(38, L+1), not reach(39, L+1), not reach(49, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(4, L), L <= 11, L>= 5, not reach(1, L+1), not reach(38, L+1), not reach(39, L+1), not reach(49, L+1), not reach(32, L+1).
reach(1,L + 1) :- reach(2, L), L <= 9, L>= 3, not reach(52, L+1), not reach(57, L+1), not reach(80, L+1), not reach(98, L+1).
reach(52,L + 1) :- reach(2, L), L <= 9, L>= 3, not reach(1, L+1), not reach(57, L+1), not reach(80, L+1), not reach(98, L+1).
reach(57,L + 1) :- reach(2, L), L <= 9, L>= 3, not reach(1, L+1), not reach(52, L+1), not reach(80, L+1), not reach(98, L+1).
reach(80,L + 1) :- reach(2, L), L <= 9, L>= 3, not reach(1, L+1), not reach(52, L+1), not reach(57, L+1), not reach(98, L+1).
reach(98,L + 1) :- reach(2, L), L <= 9, L>= 3, not reach(1, L+1), not reach(52, L+1), not reach(57, L+1), not reach(80, L+1).
reach(2,L + 1) :- reach(52, L), L <= 8, L>= 4, not reach(8, L+1), not reach(84, L+1).
reach(8,L + 1) :- reach(52, L), L <= 8, L>= 4, not reach(2, L+1), not reach(84, L+1).
reach(84,L + 1) :- reach(52, L), L <= 8, L>= 4, not reach(2, L+1), not reach(8, L+1).
reach(2,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(26, L+1), not reach(58, L+1).
reach(26,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(2, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(2, L+1), not reach(26, L+1).
reach(2,L + 1) :- reach(80, L), L <= 8, L>= 2, not reach(11, L+1), not reach(20, L+1), not reach(74, L+1), not reach(62, L+1), not reach(30, L+1).
reach(11,L + 1) :- reach(80, L), L <= 8, L>= 2, not reach(2, L+1), not reach(20, L+1), not reach(74, L+1), not reach(62, L+1), not reach(30, L+1).
reach(20,L + 1) :- reach(80, L), L <= 8, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(74, L+1), not reach(62, L+1), not reach(30, L+1).
reach(74,L + 1) :- reach(80, L), L <= 8, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(20, L+1), not reach(62, L+1), not reach(30, L+1).
reach(62,L + 1) :- reach(80, L), L <= 8, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(20, L+1), not reach(74, L+1), not reach(30, L+1).
reach(30,L + 1) :- reach(80, L), L <= 8, L>= 2, not reach(2, L+1), not reach(11, L+1), not reach(20, L+1), not reach(74, L+1), not reach(62, L+1).
reach(2,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(16, L+1), not reach(19, L+1), not reach(43, L+1), not reach(49, L+1), not reach(74, L+1), not reach(75, L+1).
reach(16,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(2, L+1), not reach(19, L+1), not reach(43, L+1), not reach(49, L+1), not reach(74, L+1), not reach(75, L+1).
reach(19,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(2, L+1), not reach(16, L+1), not reach(43, L+1), not reach(49, L+1), not reach(74, L+1), not reach(75, L+1).
reach(43,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(2, L+1), not reach(16, L+1), not reach(19, L+1), not reach(49, L+1), not reach(74, L+1), not reach(75, L+1).
reach(49,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(2, L+1), not reach(16, L+1), not reach(19, L+1), not reach(43, L+1), not reach(74, L+1), not reach(75, L+1).
reach(74,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(2, L+1), not reach(16, L+1), not reach(19, L+1), not reach(43, L+1), not reach(49, L+1), not reach(75, L+1).
reach(75,L + 1) :- reach(98, L), L <= 10, L>= 4, not reach(2, L+1), not reach(16, L+1), not reach(19, L+1), not reach(43, L+1), not reach(49, L+1), not reach(74, L+1).
reach(38,L + 1) :- reach(3, L), L <= 9, L>= 4, not reach(67, L+1), not reach(74, L+1).
reach(67,L + 1) :- reach(3, L), L <= 9, L>= 4, not reach(38, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(3, L), L <= 9, L>= 4, not reach(38, L+1), not reach(67, L+1).
reach(3,L + 1) :- reach(38, L), L <= 10, L>= 4, not reach(4, L+1), not reach(58, L+1).
reach(4,L + 1) :- reach(38, L), L <= 10, L>= 4, not reach(3, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(38, L), L <= 10, L>= 4, not reach(3, L+1), not reach(4, L+1).
reach(3,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(6, L+1), not reach(14, L+1), not reach(78, L+1), not reach(88, L+1), not reach(1, L+1), not reach(13, L+1).
reach(6,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(3, L+1), not reach(14, L+1), not reach(78, L+1), not reach(88, L+1), not reach(1, L+1), not reach(13, L+1).
reach(14,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(3, L+1), not reach(6, L+1), not reach(78, L+1), not reach(88, L+1), not reach(1, L+1), not reach(13, L+1).
reach(78,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(3, L+1), not reach(6, L+1), not reach(14, L+1), not reach(88, L+1), not reach(1, L+1), not reach(13, L+1).
reach(88,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(3, L+1), not reach(6, L+1), not reach(14, L+1), not reach(78, L+1), not reach(1, L+1), not reach(13, L+1).
reach(1,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(3, L+1), not reach(6, L+1), not reach(14, L+1), not reach(78, L+1), not reach(88, L+1), not reach(13, L+1).
reach(13,L + 1) :- reach(67, L), L <= 10, L>= 4, not reach(3, L+1), not reach(6, L+1), not reach(14, L+1), not reach(78, L+1), not reach(88, L+1), not reach(1, L+1).
reach(3,L + 1) :- reach(74, L), L <= 9, L>= 3, not reach(80, L+1), not reach(98, L+1), not reach(76, L+1).
reach(80,L + 1) :- reach(74, L), L <= 9, L>= 3, not reach(3, L+1), not reach(98, L+1), not reach(76, L+1).
reach(98,L + 1) :- reach(74, L), L <= 9, L>= 3, not reach(3, L+1), not reach(80, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(74, L), L <= 9, L>= 3, not reach(3, L+1), not reach(80, L+1), not reach(98, L+1).
reach(4,L + 1) :- reach(39, L), L <= 10, L>= 5, not reach(33, L+1), not reach(41, L+1).
reach(33,L + 1) :- reach(39, L), L <= 10, L>= 5, not reach(4, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(39, L), L <= 10, L>= 5, not reach(4, L+1), not reach(33, L+1).
reach(4,L + 1) :- reach(49, L), L <= 10, L>= 4, not reach(14, L+1), not reach(20, L+1), not reach(98, L+1).
reach(14,L + 1) :- reach(49, L), L <= 10, L>= 4, not reach(4, L+1), not reach(20, L+1), not reach(98, L+1).
reach(20,L + 1) :- reach(49, L), L <= 10, L>= 4, not reach(4, L+1), not reach(14, L+1), not reach(98, L+1).
reach(98,L + 1) :- reach(49, L), L <= 10, L>= 4, not reach(4, L+1), not reach(14, L+1), not reach(20, L+1).
reach(43,L + 1) :- reach(32, L), L <= 12, L>= 5, not reach(68, L+1), not reach(88, L+1), not reach(4, L+1), not reach(44, L+1).
reach(68,L + 1) :- reach(32, L), L <= 12, L>= 5, not reach(43, L+1), not reach(88, L+1), not reach(4, L+1), not reach(44, L+1).
reach(88,L + 1) :- reach(32, L), L <= 12, L>= 5, not reach(43, L+1), not reach(68, L+1), not reach(4, L+1), not reach(44, L+1).
reach(4,L + 1) :- reach(32, L), L <= 12, L>= 5, not reach(43, L+1), not reach(68, L+1), not reach(88, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(32, L), L <= 12, L>= 5, not reach(43, L+1), not reach(68, L+1), not reach(88, L+1), not reach(4, L+1).
reach(67,L + 1) :- reach(6, L), L <= 9, L>= 4, not reach(73, L+1), not reach(92, L+1).
reach(73,L + 1) :- reach(6, L), L <= 9, L>= 4, not reach(67, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(6, L), L <= 9, L>= 4, not reach(67, L+1), not reach(73, L+1).
reach(6,L + 1) :- reach(73, L), L <= 8, L>= 3, not reach(65, L+1), not reach(77, L+1).
reach(65,L + 1) :- reach(73, L), L <= 8, L>= 3, not reach(6, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(73, L), L <= 8, L>= 3, not reach(6, L+1), not reach(65, L+1).
reach(6,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(40, L+1), not reach(46, L+1), not reach(47, L+1), not reach(58, L+1), not reach(75, L+1), not reach(26, L+1).
reach(40,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(6, L+1), not reach(46, L+1), not reach(47, L+1), not reach(58, L+1), not reach(75, L+1), not reach(26, L+1).
reach(46,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(6, L+1), not reach(40, L+1), not reach(47, L+1), not reach(58, L+1), not reach(75, L+1), not reach(26, L+1).
reach(47,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(6, L+1), not reach(40, L+1), not reach(46, L+1), not reach(58, L+1), not reach(75, L+1), not reach(26, L+1).
reach(58,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(6, L+1), not reach(40, L+1), not reach(46, L+1), not reach(47, L+1), not reach(75, L+1), not reach(26, L+1).
reach(75,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(6, L+1), not reach(40, L+1), not reach(46, L+1), not reach(47, L+1), not reach(58, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(6, L+1), not reach(40, L+1), not reach(46, L+1), not reach(47, L+1), not reach(58, L+1), not reach(75, L+1).
reach(40,L + 1) :- reach(7, L), L <= 9, L>= 4, not reach(41, L+1), not reach(84, L+1), not reach(18, L+1).
reach(41,L + 1) :- reach(7, L), L <= 9, L>= 4, not reach(40, L+1), not reach(84, L+1), not reach(18, L+1).
reach(84,L + 1) :- reach(7, L), L <= 9, L>= 4, not reach(40, L+1), not reach(41, L+1), not reach(18, L+1).
reach(18,L + 1) :- reach(7, L), L <= 9, L>= 4, not reach(40, L+1), not reach(41, L+1), not reach(84, L+1).
reach(7,L + 1) :- reach(40, L), L <= 8, L>= 3, not reach(65, L+1), not reach(92, L+1).
reach(65,L + 1) :- reach(40, L), L <= 8, L>= 3, not reach(7, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(40, L), L <= 8, L>= 3, not reach(7, L+1), not reach(65, L+1).
reach(7,L + 1) :- reach(41, L), L <= 9, L>= 4, not reach(14, L+1), not reach(78, L+1), not reach(39, L+1).
reach(14,L + 1) :- reach(41, L), L <= 9, L>= 4, not reach(7, L+1), not reach(78, L+1), not reach(39, L+1).
reach(78,L + 1) :- reach(41, L), L <= 9, L>= 4, not reach(7, L+1), not reach(14, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(41, L), L <= 9, L>= 4, not reach(7, L+1), not reach(14, L+1), not reach(78, L+1).
reach(7,L + 1) :- reach(84, L), L <= 8, L>= 5, not reach(26, L+1), not reach(52, L+1).
reach(26,L + 1) :- reach(84, L), L <= 8, L>= 5, not reach(7, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(84, L), L <= 8, L>= 5, not reach(7, L+1), not reach(26, L+1).
reach(52,L + 1) :- reach(8, L), L <= 8, L>= 5, not reach(23, L+1), not reach(37, L+1).
reach(23,L + 1) :- reach(8, L), L <= 8, L>= 5, not reach(52, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(8, L), L <= 8, L>= 5, not reach(52, L+1), not reach(23, L+1).
reach(20,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(33, L+1), not reach(34, L+1), not reach(41, L+1), not reach(46, L+1), not reach(67, L+1), not reach(36, L+1).
reach(33,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(20, L+1), not reach(34, L+1), not reach(41, L+1), not reach(46, L+1), not reach(67, L+1), not reach(36, L+1).
reach(34,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(20, L+1), not reach(33, L+1), not reach(41, L+1), not reach(46, L+1), not reach(67, L+1), not reach(36, L+1).
reach(41,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(20, L+1), not reach(33, L+1), not reach(34, L+1), not reach(46, L+1), not reach(67, L+1), not reach(36, L+1).
reach(46,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(20, L+1), not reach(33, L+1), not reach(34, L+1), not reach(41, L+1), not reach(67, L+1), not reach(36, L+1).
reach(67,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(20, L+1), not reach(33, L+1), not reach(34, L+1), not reach(41, L+1), not reach(46, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(78, L), L <= 9, L>= 3, not reach(20, L+1), not reach(33, L+1), not reach(34, L+1), not reach(41, L+1), not reach(46, L+1), not reach(67, L+1).
reach(33,L + 1) :- reach(11, L), L <= 8, L>= 3, not reach(76, L+1), not reach(80, L+1).
reach(76,L + 1) :- reach(11, L), L <= 8, L>= 3, not reach(33, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(11, L), L <= 8, L>= 3, not reach(33, L+1), not reach(76, L+1).
reach(11,L + 1) :- reach(33, L), L <= 9, L>= 4, not reach(78, L+1), not reach(39, L+1).
reach(78,L + 1) :- reach(33, L), L <= 9, L>= 4, not reach(11, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(33, L), L <= 9, L>= 4, not reach(11, L+1), not reach(78, L+1).
reach(11,L + 1) :- reach(76, L), L <= 9, L>= 2, not reach(14, L+1), not reach(29, L+1), not reach(74, L+1).
reach(14,L + 1) :- reach(76, L), L <= 9, L>= 2, not reach(11, L+1), not reach(29, L+1), not reach(74, L+1).
reach(29,L + 1) :- reach(76, L), L <= 9, L>= 2, not reach(11, L+1), not reach(14, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(76, L), L <= 9, L>= 2, not reach(11, L+1), not reach(14, L+1), not reach(29, L+1).
reach(13,L + 1) :- reach(34, L), L <= 8, L>= 2, not reach(78, L+1), not reach(83, L+1), not reach(30, L+1).
reach(78,L + 1) :- reach(34, L), L <= 8, L>= 2, not reach(13, L+1), not reach(83, L+1), not reach(30, L+1).
reach(83,L + 1) :- reach(34, L), L <= 8, L>= 2, not reach(13, L+1), not reach(78, L+1), not reach(30, L+1).
reach(30,L + 1) :- reach(34, L), L <= 8, L>= 2, not reach(13, L+1), not reach(78, L+1), not reach(83, L+1).
reach(46,L + 1) :- reach(83, L), L <= 8, L>= 3, not reach(34, L+1), not reach(48, L+1).
reach(34,L + 1) :- reach(83, L), L <= 8, L>= 3, not reach(46, L+1), not reach(48, L+1).
reach(48,L + 1) :- reach(83, L), L <= 8, L>= 3, not reach(46, L+1), not reach(34, L+1).
reach(34,L + 1) :- reach(13, L), L <= 9, L>= 3, not reach(47, L+1), not reach(67, L+1).
reach(47,L + 1) :- reach(13, L), L <= 9, L>= 3, not reach(34, L+1), not reach(67, L+1).
reach(67,L + 1) :- reach(13, L), L <= 9, L>= 3, not reach(34, L+1), not reach(47, L+1).
reach(13,L + 1) :- reach(47, L), L <= 8, L>= 4, not reach(17, L+1), not reach(92, L+1).
reach(17,L + 1) :- reach(47, L), L <= 8, L>= 4, not reach(13, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(47, L), L <= 8, L>= 4, not reach(13, L+1), not reach(17, L+1).
reach(41,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(36, L+1), not reach(49, L+1), not reach(17, L+1), not reach(67, L+1), not reach(76, L+1), not reach(44, L+1), not reach(75, L+1).
reach(36,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(49, L+1), not reach(17, L+1), not reach(67, L+1), not reach(76, L+1), not reach(44, L+1), not reach(75, L+1).
reach(49,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(36, L+1), not reach(17, L+1), not reach(67, L+1), not reach(76, L+1), not reach(44, L+1), not reach(75, L+1).
reach(17,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(36, L+1), not reach(49, L+1), not reach(67, L+1), not reach(76, L+1), not reach(44, L+1), not reach(75, L+1).
reach(67,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(36, L+1), not reach(49, L+1), not reach(17, L+1), not reach(76, L+1), not reach(44, L+1), not reach(75, L+1).
reach(76,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(36, L+1), not reach(49, L+1), not reach(17, L+1), not reach(67, L+1), not reach(44, L+1), not reach(75, L+1).
reach(44,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(36, L+1), not reach(49, L+1), not reach(17, L+1), not reach(67, L+1), not reach(76, L+1), not reach(75, L+1).
reach(75,L + 1) :- reach(14, L), L <= 10, L>= 3, not reach(41, L+1), not reach(36, L+1), not reach(49, L+1), not reach(17, L+1), not reach(67, L+1), not reach(76, L+1), not reach(44, L+1).
reach(14,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(27, L+1), not reach(77, L+1), not reach(91, L+1), not reach(48, L+1), not reach(78, L+1), not reach(4, L+1).
reach(27,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(14, L+1), not reach(77, L+1), not reach(91, L+1), not reach(48, L+1), not reach(78, L+1), not reach(4, L+1).
reach(77,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(14, L+1), not reach(27, L+1), not reach(91, L+1), not reach(48, L+1), not reach(78, L+1), not reach(4, L+1).
reach(91,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(14, L+1), not reach(27, L+1), not reach(77, L+1), not reach(48, L+1), not reach(78, L+1), not reach(4, L+1).
reach(48,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(14, L+1), not reach(27, L+1), not reach(77, L+1), not reach(91, L+1), not reach(78, L+1), not reach(4, L+1).
reach(78,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(14, L+1), not reach(27, L+1), not reach(77, L+1), not reach(91, L+1), not reach(48, L+1), not reach(4, L+1).
reach(4,L + 1) :- reach(36, L), L <= 10, L>= 4, not reach(14, L+1), not reach(27, L+1), not reach(77, L+1), not reach(91, L+1), not reach(48, L+1), not reach(78, L+1).
reach(14,L + 1) :- reach(17, L), L <= 9, L>= 3, not reach(47, L+1), not reach(81, L+1).
reach(47,L + 1) :- reach(17, L), L <= 9, L>= 3, not reach(14, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(17, L), L <= 9, L>= 3, not reach(14, L+1), not reach(47, L+1).
reach(14,L + 1) :- reach(44, L), L <= 11, L>= 4, not reach(91, L+1), not reach(32, L+1).
reach(91,L + 1) :- reach(44, L), L <= 11, L>= 4, not reach(14, L+1), not reach(32, L+1).
reach(32,L + 1) :- reach(44, L), L <= 11, L>= 4, not reach(14, L+1), not reach(91, L+1).
reach(27,L + 1) :- reach(48, L), L <= 9, L>= 4, not reach(36, L+1), not reach(83, L+1).
reach(36,L + 1) :- reach(48, L), L <= 9, L>= 4, not reach(27, L+1), not reach(83, L+1).
reach(83,L + 1) :- reach(48, L), L <= 9, L>= 4, not reach(27, L+1), not reach(36, L+1).
reach(77,L + 1) :- reach(16, L), L <= 9, L>= 5, not reach(98, L+1), not reach(46, L+1).
reach(98,L + 1) :- reach(16, L), L <= 9, L>= 5, not reach(77, L+1), not reach(46, L+1).
reach(46,L + 1) :- reach(16, L), L <= 9, L>= 5, not reach(77, L+1), not reach(98, L+1).
reach(16,L + 1) :- reach(77, L), L <= 9, L>= 4, not reach(36, L+1), not reach(73, L+1).
reach(36,L + 1) :- reach(77, L), L <= 9, L>= 4, not reach(16, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(77, L), L <= 9, L>= 4, not reach(16, L+1), not reach(36, L+1).
reach(17,L + 1) :- reach(81, L), L <= 8, L>= 2, not reach(29, L+1), not reach(58, L+1).
reach(29,L + 1) :- reach(81, L), L <= 8, L>= 2, not reach(17, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(81, L), L <= 8, L>= 2, not reach(17, L+1), not reach(29, L+1).
reach(62,L + 1) :- reach(18, L), L <= 10, L>= 4, not reach(7, L+1), not reach(43, L+1).
reach(7,L + 1) :- reach(18, L), L <= 10, L>= 4, not reach(62, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(18, L), L <= 10, L>= 4, not reach(62, L+1), not reach(7, L+1).
reach(18,L + 1) :- reach(62, L), L <= 9, L>= 3, not reach(19, L+1), not reach(20, L+1), not reach(23, L+1), not reach(27, L+1), not reach(80, L+1).
reach(19,L + 1) :- reach(62, L), L <= 9, L>= 3, not reach(18, L+1), not reach(20, L+1), not reach(23, L+1), not reach(27, L+1), not reach(80, L+1).
reach(20,L + 1) :- reach(62, L), L <= 9, L>= 3, not reach(18, L+1), not reach(19, L+1), not reach(23, L+1), not reach(27, L+1), not reach(80, L+1).
reach(23,L + 1) :- reach(62, L), L <= 9, L>= 3, not reach(18, L+1), not reach(19, L+1), not reach(20, L+1), not reach(27, L+1), not reach(80, L+1).
reach(27,L + 1) :- reach(62, L), L <= 9, L>= 3, not reach(18, L+1), not reach(19, L+1), not reach(20, L+1), not reach(23, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(62, L), L <= 9, L>= 3, not reach(18, L+1), not reach(19, L+1), not reach(20, L+1), not reach(23, L+1), not reach(27, L+1).
reach(46,L + 1) :- reach(19, L), L <= 9, L>= 4, not reach(62, L+1), not reach(98, L+1), not reach(94, L+1).
reach(62,L + 1) :- reach(19, L), L <= 9, L>= 4, not reach(46, L+1), not reach(98, L+1), not reach(94, L+1).
reach(98,L + 1) :- reach(19, L), L <= 9, L>= 4, not reach(46, L+1), not reach(62, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(19, L), L <= 9, L>= 4, not reach(46, L+1), not reach(62, L+1), not reach(98, L+1).
reach(19,L + 1) :- reach(46, L), L <= 8, L>= 4, not reach(78, L+1), not reach(83, L+1), not reach(92, L+1), not reach(16, L+1).
reach(78,L + 1) :- reach(46, L), L <= 8, L>= 4, not reach(19, L+1), not reach(83, L+1), not reach(92, L+1), not reach(16, L+1).
reach(83,L + 1) :- reach(46, L), L <= 8, L>= 4, not reach(19, L+1), not reach(78, L+1), not reach(92, L+1), not reach(16, L+1).
reach(92,L + 1) :- reach(46, L), L <= 8, L>= 4, not reach(19, L+1), not reach(78, L+1), not reach(83, L+1), not reach(16, L+1).
reach(16,L + 1) :- reach(46, L), L <= 8, L>= 4, not reach(19, L+1), not reach(78, L+1), not reach(83, L+1), not reach(92, L+1).
reach(49,L + 1) :- reach(20, L), L <= 9, L>= 3, not reach(62, L+1), not reach(78, L+1), not reach(80, L+1), not reach(91, L+1).
reach(62,L + 1) :- reach(20, L), L <= 9, L>= 3, not reach(49, L+1), not reach(78, L+1), not reach(80, L+1), not reach(91, L+1).
reach(78,L + 1) :- reach(20, L), L <= 9, L>= 3, not reach(49, L+1), not reach(62, L+1), not reach(80, L+1), not reach(91, L+1).
reach(80,L + 1) :- reach(20, L), L <= 9, L>= 3, not reach(49, L+1), not reach(62, L+1), not reach(78, L+1), not reach(91, L+1).
reach(91,L + 1) :- reach(20, L), L <= 9, L>= 3, not reach(49, L+1), not reach(62, L+1), not reach(78, L+1), not reach(80, L+1).
reach(20,L + 1) :- reach(91, L), L <= 10, L>= 4, not reach(36, L+1), not reach(44, L+1).
reach(36,L + 1) :- reach(91, L), L <= 10, L>= 4, not reach(20, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(91, L), L <= 10, L>= 4, not reach(20, L+1), not reach(36, L+1).
reach(29,L + 1) :- reach(21, L), L <= 7, L>= 0, not reach(30, L+1).
reach(30,L + 1) :- reach(21, L), L <= 7, L>= 0, not reach(29, L+1).
reach(58,L + 1) :- reach(23, L), L <= 9, L>= 4, not reach(54, L+1), not reach(62, L+1), not reach(8, L+1).
reach(54,L + 1) :- reach(23, L), L <= 9, L>= 4, not reach(58, L+1), not reach(62, L+1), not reach(8, L+1).
reach(62,L + 1) :- reach(23, L), L <= 9, L>= 4, not reach(58, L+1), not reach(54, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(23, L), L <= 9, L>= 4, not reach(58, L+1), not reach(54, L+1), not reach(62, L+1).
reach(23,L + 1) :- reach(58, L), L <= 9, L>= 3, not reach(38, L+1), not reach(57, L+1), not reach(81, L+1), not reach(92, L+1), not reach(94, L+1).
reach(38,L + 1) :- reach(58, L), L <= 9, L>= 3, not reach(23, L+1), not reach(57, L+1), not reach(81, L+1), not reach(92, L+1), not reach(94, L+1).
reach(57,L + 1) :- reach(58, L), L <= 9, L>= 3, not reach(23, L+1), not reach(38, L+1), not reach(81, L+1), not reach(92, L+1), not reach(94, L+1).
reach(81,L + 1) :- reach(58, L), L <= 9, L>= 3, not reach(23, L+1), not reach(38, L+1), not reach(57, L+1), not reach(92, L+1), not reach(94, L+1).
reach(92,L + 1) :- reach(58, L), L <= 9, L>= 3, not reach(23, L+1), not reach(38, L+1), not reach(57, L+1), not reach(81, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(58, L), L <= 9, L>= 3, not reach(23, L+1), not reach(38, L+1), not reach(57, L+1), not reach(81, L+1), not reach(92, L+1).
reach(23,L + 1) :- reach(54, L), L <= 10, L>= 5, not reach(26, L+1), not reach(88, L+1).
reach(26,L + 1) :- reach(54, L), L <= 10, L>= 5, not reach(23, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(54, L), L <= 10, L>= 5, not reach(23, L+1), not reach(26, L+1).
reach(57,L + 1) :- reach(26, L), L <= 9, L>= 5, not reach(54, L+1), not reach(37, L+1), not reach(84, L+1), not reach(92, L+1).
reach(54,L + 1) :- reach(26, L), L <= 9, L>= 5, not reach(57, L+1), not reach(37, L+1), not reach(84, L+1), not reach(92, L+1).
reach(37,L + 1) :- reach(26, L), L <= 9, L>= 5, not reach(57, L+1), not reach(54, L+1), not reach(84, L+1), not reach(92, L+1).
reach(84,L + 1) :- reach(26, L), L <= 9, L>= 5, not reach(57, L+1), not reach(54, L+1), not reach(37, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(26, L), L <= 9, L>= 5, not reach(57, L+1), not reach(54, L+1), not reach(37, L+1), not reach(84, L+1).
reach(26,L + 1) :- reach(37, L), L <= 8, L>= 4, not reach(94, L+1), not reach(8, L+1).
reach(94,L + 1) :- reach(37, L), L <= 8, L>= 4, not reach(26, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(37, L), L <= 8, L>= 4, not reach(26, L+1), not reach(94, L+1).
reach(48,L + 1) :- reach(27, L), L <= 9, L>= 3, not reach(65, L+1), not reach(36, L+1), not reach(62, L+1).
reach(65,L + 1) :- reach(27, L), L <= 9, L>= 3, not reach(48, L+1), not reach(36, L+1), not reach(62, L+1).
reach(36,L + 1) :- reach(27, L), L <= 9, L>= 3, not reach(48, L+1), not reach(65, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(27, L), L <= 9, L>= 3, not reach(48, L+1), not reach(65, L+1), not reach(36, L+1).
reach(27,L + 1) :- reach(65, L), L <= 8, L>= 2, not reach(30, L+1), not reach(40, L+1), not reach(73, L+1), not reach(94, L+1).
reach(30,L + 1) :- reach(65, L), L <= 8, L>= 2, not reach(27, L+1), not reach(40, L+1), not reach(73, L+1), not reach(94, L+1).
reach(40,L + 1) :- reach(65, L), L <= 8, L>= 2, not reach(27, L+1), not reach(30, L+1), not reach(73, L+1), not reach(94, L+1).
reach(73,L + 1) :- reach(65, L), L <= 8, L>= 2, not reach(27, L+1), not reach(30, L+1), not reach(40, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(65, L), L <= 8, L>= 2, not reach(27, L+1), not reach(30, L+1), not reach(40, L+1), not reach(73, L+1).
reach(81,L + 1) :- reach(29, L), L <= 8, L>= 1, not reach(21, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(29, L), L <= 8, L>= 1, not reach(81, L+1), not reach(21, L+1).
reach(65,L + 1) :- reach(30, L), L <= 7, L>= 1, not reach(21, L+1), not reach(34, L+1), not reach(80, L+1).
reach(34,L + 1) :- reach(30, L), L <= 7, L>= 1, not reach(65, L+1), not reach(21, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(30, L), L <= 7, L>= 1, not reach(65, L+1), not reach(21, L+1), not reach(34, L+1).
reach(32,L + 1) :- reach(43, L), L <= 11, L>= 5, not reach(98, L+1), not reach(18, L+1).
reach(98,L + 1) :- reach(43, L), L <= 11, L>= 5, not reach(32, L+1), not reach(18, L+1).
reach(18,L + 1) :- reach(43, L), L <= 11, L>= 5, not reach(32, L+1), not reach(98, L+1).
reach(32,L + 1) :- reach(88, L), L <= 11, L>= 5, not reach(54, L+1), not reach(67, L+1).
reach(54,L + 1) :- reach(88, L), L <= 11, L>= 5, not reach(32, L+1), not reach(67, L+1).
reach(67,L + 1) :- reach(88, L), L <= 11, L>= 5, not reach(32, L+1), not reach(54, L+1).
reach(92,L + 1) :- reach(75, L), L <= 9, L>= 4, not reach(98, L+1), not reach(14, L+1).
reach(98,L + 1) :- reach(75, L), L <= 9, L>= 4, not reach(92, L+1), not reach(14, L+1).
reach(14,L + 1) :- reach(75, L), L <= 9, L>= 4, not reach(92, L+1), not reach(98, L+1).
reach(37,L + 1) :- reach(94, L), L <= 8, L>= 3, not reach(58, L+1), not reach(65, L+1), not reach(19, L+1).
reach(58,L + 1) :- reach(94, L), L <= 8, L>= 3, not reach(37, L+1), not reach(65, L+1), not reach(19, L+1).
reach(65,L + 1) :- reach(94, L), L <= 8, L>= 3, not reach(37, L+1), not reach(58, L+1), not reach(19, L+1).
reach(19,L + 1) :- reach(94, L), L <= 8, L>= 3, not reach(37, L+1), not reach(58, L+1), not reach(65, L+1).

