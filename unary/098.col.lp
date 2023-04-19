reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17), not reach(X, 18).
start(4).
goal(48).
:- reach(X, L), reach(X, L'), L != L'.
reach(2,L + 1) :- reach(1, L), L <= 14, L>= 3, not reach(9, L+1), not reach(96, L+1).
reach(9,L + 1) :- reach(1, L), L <= 14, L>= 3, not reach(2, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(1, L), L <= 14, L>= 3, not reach(2, L+1), not reach(9, L+1).
reach(1,L + 1) :- reach(2, L), L <= 13, L>= 2, not reach(3, L+1), not reach(79, L+1).
reach(3,L + 1) :- reach(2, L), L <= 13, L>= 2, not reach(1, L+1), not reach(79, L+1).
reach(79,L + 1) :- reach(2, L), L <= 13, L>= 2, not reach(1, L+1), not reach(3, L+1).
reach(1,L + 1) :- reach(9, L), L <= 15, L>= 3, not reach(10, L+1), not reach(62, L+1), not reach(90, L+1).
reach(10,L + 1) :- reach(9, L), L <= 15, L>= 3, not reach(1, L+1), not reach(62, L+1), not reach(90, L+1).
reach(62,L + 1) :- reach(9, L), L <= 15, L>= 3, not reach(1, L+1), not reach(10, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(9, L), L <= 15, L>= 3, not reach(1, L+1), not reach(10, L+1), not reach(62, L+1).
reach(1,L + 1) :- reach(96, L), L <= 14, L>= 4, not reach(31, L+1), not reach(90, L+1).
reach(31,L + 1) :- reach(96, L), L <= 14, L>= 4, not reach(1, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(96, L), L <= 14, L>= 4, not reach(1, L+1), not reach(31, L+1).
reach(2,L + 1) :- reach(3, L), L <= 14, L>= 3, not reach(11, L+1), not reach(54, L+1), not reach(47, L+1).
reach(11,L + 1) :- reach(3, L), L <= 14, L>= 3, not reach(2, L+1), not reach(54, L+1), not reach(47, L+1).
reach(54,L + 1) :- reach(3, L), L <= 14, L>= 3, not reach(2, L+1), not reach(11, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(3, L), L <= 14, L>= 3, not reach(2, L+1), not reach(11, L+1), not reach(54, L+1).
reach(2,L + 1) :- reach(79, L), L <= 13, L>= 1, not reach(47, L+1), not reach(61, L+1), not reach(62, L+1), not reach(66, L+1), not reach(13, L+1), not reach(4, L+1).
reach(47,L + 1) :- reach(79, L), L <= 13, L>= 1, not reach(2, L+1), not reach(61, L+1), not reach(62, L+1), not reach(66, L+1), not reach(13, L+1), not reach(4, L+1).
reach(61,L + 1) :- reach(79, L), L <= 13, L>= 1, not reach(2, L+1), not reach(47, L+1), not reach(62, L+1), not reach(66, L+1), not reach(13, L+1), not reach(4, L+1).
reach(62,L + 1) :- reach(79, L), L <= 13, L>= 1, not reach(2, L+1), not reach(47, L+1), not reach(61, L+1), not reach(66, L+1), not reach(13, L+1), not reach(4, L+1).
reach(66,L + 1) :- reach(79, L), L <= 13, L>= 1, not reach(2, L+1), not reach(47, L+1), not reach(61, L+1), not reach(62, L+1), not reach(13, L+1), not reach(4, L+1).
reach(13,L + 1) :- reach(79, L), L <= 13, L>= 1, not reach(2, L+1), not reach(47, L+1), not reach(61, L+1), not reach(62, L+1), not reach(66, L+1), not reach(4, L+1).
reach(3,L + 1) :- reach(11, L), L <= 13, L>= 2, not reach(17, L+1), not reach(94, L+1), not reach(6, L+1).
reach(17,L + 1) :- reach(11, L), L <= 13, L>= 2, not reach(3, L+1), not reach(94, L+1), not reach(6, L+1).
reach(94,L + 1) :- reach(11, L), L <= 13, L>= 2, not reach(3, L+1), not reach(17, L+1), not reach(6, L+1).
reach(6,L + 1) :- reach(11, L), L <= 13, L>= 2, not reach(3, L+1), not reach(17, L+1), not reach(94, L+1).
reach(3,L + 1) :- reach(54, L), L <= 15, L>= 3, not reach(32, L+1), not reach(43, L+1), not reach(47, L+1), not reach(62, L+1), not reach(10, L+1).
reach(32,L + 1) :- reach(54, L), L <= 15, L>= 3, not reach(3, L+1), not reach(43, L+1), not reach(47, L+1), not reach(62, L+1), not reach(10, L+1).
reach(43,L + 1) :- reach(54, L), L <= 15, L>= 3, not reach(3, L+1), not reach(32, L+1), not reach(47, L+1), not reach(62, L+1), not reach(10, L+1).
reach(47,L + 1) :- reach(54, L), L <= 15, L>= 3, not reach(3, L+1), not reach(32, L+1), not reach(43, L+1), not reach(62, L+1), not reach(10, L+1).
reach(62,L + 1) :- reach(54, L), L <= 15, L>= 3, not reach(3, L+1), not reach(32, L+1), not reach(43, L+1), not reach(47, L+1), not reach(10, L+1).
reach(10,L + 1) :- reach(54, L), L <= 15, L>= 3, not reach(3, L+1), not reach(32, L+1), not reach(43, L+1), not reach(47, L+1), not reach(62, L+1).
reach(6,L + 1) :- reach(4, L), L <= 12, L>= 0, not reach(79, L+1).
reach(79,L + 1) :- reach(4, L), L <= 12, L>= 0, not reach(6, L+1).
reach(11,L + 1) :- reach(6, L), L <= 12, L>= 1, not reach(4, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(6, L), L <= 12, L>= 1, not reach(4, L+1), not reach(11, L+1).
reach(38,L + 1) :- reach(43, L), L <= 14, L>= 3, not reach(54, L+1), not reach(62, L+1), not reach(94, L+1).
reach(54,L + 1) :- reach(43, L), L <= 14, L>= 3, not reach(38, L+1), not reach(62, L+1), not reach(94, L+1).
reach(62,L + 1) :- reach(43, L), L <= 14, L>= 3, not reach(38, L+1), not reach(54, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(43, L), L <= 14, L>= 3, not reach(38, L+1), not reach(54, L+1), not reach(62, L+1).
reach(11,L + 1) :- reach(94, L), L <= 13, L>= 3, not reach(14, L+1), not reach(43, L+1), not reach(76, L+1).
reach(14,L + 1) :- reach(94, L), L <= 13, L>= 3, not reach(11, L+1), not reach(43, L+1), not reach(76, L+1).
reach(43,L + 1) :- reach(94, L), L <= 13, L>= 3, not reach(11, L+1), not reach(14, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(94, L), L <= 13, L>= 3, not reach(11, L+1), not reach(14, L+1), not reach(43, L+1).
reach(14,L + 1) :- reach(8, L), L <= 14, L>= 4, not reach(18, L+1), not reach(57, L+1).
reach(18,L + 1) :- reach(8, L), L <= 14, L>= 4, not reach(14, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(8, L), L <= 14, L>= 4, not reach(14, L+1), not reach(18, L+1).
reach(8,L + 1) :- reach(14, L), L <= 13, L>= 3, not reach(90, L+1), not reach(94, L+1), not reach(62, L+1).
reach(90,L + 1) :- reach(14, L), L <= 13, L>= 3, not reach(8, L+1), not reach(94, L+1), not reach(62, L+1).
reach(94,L + 1) :- reach(14, L), L <= 13, L>= 3, not reach(8, L+1), not reach(90, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(14, L), L <= 13, L>= 3, not reach(8, L+1), not reach(90, L+1), not reach(94, L+1).
reach(8,L + 1) :- reach(18, L), L <= 15, L>= 5, not reach(30, L+1), not reach(88, L+1).
reach(30,L + 1) :- reach(18, L), L <= 15, L>= 5, not reach(8, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(18, L), L <= 15, L>= 5, not reach(8, L+1), not reach(30, L+1).
reach(9,L + 1) :- reach(10, L), L <= 16, L>= 4, not reach(15, L+1), not reach(54, L+1), not reach(74, L+1).
reach(15,L + 1) :- reach(10, L), L <= 16, L>= 4, not reach(9, L+1), not reach(54, L+1), not reach(74, L+1).
reach(54,L + 1) :- reach(10, L), L <= 16, L>= 4, not reach(9, L+1), not reach(15, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(10, L), L <= 16, L>= 4, not reach(9, L+1), not reach(15, L+1), not reach(54, L+1).
reach(9,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(15, L+1), not reach(21, L+1), not reach(43, L+1), not reach(54, L+1), not reach(79, L+1), not reach(14, L+1), not reach(31, L+1).
reach(15,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(21, L+1), not reach(43, L+1), not reach(54, L+1), not reach(79, L+1), not reach(14, L+1), not reach(31, L+1).
reach(21,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(15, L+1), not reach(43, L+1), not reach(54, L+1), not reach(79, L+1), not reach(14, L+1), not reach(31, L+1).
reach(43,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(15, L+1), not reach(21, L+1), not reach(54, L+1), not reach(79, L+1), not reach(14, L+1), not reach(31, L+1).
reach(54,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(15, L+1), not reach(21, L+1), not reach(43, L+1), not reach(79, L+1), not reach(14, L+1), not reach(31, L+1).
reach(79,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(15, L+1), not reach(21, L+1), not reach(43, L+1), not reach(54, L+1), not reach(14, L+1), not reach(31, L+1).
reach(14,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(15, L+1), not reach(21, L+1), not reach(43, L+1), not reach(54, L+1), not reach(79, L+1), not reach(31, L+1).
reach(31,L + 1) :- reach(62, L), L <= 14, L>= 2, not reach(9, L+1), not reach(15, L+1), not reach(21, L+1), not reach(43, L+1), not reach(54, L+1), not reach(79, L+1), not reach(14, L+1).
reach(9,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(14, L+1), not reach(41, L+1), not reach(58, L+1), not reach(64, L+1), not reach(65, L+1), not reach(66, L+1), not reach(87, L+1), not reach(96, L+1).
reach(14,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(41, L+1), not reach(58, L+1), not reach(64, L+1), not reach(65, L+1), not reach(66, L+1), not reach(87, L+1), not reach(96, L+1).
reach(41,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(58, L+1), not reach(64, L+1), not reach(65, L+1), not reach(66, L+1), not reach(87, L+1), not reach(96, L+1).
reach(58,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(41, L+1), not reach(64, L+1), not reach(65, L+1), not reach(66, L+1), not reach(87, L+1), not reach(96, L+1).
reach(64,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(41, L+1), not reach(58, L+1), not reach(65, L+1), not reach(66, L+1), not reach(87, L+1), not reach(96, L+1).
reach(65,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(41, L+1), not reach(58, L+1), not reach(64, L+1), not reach(66, L+1), not reach(87, L+1), not reach(96, L+1).
reach(66,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(41, L+1), not reach(58, L+1), not reach(64, L+1), not reach(65, L+1), not reach(87, L+1), not reach(96, L+1).
reach(87,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(41, L+1), not reach(58, L+1), not reach(64, L+1), not reach(65, L+1), not reach(66, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(90, L), L <= 14, L>= 3, not reach(9, L+1), not reach(14, L+1), not reach(41, L+1), not reach(58, L+1), not reach(64, L+1), not reach(65, L+1), not reach(66, L+1), not reach(87, L+1).
reach(10,L + 1) :- reach(15, L), L <= 15, L>= 3, not reach(62, L+1), not reach(35, L+1).
reach(62,L + 1) :- reach(15, L), L <= 15, L>= 3, not reach(10, L+1), not reach(35, L+1).
reach(35,L + 1) :- reach(15, L), L <= 15, L>= 3, not reach(10, L+1), not reach(62, L+1).
reach(11,L + 1) :- reach(17, L), L <= 12, L>= 3, not reach(21, L+1), not reach(57, L+1).
reach(21,L + 1) :- reach(17, L), L <= 12, L>= 3, not reach(11, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(17, L), L <= 12, L>= 3, not reach(11, L+1), not reach(21, L+1).
reach(21,L + 1) :- reach(13, L), L <= 14, L>= 2, not reach(60, L+1), not reach(85, L+1), not reach(79, L+1).
reach(60,L + 1) :- reach(13, L), L <= 14, L>= 2, not reach(21, L+1), not reach(85, L+1), not reach(79, L+1).
reach(85,L + 1) :- reach(13, L), L <= 14, L>= 2, not reach(21, L+1), not reach(60, L+1), not reach(79, L+1).
reach(79,L + 1) :- reach(13, L), L <= 14, L>= 2, not reach(21, L+1), not reach(60, L+1), not reach(85, L+1).
reach(13,L + 1) :- reach(21, L), L <= 13, L>= 3, not reach(17, L+1), not reach(62, L+1), not reach(42, L+1), not reach(26, L+1).
reach(17,L + 1) :- reach(21, L), L <= 13, L>= 3, not reach(13, L+1), not reach(62, L+1), not reach(42, L+1), not reach(26, L+1).
reach(62,L + 1) :- reach(21, L), L <= 13, L>= 3, not reach(13, L+1), not reach(17, L+1), not reach(42, L+1), not reach(26, L+1).
reach(42,L + 1) :- reach(21, L), L <= 13, L>= 3, not reach(13, L+1), not reach(17, L+1), not reach(62, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(21, L), L <= 13, L>= 3, not reach(13, L+1), not reach(17, L+1), not reach(62, L+1), not reach(42, L+1).
reach(13,L + 1) :- reach(60, L), L <= 15, L>= 3, not reach(55, L+1), not reach(58, L+1).
reach(55,L + 1) :- reach(60, L), L <= 15, L>= 3, not reach(13, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(60, L), L <= 15, L>= 3, not reach(13, L+1), not reach(55, L+1).
reach(13,L + 1) :- reach(85, L), L <= 15, L>= 3, not reach(51, L+1), not reach(88, L+1).
reach(51,L + 1) :- reach(85, L), L <= 15, L>= 3, not reach(13, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(85, L), L <= 15, L>= 3, not reach(13, L+1), not reach(51, L+1).
reach(15,L + 1) :- reach(35, L), L <= 15, L>= 4, not reach(30, L+1), not reach(88, L+1).
reach(30,L + 1) :- reach(35, L), L <= 15, L>= 4, not reach(15, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(35, L), L <= 15, L>= 4, not reach(15, L+1), not reach(30, L+1).
reach(17,L + 1) :- reach(57, L), L <= 13, L>= 4, not reach(50, L+1), not reach(8, L+1).
reach(50,L + 1) :- reach(57, L), L <= 13, L>= 4, not reach(17, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(57, L), L <= 13, L>= 4, not reach(17, L+1), not reach(50, L+1).
reach(18,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(55, L+1), not reach(75, L+1), not reach(50, L+1), not reach(33, L+1), not reach(41, L+1).
reach(35,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(37, L+1), not reach(38, L+1), not reach(55, L+1), not reach(75, L+1), not reach(50, L+1), not reach(33, L+1), not reach(41, L+1).
reach(37,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(38, L+1), not reach(55, L+1), not reach(75, L+1), not reach(50, L+1), not reach(33, L+1), not reach(41, L+1).
reach(38,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(37, L+1), not reach(55, L+1), not reach(75, L+1), not reach(50, L+1), not reach(33, L+1), not reach(41, L+1).
reach(55,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(75, L+1), not reach(50, L+1), not reach(33, L+1), not reach(41, L+1).
reach(75,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(55, L+1), not reach(50, L+1), not reach(33, L+1), not reach(41, L+1).
reach(50,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(55, L+1), not reach(75, L+1), not reach(33, L+1), not reach(41, L+1).
reach(33,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(55, L+1), not reach(75, L+1), not reach(50, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(30, L), L <= 15, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(37, L+1), not reach(38, L+1), not reach(55, L+1), not reach(75, L+1), not reach(50, L+1), not reach(33, L+1).
reach(18,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(25, L+1), not reach(35, L+1), not reach(83, L+1), not reach(85, L+1), not reach(41, L+1), not reach(51, L+1), not reach(74, L+1).
reach(25,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(35, L+1), not reach(83, L+1), not reach(85, L+1), not reach(41, L+1), not reach(51, L+1), not reach(74, L+1).
reach(35,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(25, L+1), not reach(83, L+1), not reach(85, L+1), not reach(41, L+1), not reach(51, L+1), not reach(74, L+1).
reach(83,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(25, L+1), not reach(35, L+1), not reach(85, L+1), not reach(41, L+1), not reach(51, L+1), not reach(74, L+1).
reach(85,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(25, L+1), not reach(35, L+1), not reach(83, L+1), not reach(41, L+1), not reach(51, L+1), not reach(74, L+1).
reach(41,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(25, L+1), not reach(35, L+1), not reach(83, L+1), not reach(85, L+1), not reach(51, L+1), not reach(74, L+1).
reach(51,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(25, L+1), not reach(35, L+1), not reach(83, L+1), not reach(85, L+1), not reach(41, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(18, L+1), not reach(25, L+1), not reach(35, L+1), not reach(83, L+1), not reach(85, L+1), not reach(41, L+1), not reach(51, L+1).
reach(20,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(21, L+1), not reach(32, L+1), not reach(6, L+1).
reach(21,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(20, L+1), not reach(32, L+1), not reach(6, L+1).
reach(32,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(20, L+1), not reach(21, L+1), not reach(6, L+1).
reach(6,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(20, L+1), not reach(21, L+1), not reach(32, L+1).
reach(25,L + 1) :- reach(20, L), L <= 14, L>= 3, not reach(26, L+1), not reach(32, L+1).
reach(26,L + 1) :- reach(20, L), L <= 14, L>= 3, not reach(25, L+1), not reach(32, L+1).
reach(32,L + 1) :- reach(20, L), L <= 14, L>= 3, not reach(25, L+1), not reach(26, L+1).
reach(20,L + 1) :- reach(25, L), L <= 15, L>= 4, not reach(31, L+1), not reach(50, L+1), not reach(88, L+1).
reach(31,L + 1) :- reach(25, L), L <= 15, L>= 4, not reach(20, L+1), not reach(50, L+1), not reach(88, L+1).
reach(50,L + 1) :- reach(25, L), L <= 15, L>= 4, not reach(20, L+1), not reach(31, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(25, L), L <= 15, L>= 4, not reach(20, L+1), not reach(31, L+1), not reach(50, L+1).
reach(21,L + 1) :- reach(42, L), L <= 13, L>= 4, not reach(33, L+1), not reach(37, L+1).
reach(33,L + 1) :- reach(42, L), L <= 13, L>= 4, not reach(21, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(42, L), L <= 13, L>= 4, not reach(21, L+1), not reach(33, L+1).
reach(33,L + 1) :- reach(29, L), L <= 13, L>= 3, not reach(66, L+1), not reach(76, L+1).
reach(66,L + 1) :- reach(29, L), L <= 13, L>= 3, not reach(33, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(29, L), L <= 13, L>= 3, not reach(33, L+1), not reach(66, L+1).
reach(25,L + 1) :- reach(31, L), L <= 15, L>= 3, not reach(55, L+1), not reach(75, L+1), not reach(96, L+1), not reach(62, L+1).
reach(55,L + 1) :- reach(31, L), L <= 15, L>= 3, not reach(25, L+1), not reach(75, L+1), not reach(96, L+1), not reach(62, L+1).
reach(75,L + 1) :- reach(31, L), L <= 15, L>= 3, not reach(25, L+1), not reach(55, L+1), not reach(96, L+1), not reach(62, L+1).
reach(96,L + 1) :- reach(31, L), L <= 15, L>= 3, not reach(25, L+1), not reach(55, L+1), not reach(75, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(31, L), L <= 15, L>= 3, not reach(25, L+1), not reach(55, L+1), not reach(75, L+1), not reach(96, L+1).
reach(25,L + 1) :- reach(50, L), L <= 14, L>= 3, not reach(30, L+1), not reach(57, L+1), not reach(61, L+1).
reach(30,L + 1) :- reach(50, L), L <= 14, L>= 3, not reach(25, L+1), not reach(57, L+1), not reach(61, L+1).
reach(57,L + 1) :- reach(50, L), L <= 14, L>= 3, not reach(25, L+1), not reach(30, L+1), not reach(61, L+1).
reach(61,L + 1) :- reach(50, L), L <= 14, L>= 3, not reach(25, L+1), not reach(30, L+1), not reach(57, L+1).
reach(26,L + 1) :- reach(32, L), L <= 14, L>= 3, not reach(54, L+1), not reach(20, L+1), not reach(33, L+1).
reach(54,L + 1) :- reach(32, L), L <= 14, L>= 3, not reach(26, L+1), not reach(20, L+1), not reach(33, L+1).
reach(20,L + 1) :- reach(32, L), L <= 14, L>= 3, not reach(26, L+1), not reach(54, L+1), not reach(33, L+1).
reach(33,L + 1) :- reach(32, L), L <= 14, L>= 3, not reach(26, L+1), not reach(54, L+1), not reach(20, L+1).
reach(29,L + 1) :- reach(33, L), L <= 14, L>= 4, not reach(42, L+1), not reach(30, L+1), not reach(32, L+1).
reach(42,L + 1) :- reach(33, L), L <= 14, L>= 4, not reach(29, L+1), not reach(30, L+1), not reach(32, L+1).
reach(30,L + 1) :- reach(33, L), L <= 14, L>= 4, not reach(29, L+1), not reach(42, L+1), not reach(32, L+1).
reach(32,L + 1) :- reach(33, L), L <= 14, L>= 4, not reach(29, L+1), not reach(42, L+1), not reach(30, L+1).
reach(29,L + 1) :- reach(66, L), L <= 13, L>= 2, not reach(38, L+1), not reach(79, L+1), not reach(90, L+1), not reach(49, L+1).
reach(38,L + 1) :- reach(66, L), L <= 13, L>= 2, not reach(29, L+1), not reach(79, L+1), not reach(90, L+1), not reach(49, L+1).
reach(79,L + 1) :- reach(66, L), L <= 13, L>= 2, not reach(29, L+1), not reach(38, L+1), not reach(90, L+1), not reach(49, L+1).
reach(90,L + 1) :- reach(66, L), L <= 13, L>= 2, not reach(29, L+1), not reach(38, L+1), not reach(79, L+1), not reach(49, L+1).
reach(49,L + 1) :- reach(66, L), L <= 13, L>= 2, not reach(29, L+1), not reach(38, L+1), not reach(79, L+1), not reach(90, L+1).
reach(29,L + 1) :- reach(76, L), L <= 13, L>= 4, not reach(58, L+1), not reach(94, L+1).
reach(58,L + 1) :- reach(76, L), L <= 13, L>= 4, not reach(29, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(76, L), L <= 13, L>= 4, not reach(29, L+1), not reach(58, L+1).
reach(30,L + 1) :- reach(37, L), L <= 14, L>= 5, not reach(42, L+1), not reach(64, L+1), not reach(63, L+1).
reach(42,L + 1) :- reach(37, L), L <= 14, L>= 5, not reach(30, L+1), not reach(64, L+1), not reach(63, L+1).
reach(64,L + 1) :- reach(37, L), L <= 14, L>= 5, not reach(30, L+1), not reach(42, L+1), not reach(63, L+1).
reach(63,L + 1) :- reach(37, L), L <= 14, L>= 5, not reach(30, L+1), not reach(42, L+1), not reach(64, L+1).
reach(30,L + 1) :- reach(38, L), L <= 14, L>= 3, not reach(43, L+1), not reach(66, L+1), not reach(65, L+1).
reach(43,L + 1) :- reach(38, L), L <= 14, L>= 3, not reach(30, L+1), not reach(66, L+1), not reach(65, L+1).
reach(66,L + 1) :- reach(38, L), L <= 14, L>= 3, not reach(30, L+1), not reach(43, L+1), not reach(65, L+1).
reach(65,L + 1) :- reach(38, L), L <= 14, L>= 3, not reach(30, L+1), not reach(43, L+1), not reach(66, L+1).
reach(30,L + 1) :- reach(55, L), L <= 16, L>= 4, not reach(31, L+1), not reach(60, L+1), not reach(74, L+1).
reach(31,L + 1) :- reach(55, L), L <= 16, L>= 4, not reach(30, L+1), not reach(60, L+1), not reach(74, L+1).
reach(60,L + 1) :- reach(55, L), L <= 16, L>= 4, not reach(30, L+1), not reach(31, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(55, L), L <= 16, L>= 4, not reach(30, L+1), not reach(31, L+1), not reach(60, L+1).
reach(30,L + 1) :- reach(75, L), L <= 14, L>= 4, not reach(31, L+1), not reach(83, L+1).
reach(31,L + 1) :- reach(75, L), L <= 14, L>= 4, not reach(30, L+1), not reach(83, L+1).
reach(83,L + 1) :- reach(75, L), L <= 14, L>= 4, not reach(30, L+1), not reach(31, L+1).
reach(44,L + 1) :- reach(41, L), L <= 15, L>= 4, not reach(90, L+1), not reach(30, L+1), not reach(88, L+1), not reach(58, L+1).
reach(90,L + 1) :- reach(41, L), L <= 15, L>= 4, not reach(44, L+1), not reach(30, L+1), not reach(88, L+1), not reach(58, L+1).
reach(30,L + 1) :- reach(41, L), L <= 15, L>= 4, not reach(44, L+1), not reach(90, L+1), not reach(88, L+1), not reach(58, L+1).
reach(88,L + 1) :- reach(41, L), L <= 15, L>= 4, not reach(44, L+1), not reach(90, L+1), not reach(30, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(41, L), L <= 15, L>= 4, not reach(44, L+1), not reach(90, L+1), not reach(30, L+1), not reach(88, L+1).
reach(37,L + 1) :- reach(64, L), L <= 13, L>= 4, not reach(87, L+1), not reach(90, L+1).
reach(87,L + 1) :- reach(64, L), L <= 13, L>= 4, not reach(37, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(64, L), L <= 13, L>= 4, not reach(37, L+1), not reach(87, L+1).
reach(37,L + 1) :- reach(63, L), L <= 13, L>= 5, not reach(44, L+1), not reach(87, L+1).
reach(44,L + 1) :- reach(63, L), L <= 13, L>= 5, not reach(37, L+1), not reach(87, L+1).
reach(87,L + 1) :- reach(63, L), L <= 13, L>= 5, not reach(37, L+1), not reach(44, L+1).
reach(38,L + 1) :- reach(65, L), L <= 13, L>= 4, not reach(44, L+1), not reach(90, L+1).
reach(44,L + 1) :- reach(65, L), L <= 13, L>= 4, not reach(38, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(65, L), L <= 13, L>= 4, not reach(38, L+1), not reach(44, L+1).
reach(41,L + 1) :- reach(44, L), L <= 14, L>= 5, not reach(65, L+1), not reach(63, L+1), not reach(83, L+1).
reach(65,L + 1) :- reach(44, L), L <= 14, L>= 5, not reach(41, L+1), not reach(63, L+1), not reach(83, L+1).
reach(63,L + 1) :- reach(44, L), L <= 14, L>= 5, not reach(41, L+1), not reach(65, L+1), not reach(83, L+1).
reach(83,L + 1) :- reach(44, L), L <= 14, L>= 5, not reach(41, L+1), not reach(65, L+1), not reach(63, L+1).
reach(75,L + 1) :- reach(83, L), L <= 15, L>= 5, not reach(88, L+1), not reach(44, L+1).
reach(88,L + 1) :- reach(83, L), L <= 15, L>= 5, not reach(75, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(83, L), L <= 15, L>= 5, not reach(75, L+1), not reach(88, L+1).
reach(54,L + 1) :- reach(47, L), L <= 14, L>= 2, not reach(79, L+1), not reach(3, L+1).
reach(79,L + 1) :- reach(47, L), L <= 14, L>= 2, not reach(54, L+1), not reach(3, L+1).
reach(3,L + 1) :- reach(47, L), L <= 14, L>= 2, not reach(54, L+1), not reach(79, L+1).
reach(48,L + 1) :- reach(74, L), L <= 17, L>= 5, not reach(55, L+1), not reach(10, L+1), not reach(88, L+1).
reach(55,L + 1) :- reach(74, L), L <= 17, L>= 5, not reach(48, L+1), not reach(10, L+1), not reach(88, L+1).
reach(10,L + 1) :- reach(74, L), L <= 17, L>= 5, not reach(48, L+1), not reach(55, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(74, L), L <= 17, L>= 5, not reach(48, L+1), not reach(55, L+1), not reach(10, L+1).
reach(51,L + 1) :- reach(49, L), L <= 14, L>= 3, not reach(61, L+1), not reach(66, L+1).
reach(61,L + 1) :- reach(49, L), L <= 14, L>= 3, not reach(51, L+1), not reach(66, L+1).
reach(66,L + 1) :- reach(49, L), L <= 14, L>= 3, not reach(51, L+1), not reach(61, L+1).
reach(49,L + 1) :- reach(51, L), L <= 15, L>= 4, not reach(85, L+1), not reach(88, L+1).
reach(85,L + 1) :- reach(51, L), L <= 15, L>= 4, not reach(49, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(51, L), L <= 15, L>= 4, not reach(49, L+1), not reach(85, L+1).
reach(49,L + 1) :- reach(61, L), L <= 13, L>= 2, not reach(50, L+1), not reach(79, L+1).
reach(50,L + 1) :- reach(61, L), L <= 13, L>= 2, not reach(49, L+1), not reach(79, L+1).
reach(79,L + 1) :- reach(61, L), L <= 13, L>= 2, not reach(49, L+1), not reach(50, L+1).
reach(60,L + 1) :- reach(58, L), L <= 14, L>= 4, not reach(76, L+1), not reach(90, L+1), not reach(41, L+1).
reach(76,L + 1) :- reach(58, L), L <= 14, L>= 4, not reach(60, L+1), not reach(90, L+1), not reach(41, L+1).
reach(90,L + 1) :- reach(58, L), L <= 14, L>= 4, not reach(60, L+1), not reach(76, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(58, L), L <= 14, L>= 4, not reach(60, L+1), not reach(76, L+1), not reach(90, L+1).
reach(63,L + 1) :- reach(87, L), L <= 13, L>= 4, not reach(64, L+1), not reach(90, L+1).
reach(64,L + 1) :- reach(87, L), L <= 13, L>= 4, not reach(63, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(87, L), L <= 13, L>= 4, not reach(63, L+1), not reach(64, L+1).

