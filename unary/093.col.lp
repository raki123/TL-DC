reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12).
start(46).
goal(70).
:- reach(X, L), reach(X, L'), L != L'.
reach(2,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(7, L+1), not reach(5, L+1).
reach(7,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(2, L+1), not reach(5, L+1).
reach(5,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(2, L+1), not reach(7, L+1).
reach(49,L + 1) :- reach(5, L), L <= 7, L>= 4, not reach(34, L+1), not reach(48, L+1).
reach(34,L + 1) :- reach(5, L), L <= 7, L>= 4, not reach(49, L+1), not reach(48, L+1).
reach(48,L + 1) :- reach(5, L), L <= 7, L>= 4, not reach(49, L+1), not reach(34, L+1).
reach(34,L + 1) :- reach(2, L), L <= 9, L>= 5, not reach(35, L+1), not reach(90, L+1), not reach(96, L+1), not reach(62, L+1).
reach(35,L + 1) :- reach(2, L), L <= 9, L>= 5, not reach(34, L+1), not reach(90, L+1), not reach(96, L+1), not reach(62, L+1).
reach(90,L + 1) :- reach(2, L), L <= 9, L>= 5, not reach(34, L+1), not reach(35, L+1), not reach(96, L+1), not reach(62, L+1).
reach(96,L + 1) :- reach(2, L), L <= 9, L>= 5, not reach(34, L+1), not reach(35, L+1), not reach(90, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(2, L), L <= 9, L>= 5, not reach(34, L+1), not reach(35, L+1), not reach(90, L+1), not reach(96, L+1).
reach(2,L + 1) :- reach(35, L), L <= 8, L>= 5, not reach(3, L+1), not reach(88, L+1).
reach(3,L + 1) :- reach(35, L), L <= 8, L>= 5, not reach(2, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(35, L), L <= 8, L>= 5, not reach(2, L+1), not reach(3, L+1).
reach(2,L + 1) :- reach(90, L), L <= 10, L>= 5, not reach(71, L+1), not reach(78, L+1), not reach(96, L+1).
reach(71,L + 1) :- reach(90, L), L <= 10, L>= 5, not reach(2, L+1), not reach(78, L+1), not reach(96, L+1).
reach(78,L + 1) :- reach(90, L), L <= 10, L>= 5, not reach(2, L+1), not reach(71, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(90, L), L <= 10, L>= 5, not reach(2, L+1), not reach(71, L+1), not reach(78, L+1).
reach(2,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(23, L+1), not reach(26, L+1), not reach(32, L+1), not reach(53, L+1), not reach(57, L+1), not reach(61, L+1), not reach(90, L+1).
reach(23,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(26, L+1), not reach(32, L+1), not reach(53, L+1), not reach(57, L+1), not reach(61, L+1), not reach(90, L+1).
reach(26,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(23, L+1), not reach(32, L+1), not reach(53, L+1), not reach(57, L+1), not reach(61, L+1), not reach(90, L+1).
reach(32,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(23, L+1), not reach(26, L+1), not reach(53, L+1), not reach(57, L+1), not reach(61, L+1), not reach(90, L+1).
reach(53,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(23, L+1), not reach(26, L+1), not reach(32, L+1), not reach(57, L+1), not reach(61, L+1), not reach(90, L+1).
reach(57,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(23, L+1), not reach(26, L+1), not reach(32, L+1), not reach(53, L+1), not reach(61, L+1), not reach(90, L+1).
reach(61,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(23, L+1), not reach(26, L+1), not reach(32, L+1), not reach(53, L+1), not reach(57, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(96, L), L <= 9, L>= 4, not reach(2, L+1), not reach(23, L+1), not reach(26, L+1), not reach(32, L+1), not reach(53, L+1), not reach(57, L+1), not reach(61, L+1).
reach(35,L + 1) :- reach(3, L), L <= 7, L>= 5, not reach(42, L+1), not reach(18, L+1), not reach(88, L+1).
reach(42,L + 1) :- reach(3, L), L <= 7, L>= 5, not reach(35, L+1), not reach(18, L+1), not reach(88, L+1).
reach(18,L + 1) :- reach(3, L), L <= 7, L>= 5, not reach(35, L+1), not reach(42, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(3, L), L <= 7, L>= 5, not reach(35, L+1), not reach(42, L+1), not reach(18, L+1).
reach(3,L + 1) :- reach(42, L), L <= 6, L>= 6, not reach(98, L+1), not reach(60, L+1).
reach(98,L + 1) :- reach(42, L), L <= 6, L>= 6, not reach(3, L+1), not reach(60, L+1).
reach(60,L + 1) :- reach(42, L), L <= 6, L>= 6, not reach(3, L+1), not reach(98, L+1).
reach(3,L + 1) :- reach(18, L), L <= 7, L>= 5, not reach(7, L+1), not reach(57, L+1).
reach(7,L + 1) :- reach(18, L), L <= 7, L>= 5, not reach(3, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(18, L), L <= 7, L>= 5, not reach(3, L+1), not reach(7, L+1).
reach(3,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(35, L+1), not reach(39, L+1), not reach(61, L+1), not reach(66, L+1), not reach(74, L+1), not reach(75, L+1).
reach(35,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(3, L+1), not reach(39, L+1), not reach(61, L+1), not reach(66, L+1), not reach(74, L+1), not reach(75, L+1).
reach(39,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(3, L+1), not reach(35, L+1), not reach(61, L+1), not reach(66, L+1), not reach(74, L+1), not reach(75, L+1).
reach(61,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(3, L+1), not reach(35, L+1), not reach(39, L+1), not reach(66, L+1), not reach(74, L+1), not reach(75, L+1).
reach(66,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(3, L+1), not reach(35, L+1), not reach(39, L+1), not reach(61, L+1), not reach(74, L+1), not reach(75, L+1).
reach(74,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(3, L+1), not reach(35, L+1), not reach(39, L+1), not reach(61, L+1), not reach(66, L+1), not reach(75, L+1).
reach(75,L + 1) :- reach(88, L), L <= 7, L>= 4, not reach(3, L+1), not reach(35, L+1), not reach(39, L+1), not reach(61, L+1), not reach(66, L+1), not reach(74, L+1).
reach(39,L + 1) :- reach(66, L), L <= 8, L>= 5, not reach(53, L+1), not reach(62, L+1), not reach(88, L+1), not reach(26, L+1).
reach(53,L + 1) :- reach(66, L), L <= 8, L>= 5, not reach(39, L+1), not reach(62, L+1), not reach(88, L+1), not reach(26, L+1).
reach(62,L + 1) :- reach(66, L), L <= 8, L>= 5, not reach(39, L+1), not reach(53, L+1), not reach(88, L+1), not reach(26, L+1).
reach(88,L + 1) :- reach(66, L), L <= 8, L>= 5, not reach(39, L+1), not reach(53, L+1), not reach(62, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(66, L), L <= 8, L>= 5, not reach(39, L+1), not reach(53, L+1), not reach(62, L+1), not reach(88, L+1).
reach(48,L + 1) :- reach(26, L), L <= 8, L>= 5, not reach(96, L+1), not reach(66, L+1).
reach(96,L + 1) :- reach(26, L), L <= 8, L>= 5, not reach(48, L+1), not reach(66, L+1).
reach(66,L + 1) :- reach(26, L), L <= 8, L>= 5, not reach(48, L+1), not reach(96, L+1).
reach(5,L + 1) :- reach(49, L), L <= 7, L>= 3, not reach(47, L+1), not reach(55, L+1), not reach(92, L+1).
reach(47,L + 1) :- reach(49, L), L <= 7, L>= 3, not reach(5, L+1), not reach(55, L+1), not reach(92, L+1).
reach(55,L + 1) :- reach(49, L), L <= 7, L>= 3, not reach(5, L+1), not reach(47, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(49, L), L <= 7, L>= 3, not reach(5, L+1), not reach(47, L+1), not reach(55, L+1).
reach(37,L + 1) :- reach(6, L), L <= 6, L>= 5, not reach(45, L+1), not reach(60, L+1).
reach(45,L + 1) :- reach(6, L), L <= 6, L>= 5, not reach(37, L+1), not reach(60, L+1).
reach(60,L + 1) :- reach(6, L), L <= 6, L>= 5, not reach(37, L+1), not reach(45, L+1).
reach(6,L + 1) :- reach(37, L), L <= 7, L>= 5, not reach(57, L+1), not reach(55, L+1).
reach(57,L + 1) :- reach(37, L), L <= 7, L>= 5, not reach(6, L+1), not reach(55, L+1).
reach(55,L + 1) :- reach(37, L), L <= 7, L>= 5, not reach(6, L+1), not reach(57, L+1).
reach(6,L + 1) :- reach(45, L), L <= 7, L>= 4, not reach(28, L+1), not reach(61, L+1), not reach(7, L+1).
reach(28,L + 1) :- reach(45, L), L <= 7, L>= 4, not reach(6, L+1), not reach(61, L+1), not reach(7, L+1).
reach(61,L + 1) :- reach(45, L), L <= 7, L>= 4, not reach(6, L+1), not reach(28, L+1), not reach(7, L+1).
reach(7,L + 1) :- reach(45, L), L <= 7, L>= 4, not reach(6, L+1), not reach(28, L+1), not reach(61, L+1).
reach(6,L + 1) :- reach(60, L), L <= 7, L>= 6, not reach(21, L+1), not reach(53, L+1), not reach(42, L+1).
reach(21,L + 1) :- reach(60, L), L <= 7, L>= 6, not reach(6, L+1), not reach(53, L+1), not reach(42, L+1).
reach(53,L + 1) :- reach(60, L), L <= 7, L>= 6, not reach(6, L+1), not reach(21, L+1), not reach(42, L+1).
reach(42,L + 1) :- reach(60, L), L <= 7, L>= 6, not reach(6, L+1), not reach(21, L+1), not reach(53, L+1).
reach(34,L + 1) :- reach(7, L), L <= 7, L>= 5, not reach(18, L+1), not reach(8, L+1), not reach(23, L+1), not reach(45, L+1).
reach(18,L + 1) :- reach(7, L), L <= 7, L>= 5, not reach(34, L+1), not reach(8, L+1), not reach(23, L+1), not reach(45, L+1).
reach(8,L + 1) :- reach(7, L), L <= 7, L>= 5, not reach(34, L+1), not reach(18, L+1), not reach(23, L+1), not reach(45, L+1).
reach(23,L + 1) :- reach(7, L), L <= 7, L>= 5, not reach(34, L+1), not reach(18, L+1), not reach(8, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(7, L), L <= 7, L>= 5, not reach(34, L+1), not reach(18, L+1), not reach(8, L+1), not reach(23, L+1).
reach(7,L + 1) :- reach(8, L), L <= 8, L>= 6, not reach(9, L+1), not reach(11, L+1), not reach(12, L+1).
reach(9,L + 1) :- reach(8, L), L <= 8, L>= 6, not reach(7, L+1), not reach(11, L+1), not reach(12, L+1).
reach(11,L + 1) :- reach(8, L), L <= 8, L>= 6, not reach(7, L+1), not reach(9, L+1), not reach(12, L+1).
reach(12,L + 1) :- reach(8, L), L <= 8, L>= 6, not reach(7, L+1), not reach(9, L+1), not reach(11, L+1).
reach(8,L + 1) :- reach(9, L), L <= 7, L>= 5, not reach(36, L+1), not reach(98, L+1).
reach(36,L + 1) :- reach(9, L), L <= 7, L>= 5, not reach(8, L+1), not reach(98, L+1).
reach(98,L + 1) :- reach(9, L), L <= 7, L>= 5, not reach(8, L+1), not reach(36, L+1).
reach(8,L + 1) :- reach(11, L), L <= 7, L>= 5, not reach(57, L+1), not reach(14, L+1), not reach(98, L+1), not reach(23, L+1).
reach(57,L + 1) :- reach(11, L), L <= 7, L>= 5, not reach(8, L+1), not reach(14, L+1), not reach(98, L+1), not reach(23, L+1).
reach(14,L + 1) :- reach(11, L), L <= 7, L>= 5, not reach(8, L+1), not reach(57, L+1), not reach(98, L+1), not reach(23, L+1).
reach(98,L + 1) :- reach(11, L), L <= 7, L>= 5, not reach(8, L+1), not reach(57, L+1), not reach(14, L+1), not reach(23, L+1).
reach(23,L + 1) :- reach(11, L), L <= 7, L>= 5, not reach(8, L+1), not reach(57, L+1), not reach(14, L+1), not reach(98, L+1).
reach(9,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(12, L+1), not reach(74, L+1).
reach(12,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(9, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(9, L+1), not reach(12, L+1).
reach(9,L + 1) :- reach(98, L), L <= 6, L>= 6, not reach(42, L+1), not reach(11, L+1).
reach(42,L + 1) :- reach(98, L), L <= 6, L>= 6, not reach(9, L+1), not reach(11, L+1).
reach(11,L + 1) :- reach(98, L), L <= 6, L>= 6, not reach(9, L+1), not reach(42, L+1).
reach(11,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(18, L+1), not reach(30, L+1), not reach(37, L+1), not reach(48, L+1), not reach(62, L+1), not reach(96, L+1).
reach(18,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(11, L+1), not reach(30, L+1), not reach(37, L+1), not reach(48, L+1), not reach(62, L+1), not reach(96, L+1).
reach(30,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(11, L+1), not reach(18, L+1), not reach(37, L+1), not reach(48, L+1), not reach(62, L+1), not reach(96, L+1).
reach(37,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(11, L+1), not reach(18, L+1), not reach(30, L+1), not reach(48, L+1), not reach(62, L+1), not reach(96, L+1).
reach(48,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(11, L+1), not reach(18, L+1), not reach(30, L+1), not reach(37, L+1), not reach(62, L+1), not reach(96, L+1).
reach(62,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(11, L+1), not reach(18, L+1), not reach(30, L+1), not reach(37, L+1), not reach(48, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(57, L), L <= 8, L>= 4, not reach(11, L+1), not reach(18, L+1), not reach(30, L+1), not reach(37, L+1), not reach(48, L+1), not reach(62, L+1).
reach(64,L + 1) :- reach(12, L), L <= 9, L>= 5, not reach(8, L+1), not reach(36, L+1).
reach(8,L + 1) :- reach(12, L), L <= 9, L>= 5, not reach(64, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(12, L), L <= 9, L>= 5, not reach(64, L+1), not reach(8, L+1).
reach(12,L + 1) :- reach(64, L), L <= 10, L>= 5, not reach(71, L+1), not reach(84, L+1).
reach(71,L + 1) :- reach(64, L), L <= 10, L>= 5, not reach(12, L+1), not reach(84, L+1).
reach(84,L + 1) :- reach(64, L), L <= 10, L>= 5, not reach(12, L+1), not reach(71, L+1).
reach(55,L + 1) :- reach(14, L), L <= 7, L>= 5, not reach(39, L+1), not reach(75, L+1), not reach(11, L+1), not reach(28, L+1).
reach(39,L + 1) :- reach(14, L), L <= 7, L>= 5, not reach(55, L+1), not reach(75, L+1), not reach(11, L+1), not reach(28, L+1).
reach(75,L + 1) :- reach(14, L), L <= 7, L>= 5, not reach(55, L+1), not reach(39, L+1), not reach(11, L+1), not reach(28, L+1).
reach(11,L + 1) :- reach(14, L), L <= 7, L>= 5, not reach(55, L+1), not reach(39, L+1), not reach(75, L+1), not reach(28, L+1).
reach(28,L + 1) :- reach(14, L), L <= 7, L>= 5, not reach(55, L+1), not reach(39, L+1), not reach(75, L+1), not reach(11, L+1).
reach(14,L + 1) :- reach(55, L), L <= 6, L>= 4, not reach(37, L+1), not reach(49, L+1).
reach(37,L + 1) :- reach(55, L), L <= 6, L>= 4, not reach(14, L+1), not reach(49, L+1).
reach(49,L + 1) :- reach(55, L), L <= 6, L>= 4, not reach(14, L+1), not reach(37, L+1).
reach(14,L + 1) :- reach(39, L), L <= 7, L>= 4, not reach(66, L+1), not reach(76, L+1), not reach(88, L+1), not reach(69, L+1).
reach(66,L + 1) :- reach(39, L), L <= 7, L>= 4, not reach(14, L+1), not reach(76, L+1), not reach(88, L+1), not reach(69, L+1).
reach(76,L + 1) :- reach(39, L), L <= 7, L>= 4, not reach(14, L+1), not reach(66, L+1), not reach(88, L+1), not reach(69, L+1).
reach(88,L + 1) :- reach(39, L), L <= 7, L>= 4, not reach(14, L+1), not reach(66, L+1), not reach(76, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(39, L), L <= 7, L>= 4, not reach(14, L+1), not reach(66, L+1), not reach(76, L+1), not reach(88, L+1).
reach(14,L + 1) :- reach(75, L), L <= 8, L>= 5, not reach(19, L+1), not reach(32, L+1), not reach(88, L+1), not reach(62, L+1).
reach(19,L + 1) :- reach(75, L), L <= 8, L>= 5, not reach(14, L+1), not reach(32, L+1), not reach(88, L+1), not reach(62, L+1).
reach(32,L + 1) :- reach(75, L), L <= 8, L>= 5, not reach(14, L+1), not reach(19, L+1), not reach(88, L+1), not reach(62, L+1).
reach(88,L + 1) :- reach(75, L), L <= 8, L>= 5, not reach(14, L+1), not reach(19, L+1), not reach(32, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(75, L), L <= 8, L>= 5, not reach(14, L+1), not reach(19, L+1), not reach(32, L+1), not reach(88, L+1).
reach(61,L + 1) :- reach(16, L), L <= 7, L>= 2, not reach(74, L+1), not reach(29, L+1).
reach(74,L + 1) :- reach(16, L), L <= 7, L>= 2, not reach(61, L+1), not reach(29, L+1).
reach(29,L + 1) :- reach(16, L), L <= 7, L>= 2, not reach(61, L+1), not reach(74, L+1).
reach(16,L + 1) :- reach(61, L), L <= 8, L>= 3, not reach(45, L+1), not reach(56, L+1), not reach(69, L+1), not reach(88, L+1), not reach(96, L+1).
reach(45,L + 1) :- reach(61, L), L <= 8, L>= 3, not reach(16, L+1), not reach(56, L+1), not reach(69, L+1), not reach(88, L+1), not reach(96, L+1).
reach(56,L + 1) :- reach(61, L), L <= 8, L>= 3, not reach(16, L+1), not reach(45, L+1), not reach(69, L+1), not reach(88, L+1), not reach(96, L+1).
reach(69,L + 1) :- reach(61, L), L <= 8, L>= 3, not reach(16, L+1), not reach(45, L+1), not reach(56, L+1), not reach(88, L+1), not reach(96, L+1).
reach(88,L + 1) :- reach(61, L), L <= 8, L>= 3, not reach(16, L+1), not reach(45, L+1), not reach(56, L+1), not reach(69, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(61, L), L <= 8, L>= 3, not reach(16, L+1), not reach(45, L+1), not reach(56, L+1), not reach(69, L+1), not reach(88, L+1).
reach(16,L + 1) :- reach(74, L), L <= 8, L>= 3, not reach(84, L+1), not reach(88, L+1), not reach(92, L+1), not reach(36, L+1).
reach(84,L + 1) :- reach(74, L), L <= 8, L>= 3, not reach(16, L+1), not reach(88, L+1), not reach(92, L+1), not reach(36, L+1).
reach(88,L + 1) :- reach(74, L), L <= 8, L>= 3, not reach(16, L+1), not reach(84, L+1), not reach(92, L+1), not reach(36, L+1).
reach(92,L + 1) :- reach(74, L), L <= 8, L>= 3, not reach(16, L+1), not reach(84, L+1), not reach(88, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(74, L), L <= 8, L>= 3, not reach(16, L+1), not reach(84, L+1), not reach(88, L+1), not reach(92, L+1).
reach(51,L + 1) :- reach(19, L), L <= 9, L>= 6, not reach(75, L+1), not reach(78, L+1).
reach(75,L + 1) :- reach(19, L), L <= 9, L>= 6, not reach(51, L+1), not reach(78, L+1).
reach(78,L + 1) :- reach(19, L), L <= 9, L>= 6, not reach(51, L+1), not reach(75, L+1).
reach(19,L + 1) :- reach(51, L), L <= 10, L>= 6, not reach(62, L+1), not reach(85, L+1), not reach(71, L+1).
reach(62,L + 1) :- reach(51, L), L <= 10, L>= 6, not reach(19, L+1), not reach(85, L+1), not reach(71, L+1).
reach(85,L + 1) :- reach(51, L), L <= 10, L>= 6, not reach(19, L+1), not reach(62, L+1), not reach(71, L+1).
reach(71,L + 1) :- reach(51, L), L <= 10, L>= 6, not reach(19, L+1), not reach(62, L+1), not reach(85, L+1).
reach(19,L + 1) :- reach(78, L), L <= 9, L>= 6, not reach(43, L+1), not reach(90, L+1).
reach(43,L + 1) :- reach(78, L), L <= 9, L>= 6, not reach(19, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(78, L), L <= 9, L>= 6, not reach(19, L+1), not reach(43, L+1).
reach(52,L + 1) :- reach(56, L), L <= 7, L>= 4, not reach(61, L+1), not reach(43, L+1).
reach(61,L + 1) :- reach(56, L), L <= 7, L>= 4, not reach(52, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(56, L), L <= 7, L>= 4, not reach(52, L+1), not reach(61, L+1).
reach(60,L + 1) :- reach(21, L), L <= 7, L>= 5, not reach(76, L+1), not reach(94, L+1).
reach(76,L + 1) :- reach(21, L), L <= 7, L>= 5, not reach(60, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(21, L), L <= 7, L>= 5, not reach(60, L+1), not reach(76, L+1).
reach(21,L + 1) :- reach(76, L), L <= 8, L>= 5, not reach(39, L+1), not reach(92, L+1), not reach(84, L+1).
reach(39,L + 1) :- reach(76, L), L <= 8, L>= 5, not reach(21, L+1), not reach(92, L+1), not reach(84, L+1).
reach(92,L + 1) :- reach(76, L), L <= 8, L>= 5, not reach(21, L+1), not reach(39, L+1), not reach(84, L+1).
reach(84,L + 1) :- reach(76, L), L <= 8, L>= 5, not reach(21, L+1), not reach(39, L+1), not reach(92, L+1).
reach(21,L + 1) :- reach(94, L), L <= 8, L>= 4, not reach(30, L+1), not reach(62, L+1).
reach(30,L + 1) :- reach(94, L), L <= 8, L>= 4, not reach(21, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(94, L), L <= 8, L>= 4, not reach(21, L+1), not reach(30, L+1).
reach(23,L + 1) :- reach(53, L), L <= 8, L>= 5, not reach(66, L+1), not reach(96, L+1), not reach(60, L+1).
reach(66,L + 1) :- reach(53, L), L <= 8, L>= 5, not reach(23, L+1), not reach(96, L+1), not reach(60, L+1).
reach(96,L + 1) :- reach(53, L), L <= 8, L>= 5, not reach(23, L+1), not reach(66, L+1), not reach(60, L+1).
reach(60,L + 1) :- reach(53, L), L <= 8, L>= 5, not reach(23, L+1), not reach(66, L+1), not reach(96, L+1).
reach(53,L + 1) :- reach(23, L), L <= 8, L>= 5, not reach(96, L+1), not reach(7, L+1), not reach(11, L+1).
reach(96,L + 1) :- reach(23, L), L <= 8, L>= 5, not reach(53, L+1), not reach(7, L+1), not reach(11, L+1).
reach(7,L + 1) :- reach(23, L), L <= 8, L>= 5, not reach(53, L+1), not reach(96, L+1), not reach(11, L+1).
reach(11,L + 1) :- reach(23, L), L <= 8, L>= 5, not reach(53, L+1), not reach(96, L+1), not reach(7, L+1).
reach(48,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(51, L+1), not reach(57, L+1), not reach(66, L+1), not reach(94, L+1), not reach(2, L+1), not reach(75, L+1).
reach(51,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(48, L+1), not reach(57, L+1), not reach(66, L+1), not reach(94, L+1), not reach(2, L+1), not reach(75, L+1).
reach(57,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(48, L+1), not reach(51, L+1), not reach(66, L+1), not reach(94, L+1), not reach(2, L+1), not reach(75, L+1).
reach(66,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(48, L+1), not reach(51, L+1), not reach(57, L+1), not reach(94, L+1), not reach(2, L+1), not reach(75, L+1).
reach(94,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(48, L+1), not reach(51, L+1), not reach(57, L+1), not reach(66, L+1), not reach(2, L+1), not reach(75, L+1).
reach(2,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(48, L+1), not reach(51, L+1), not reach(57, L+1), not reach(66, L+1), not reach(94, L+1), not reach(75, L+1).
reach(75,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(48, L+1), not reach(51, L+1), not reach(57, L+1), not reach(66, L+1), not reach(94, L+1), not reach(2, L+1).
reach(26,L + 1) :- reach(48, L), L <= 8, L>= 5, not reach(28, L+1), not reach(57, L+1), not reach(62, L+1), not reach(92, L+1), not reach(5, L+1).
reach(28,L + 1) :- reach(48, L), L <= 8, L>= 5, not reach(26, L+1), not reach(57, L+1), not reach(62, L+1), not reach(92, L+1), not reach(5, L+1).
reach(57,L + 1) :- reach(48, L), L <= 8, L>= 5, not reach(26, L+1), not reach(28, L+1), not reach(62, L+1), not reach(92, L+1), not reach(5, L+1).
reach(62,L + 1) :- reach(48, L), L <= 8, L>= 5, not reach(26, L+1), not reach(28, L+1), not reach(57, L+1), not reach(92, L+1), not reach(5, L+1).
reach(92,L + 1) :- reach(48, L), L <= 8, L>= 5, not reach(26, L+1), not reach(28, L+1), not reach(57, L+1), not reach(62, L+1), not reach(5, L+1).
reach(5,L + 1) :- reach(48, L), L <= 8, L>= 5, not reach(26, L+1), not reach(28, L+1), not reach(57, L+1), not reach(62, L+1), not reach(92, L+1).
reach(45,L + 1) :- reach(28, L), L <= 7, L>= 5, not reach(48, L+1), not reach(14, L+1).
reach(48,L + 1) :- reach(28, L), L <= 7, L>= 5, not reach(45, L+1), not reach(14, L+1).
reach(14,L + 1) :- reach(28, L), L <= 7, L>= 5, not reach(45, L+1), not reach(48, L+1).
reach(47,L + 1) :- reach(29, L), L <= 6, L>= 1, not reach(46, L+1), not reach(16, L+1).
reach(16,L + 1) :- reach(29, L), L <= 6, L>= 1, not reach(46, L+1), not reach(47, L+1).
reach(29,L + 1) :- reach(46, L), L <= 5, L>= 0.
reach(29,L + 1) :- reach(47, L), L <= 6, L>= 2, not reach(30, L+1), not reach(69, L+1), not reach(49, L+1).
reach(30,L + 1) :- reach(47, L), L <= 6, L>= 2, not reach(29, L+1), not reach(69, L+1), not reach(49, L+1).
reach(69,L + 1) :- reach(47, L), L <= 6, L>= 2, not reach(29, L+1), not reach(30, L+1), not reach(49, L+1).
reach(49,L + 1) :- reach(47, L), L <= 6, L>= 2, not reach(29, L+1), not reach(30, L+1), not reach(69, L+1).
reach(47,L + 1) :- reach(30, L), L <= 7, L>= 3, not reach(57, L+1), not reach(94, L+1).
reach(57,L + 1) :- reach(30, L), L <= 7, L>= 3, not reach(47, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(30, L), L <= 7, L>= 3, not reach(47, L+1), not reach(57, L+1).
reach(75,L + 1) :- reach(32, L), L <= 8, L>= 5, not reach(85, L+1), not reach(96, L+1), not reach(52, L+1).
reach(85,L + 1) :- reach(32, L), L <= 8, L>= 5, not reach(75, L+1), not reach(96, L+1), not reach(52, L+1).
reach(96,L + 1) :- reach(32, L), L <= 8, L>= 5, not reach(75, L+1), not reach(85, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(32, L), L <= 8, L>= 5, not reach(75, L+1), not reach(85, L+1), not reach(96, L+1).
reach(32,L + 1) :- reach(85, L), L <= 9, L>= 5, not reach(51, L+1), not reach(92, L+1).
reach(51,L + 1) :- reach(85, L), L <= 9, L>= 5, not reach(32, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(85, L), L <= 9, L>= 5, not reach(32, L+1), not reach(51, L+1).
reach(56,L + 1) :- reach(52, L), L <= 7, L>= 5, not reach(92, L+1), not reach(32, L+1), not reach(43, L+1).
reach(92,L + 1) :- reach(52, L), L <= 7, L>= 5, not reach(56, L+1), not reach(32, L+1), not reach(43, L+1).
reach(32,L + 1) :- reach(52, L), L <= 7, L>= 5, not reach(56, L+1), not reach(92, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(52, L), L <= 7, L>= 5, not reach(56, L+1), not reach(92, L+1), not reach(32, L+1).
reach(78,L + 1) :- reach(43, L), L <= 8, L>= 5, not reach(56, L+1), not reach(52, L+1).
reach(56,L + 1) :- reach(43, L), L <= 8, L>= 5, not reach(78, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(43, L), L <= 8, L>= 5, not reach(78, L+1), not reach(56, L+1).
reach(39,L + 1) :- reach(69, L), L <= 7, L>= 3, not reach(47, L+1), not reach(61, L+1).
reach(47,L + 1) :- reach(69, L), L <= 7, L>= 3, not reach(39, L+1), not reach(61, L+1).
reach(61,L + 1) :- reach(69, L), L <= 7, L>= 3, not reach(39, L+1), not reach(47, L+1).
reach(48,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(49, L+1), not reach(52, L+1), not reach(74, L+1), not reach(76, L+1), not reach(85, L+1).
reach(49,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(48, L+1), not reach(52, L+1), not reach(74, L+1), not reach(76, L+1), not reach(85, L+1).
reach(52,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(48, L+1), not reach(49, L+1), not reach(74, L+1), not reach(76, L+1), not reach(85, L+1).
reach(74,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(48, L+1), not reach(49, L+1), not reach(52, L+1), not reach(76, L+1), not reach(85, L+1).
reach(76,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(48, L+1), not reach(49, L+1), not reach(52, L+1), not reach(74, L+1), not reach(85, L+1).
reach(85,L + 1) :- reach(92, L), L <= 8, L>= 4, not reach(48, L+1), not reach(49, L+1), not reach(52, L+1), not reach(74, L+1), not reach(76, L+1).
reach(51,L + 1) :- reach(71, L), L <= 11, L>= 6, not reach(64, L+1), not reach(90, L+1), not reach(70, L+1).
reach(64,L + 1) :- reach(71, L), L <= 11, L>= 6, not reach(51, L+1), not reach(90, L+1), not reach(70, L+1).
reach(90,L + 1) :- reach(71, L), L <= 11, L>= 6, not reach(51, L+1), not reach(64, L+1), not reach(70, L+1).
reach(70,L + 1) :- reach(71, L), L <= 11, L>= 6, not reach(51, L+1), not reach(64, L+1), not reach(90, L+1).
reach(64,L + 1) :- reach(84, L), L <= 9, L>= 4, not reach(74, L+1), not reach(76, L+1).
reach(74,L + 1) :- reach(84, L), L <= 9, L>= 4, not reach(64, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(84, L), L <= 9, L>= 4, not reach(64, L+1), not reach(74, L+1).

