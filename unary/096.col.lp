reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17).
start(7).
goal(45).
:- reach(X, L), reach(X, L'), L != L'.
reach(3,L + 1) :- reach(2, L), L <= 13, L>= 4, not reach(50, L+1), not reach(91, L+1).
reach(50,L + 1) :- reach(2, L), L <= 13, L>= 4, not reach(3, L+1), not reach(91, L+1).
reach(91,L + 1) :- reach(2, L), L <= 13, L>= 4, not reach(3, L+1), not reach(50, L+1).
reach(3,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(8, L+1), not reach(30, L+1), not reach(44, L+1), not reach(71, L+1), not reach(87, L+1), not reach(2, L+1), not reach(77, L+1).
reach(8,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(30, L+1), not reach(44, L+1), not reach(71, L+1), not reach(87, L+1), not reach(2, L+1), not reach(77, L+1).
reach(30,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(8, L+1), not reach(44, L+1), not reach(71, L+1), not reach(87, L+1), not reach(2, L+1), not reach(77, L+1).
reach(44,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(8, L+1), not reach(30, L+1), not reach(71, L+1), not reach(87, L+1), not reach(2, L+1), not reach(77, L+1).
reach(71,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(8, L+1), not reach(30, L+1), not reach(44, L+1), not reach(87, L+1), not reach(2, L+1), not reach(77, L+1).
reach(87,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(8, L+1), not reach(30, L+1), not reach(44, L+1), not reach(71, L+1), not reach(2, L+1), not reach(77, L+1).
reach(2,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(8, L+1), not reach(30, L+1), not reach(44, L+1), not reach(71, L+1), not reach(87, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(91, L), L <= 14, L>= 4, not reach(3, L+1), not reach(8, L+1), not reach(30, L+1), not reach(44, L+1), not reach(71, L+1), not reach(87, L+1), not reach(2, L+1).
reach(2,L + 1) :- reach(3, L), L <= 13, L>= 5, not reach(56, L+1), not reach(91, L+1), not reach(5, L+1).
reach(56,L + 1) :- reach(3, L), L <= 13, L>= 5, not reach(2, L+1), not reach(91, L+1), not reach(5, L+1).
reach(91,L + 1) :- reach(3, L), L <= 13, L>= 5, not reach(2, L+1), not reach(56, L+1), not reach(5, L+1).
reach(5,L + 1) :- reach(3, L), L <= 13, L>= 5, not reach(2, L+1), not reach(56, L+1), not reach(91, L+1).
reach(2,L + 1) :- reach(50, L), L <= 13, L>= 3, not reach(51, L+1), not reach(52, L+1), not reach(77, L+1).
reach(51,L + 1) :- reach(50, L), L <= 13, L>= 3, not reach(2, L+1), not reach(52, L+1), not reach(77, L+1).
reach(52,L + 1) :- reach(50, L), L <= 13, L>= 3, not reach(2, L+1), not reach(51, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(50, L), L <= 13, L>= 3, not reach(2, L+1), not reach(51, L+1), not reach(52, L+1).
reach(3,L + 1) :- reach(56, L), L <= 13, L>= 4, not reach(34, L+1), not reach(43, L+1).
reach(34,L + 1) :- reach(56, L), L <= 13, L>= 4, not reach(3, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(56, L), L <= 13, L>= 4, not reach(3, L+1), not reach(34, L+1).
reach(6,L + 1) :- reach(4, L), L <= 14, L>= 2, not reach(71, L+1), not reach(73, L+1), not reach(48, L+1).
reach(71,L + 1) :- reach(4, L), L <= 14, L>= 2, not reach(6, L+1), not reach(73, L+1), not reach(48, L+1).
reach(73,L + 1) :- reach(4, L), L <= 14, L>= 2, not reach(6, L+1), not reach(71, L+1), not reach(48, L+1).
reach(48,L + 1) :- reach(4, L), L <= 14, L>= 2, not reach(6, L+1), not reach(71, L+1), not reach(73, L+1).
reach(4,L + 1) :- reach(6, L), L <= 13, L>= 1, not reach(7, L+1), not reach(51, L+1), not reach(81, L+1), not reach(29, L+1).
reach(51,L + 1) :- reach(6, L), L <= 13, L>= 1, not reach(4, L+1), not reach(7, L+1), not reach(81, L+1), not reach(29, L+1).
reach(81,L + 1) :- reach(6, L), L <= 13, L>= 1, not reach(4, L+1), not reach(7, L+1), not reach(51, L+1), not reach(29, L+1).
reach(29,L + 1) :- reach(6, L), L <= 13, L>= 1, not reach(4, L+1), not reach(7, L+1), not reach(51, L+1), not reach(81, L+1).
reach(4,L + 1) :- reach(71, L), L <= 13, L>= 3, not reach(22, L+1), not reach(39, L+1), not reach(82, L+1), not reach(91, L+1), not reach(8, L+1).
reach(22,L + 1) :- reach(71, L), L <= 13, L>= 3, not reach(4, L+1), not reach(39, L+1), not reach(82, L+1), not reach(91, L+1), not reach(8, L+1).
reach(39,L + 1) :- reach(71, L), L <= 13, L>= 3, not reach(4, L+1), not reach(22, L+1), not reach(82, L+1), not reach(91, L+1), not reach(8, L+1).
reach(82,L + 1) :- reach(71, L), L <= 13, L>= 3, not reach(4, L+1), not reach(22, L+1), not reach(39, L+1), not reach(91, L+1), not reach(8, L+1).
reach(91,L + 1) :- reach(71, L), L <= 13, L>= 3, not reach(4, L+1), not reach(22, L+1), not reach(39, L+1), not reach(82, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(71, L), L <= 13, L>= 3, not reach(4, L+1), not reach(22, L+1), not reach(39, L+1), not reach(82, L+1), not reach(91, L+1).
reach(4,L + 1) :- reach(73, L), L <= 13, L>= 3, not reach(46, L+1), not reach(89, L+1).
reach(46,L + 1) :- reach(73, L), L <= 13, L>= 3, not reach(4, L+1), not reach(89, L+1).
reach(89,L + 1) :- reach(73, L), L <= 13, L>= 3, not reach(4, L+1), not reach(46, L+1).
reach(25,L + 1) :- reach(5, L), L <= 12, L>= 4, not reach(72, L+1), not reach(3, L+1).
reach(72,L + 1) :- reach(5, L), L <= 12, L>= 4, not reach(25, L+1), not reach(3, L+1).
reach(3,L + 1) :- reach(5, L), L <= 12, L>= 4, not reach(25, L+1), not reach(72, L+1).
reach(5,L + 1) :- reach(25, L), L <= 13, L>= 3, not reach(24, L+1), not reach(81, L+1), not reach(82, L+1), not reach(39, L+1).
reach(24,L + 1) :- reach(25, L), L <= 13, L>= 3, not reach(5, L+1), not reach(81, L+1), not reach(82, L+1), not reach(39, L+1).
reach(81,L + 1) :- reach(25, L), L <= 13, L>= 3, not reach(5, L+1), not reach(24, L+1), not reach(82, L+1), not reach(39, L+1).
reach(82,L + 1) :- reach(25, L), L <= 13, L>= 3, not reach(5, L+1), not reach(24, L+1), not reach(81, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(25, L), L <= 13, L>= 3, not reach(5, L+1), not reach(24, L+1), not reach(81, L+1), not reach(82, L+1).
reach(5,L + 1) :- reach(72, L), L <= 12, L>= 3, not reach(22, L+1), not reach(69, L+1), not reach(87, L+1).
reach(22,L + 1) :- reach(72, L), L <= 12, L>= 3, not reach(5, L+1), not reach(69, L+1), not reach(87, L+1).
reach(69,L + 1) :- reach(72, L), L <= 12, L>= 3, not reach(5, L+1), not reach(22, L+1), not reach(87, L+1).
reach(87,L + 1) :- reach(72, L), L <= 12, L>= 3, not reach(5, L+1), not reach(22, L+1), not reach(69, L+1).
reach(6,L + 1) :- reach(7, L), L <= 12, L>= 0, not reach(13, L+1), not reach(31, L+1).
reach(13,L + 1) :- reach(7, L), L <= 12, L>= 0, not reach(6, L+1), not reach(31, L+1).
reach(31,L + 1) :- reach(7, L), L <= 12, L>= 0, not reach(6, L+1), not reach(13, L+1).
reach(6,L + 1) :- reach(51, L), L <= 12, L>= 2, not reach(49, L+1), not reach(50, L+1), not reach(87, L+1).
reach(49,L + 1) :- reach(51, L), L <= 12, L>= 2, not reach(6, L+1), not reach(50, L+1), not reach(87, L+1).
reach(50,L + 1) :- reach(51, L), L <= 12, L>= 2, not reach(6, L+1), not reach(49, L+1), not reach(87, L+1).
reach(87,L + 1) :- reach(51, L), L <= 12, L>= 2, not reach(6, L+1), not reach(49, L+1), not reach(50, L+1).
reach(14,L + 1) :- reach(8, L), L <= 13, L>= 3, not reach(77, L+1), not reach(91, L+1), not reach(71, L+1).
reach(77,L + 1) :- reach(8, L), L <= 13, L>= 3, not reach(14, L+1), not reach(91, L+1), not reach(71, L+1).
reach(91,L + 1) :- reach(8, L), L <= 13, L>= 3, not reach(14, L+1), not reach(77, L+1), not reach(71, L+1).
reach(71,L + 1) :- reach(8, L), L <= 13, L>= 3, not reach(14, L+1), not reach(77, L+1), not reach(91, L+1).
reach(8,L + 1) :- reach(14, L), L <= 12, L>= 2, not reach(13, L+1), not reach(26, L+1).
reach(13,L + 1) :- reach(14, L), L <= 12, L>= 2, not reach(8, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(14, L), L <= 12, L>= 2, not reach(8, L+1), not reach(13, L+1).
reach(8,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(36, L+1), not reach(47, L+1), not reach(62, L+1), not reach(67, L+1), not reach(91, L+1), not reach(50, L+1), not reach(10, L+1).
reach(36,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(47, L+1), not reach(62, L+1), not reach(67, L+1), not reach(91, L+1), not reach(50, L+1), not reach(10, L+1).
reach(47,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(36, L+1), not reach(62, L+1), not reach(67, L+1), not reach(91, L+1), not reach(50, L+1), not reach(10, L+1).
reach(62,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(36, L+1), not reach(47, L+1), not reach(67, L+1), not reach(91, L+1), not reach(50, L+1), not reach(10, L+1).
reach(67,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(36, L+1), not reach(47, L+1), not reach(62, L+1), not reach(91, L+1), not reach(50, L+1), not reach(10, L+1).
reach(91,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(36, L+1), not reach(47, L+1), not reach(62, L+1), not reach(67, L+1), not reach(50, L+1), not reach(10, L+1).
reach(50,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(36, L+1), not reach(47, L+1), not reach(62, L+1), not reach(67, L+1), not reach(91, L+1), not reach(10, L+1).
reach(10,L + 1) :- reach(77, L), L <= 14, L>= 4, not reach(8, L+1), not reach(36, L+1), not reach(47, L+1), not reach(62, L+1), not reach(67, L+1), not reach(91, L+1), not reach(50, L+1).
reach(15,L + 1) :- reach(10, L), L <= 13, L>= 4, not reach(59, L+1), not reach(77, L+1).
reach(59,L + 1) :- reach(10, L), L <= 13, L>= 4, not reach(15, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(10, L), L <= 13, L>= 4, not reach(15, L+1), not reach(59, L+1).
reach(10,L + 1) :- reach(15, L), L <= 13, L>= 3, not reach(66, L+1), not reach(29, L+1).
reach(66,L + 1) :- reach(15, L), L <= 13, L>= 3, not reach(10, L+1), not reach(29, L+1).
reach(29,L + 1) :- reach(15, L), L <= 13, L>= 3, not reach(10, L+1), not reach(66, L+1).
reach(10,L + 1) :- reach(59, L), L <= 12, L>= 3, not reach(22, L+1), not reach(52, L+1).
reach(22,L + 1) :- reach(59, L), L <= 12, L>= 3, not reach(10, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(59, L), L <= 12, L>= 3, not reach(10, L+1), not reach(22, L+1).
reach(52,L + 1) :- reach(13, L), L <= 12, L>= 1, not reach(14, L+1), not reach(7, L+1), not reach(69, L+1).
reach(14,L + 1) :- reach(13, L), L <= 12, L>= 1, not reach(52, L+1), not reach(7, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(13, L), L <= 12, L>= 1, not reach(52, L+1), not reach(14, L+1), not reach(7, L+1).
reach(13,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(50, L+1), not reach(66, L+1), not reach(59, L+1), not reach(81, L+1), not reach(89, L+1), not reach(67, L+1).
reach(50,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(13, L+1), not reach(66, L+1), not reach(59, L+1), not reach(81, L+1), not reach(89, L+1), not reach(67, L+1).
reach(66,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(13, L+1), not reach(50, L+1), not reach(59, L+1), not reach(81, L+1), not reach(89, L+1), not reach(67, L+1).
reach(59,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(13, L+1), not reach(50, L+1), not reach(66, L+1), not reach(81, L+1), not reach(89, L+1), not reach(67, L+1).
reach(81,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(13, L+1), not reach(50, L+1), not reach(66, L+1), not reach(59, L+1), not reach(89, L+1), not reach(67, L+1).
reach(89,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(13, L+1), not reach(50, L+1), not reach(66, L+1), not reach(59, L+1), not reach(81, L+1), not reach(67, L+1).
reach(67,L + 1) :- reach(52, L), L <= 13, L>= 2, not reach(13, L+1), not reach(50, L+1), not reach(66, L+1), not reach(59, L+1), not reach(81, L+1), not reach(89, L+1).
reach(15,L + 1) :- reach(66, L), L <= 14, L>= 3, not reach(48, L+1), not reach(52, L+1), not reach(81, L+1), not reach(39, L+1).
reach(48,L + 1) :- reach(66, L), L <= 14, L>= 3, not reach(15, L+1), not reach(52, L+1), not reach(81, L+1), not reach(39, L+1).
reach(52,L + 1) :- reach(66, L), L <= 14, L>= 3, not reach(15, L+1), not reach(48, L+1), not reach(81, L+1), not reach(39, L+1).
reach(81,L + 1) :- reach(66, L), L <= 14, L>= 3, not reach(15, L+1), not reach(48, L+1), not reach(52, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(66, L), L <= 14, L>= 3, not reach(15, L+1), not reach(48, L+1), not reach(52, L+1), not reach(81, L+1).
reach(53,L + 1) :- reach(18, L), L <= 12, L>= 3, not reach(32, L+1), not reach(81, L+1).
reach(32,L + 1) :- reach(18, L), L <= 12, L>= 3, not reach(53, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(18, L), L <= 12, L>= 3, not reach(53, L+1), not reach(32, L+1).
reach(18,L + 1) :- reach(53, L), L <= 12, L>= 4, not reach(32, L+1), not reach(93, L+1).
reach(32,L + 1) :- reach(53, L), L <= 12, L>= 4, not reach(18, L+1), not reach(93, L+1).
reach(93,L + 1) :- reach(53, L), L <= 12, L>= 4, not reach(18, L+1), not reach(32, L+1).
reach(18,L + 1) :- reach(32, L), L <= 13, L>= 4, not reach(40, L+1), not reach(53, L+1).
reach(40,L + 1) :- reach(32, L), L <= 13, L>= 4, not reach(18, L+1), not reach(53, L+1).
reach(53,L + 1) :- reach(32, L), L <= 13, L>= 4, not reach(18, L+1), not reach(40, L+1).
reach(18,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(25, L+1), not reach(36, L+1), not reach(40, L+1), not reach(52, L+1), not reach(66, L+1), not reach(69, L+1), not reach(89, L+1), not reach(6, L+1).
reach(25,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(36, L+1), not reach(40, L+1), not reach(52, L+1), not reach(66, L+1), not reach(69, L+1), not reach(89, L+1), not reach(6, L+1).
reach(36,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(40, L+1), not reach(52, L+1), not reach(66, L+1), not reach(69, L+1), not reach(89, L+1), not reach(6, L+1).
reach(40,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(36, L+1), not reach(52, L+1), not reach(66, L+1), not reach(69, L+1), not reach(89, L+1), not reach(6, L+1).
reach(52,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(36, L+1), not reach(40, L+1), not reach(66, L+1), not reach(69, L+1), not reach(89, L+1), not reach(6, L+1).
reach(66,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(36, L+1), not reach(40, L+1), not reach(52, L+1), not reach(69, L+1), not reach(89, L+1), not reach(6, L+1).
reach(69,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(36, L+1), not reach(40, L+1), not reach(52, L+1), not reach(66, L+1), not reach(89, L+1), not reach(6, L+1).
reach(89,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(36, L+1), not reach(40, L+1), not reach(52, L+1), not reach(66, L+1), not reach(69, L+1), not reach(6, L+1).
reach(6,L + 1) :- reach(81, L), L <= 13, L>= 2, not reach(18, L+1), not reach(25, L+1), not reach(36, L+1), not reach(40, L+1), not reach(52, L+1), not reach(66, L+1), not reach(69, L+1), not reach(89, L+1).
reach(71,L + 1) :- reach(82, L), L <= 14, L>= 4, not reach(25, L+1), not reach(35, L+1), not reach(57, L+1).
reach(25,L + 1) :- reach(82, L), L <= 14, L>= 4, not reach(71, L+1), not reach(35, L+1), not reach(57, L+1).
reach(35,L + 1) :- reach(82, L), L <= 14, L>= 4, not reach(71, L+1), not reach(25, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(82, L), L <= 14, L>= 4, not reach(71, L+1), not reach(25, L+1), not reach(35, L+1).
reach(26,L + 1) :- reach(31, L), L <= 12, L>= 1, not reach(39, L+1), not reach(49, L+1), not reach(7, L+1).
reach(39,L + 1) :- reach(31, L), L <= 12, L>= 1, not reach(26, L+1), not reach(49, L+1), not reach(7, L+1).
reach(49,L + 1) :- reach(31, L), L <= 12, L>= 1, not reach(26, L+1), not reach(39, L+1), not reach(7, L+1).
reach(59,L + 1) :- reach(22, L), L <= 13, L>= 4, not reach(71, L+1), not reach(72, L+1), not reach(80, L+1), not reach(24, L+1).
reach(71,L + 1) :- reach(22, L), L <= 13, L>= 4, not reach(59, L+1), not reach(72, L+1), not reach(80, L+1), not reach(24, L+1).
reach(72,L + 1) :- reach(22, L), L <= 13, L>= 4, not reach(59, L+1), not reach(71, L+1), not reach(80, L+1), not reach(24, L+1).
reach(80,L + 1) :- reach(22, L), L <= 13, L>= 4, not reach(59, L+1), not reach(71, L+1), not reach(72, L+1), not reach(24, L+1).
reach(24,L + 1) :- reach(22, L), L <= 13, L>= 4, not reach(59, L+1), not reach(71, L+1), not reach(72, L+1), not reach(80, L+1).
reach(22,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(41, L+1), not reach(43, L+1), not reach(47, L+1), not reach(48, L+1), not reach(49, L+1), not reach(57, L+1).
reach(41,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(22, L+1), not reach(43, L+1), not reach(47, L+1), not reach(48, L+1), not reach(49, L+1), not reach(57, L+1).
reach(43,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(22, L+1), not reach(41, L+1), not reach(47, L+1), not reach(48, L+1), not reach(49, L+1), not reach(57, L+1).
reach(47,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(22, L+1), not reach(41, L+1), not reach(43, L+1), not reach(48, L+1), not reach(49, L+1), not reach(57, L+1).
reach(48,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(22, L+1), not reach(41, L+1), not reach(43, L+1), not reach(47, L+1), not reach(49, L+1), not reach(57, L+1).
reach(49,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(22, L+1), not reach(41, L+1), not reach(43, L+1), not reach(47, L+1), not reach(48, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(80, L), L <= 14, L>= 3, not reach(22, L+1), not reach(41, L+1), not reach(43, L+1), not reach(47, L+1), not reach(48, L+1), not reach(49, L+1).
reach(25,L + 1) :- reach(24, L), L <= 13, L>= 4, not reach(40, L+1), not reach(22, L+1).
reach(40,L + 1) :- reach(24, L), L <= 13, L>= 4, not reach(25, L+1), not reach(22, L+1).
reach(22,L + 1) :- reach(24, L), L <= 13, L>= 4, not reach(25, L+1), not reach(40, L+1).
reach(24,L + 1) :- reach(40, L), L <= 14, L>= 3, not reach(32, L+1), not reach(57, L+1), not reach(81, L+1).
reach(32,L + 1) :- reach(40, L), L <= 14, L>= 3, not reach(24, L+1), not reach(57, L+1), not reach(81, L+1).
reach(57,L + 1) :- reach(40, L), L <= 14, L>= 3, not reach(24, L+1), not reach(32, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(40, L), L <= 14, L>= 3, not reach(24, L+1), not reach(32, L+1), not reach(57, L+1).
reach(31,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(34, L+1), not reach(35, L+1), not reach(14, L+1), not reach(49, L+1).
reach(34,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(31, L+1), not reach(35, L+1), not reach(14, L+1), not reach(49, L+1).
reach(35,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(31, L+1), not reach(34, L+1), not reach(14, L+1), not reach(49, L+1).
reach(14,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(31, L+1), not reach(34, L+1), not reach(35, L+1), not reach(49, L+1).
reach(49,L + 1) :- reach(26, L), L <= 13, L>= 2, not reach(31, L+1), not reach(34, L+1), not reach(35, L+1), not reach(14, L+1).
reach(26,L + 1) :- reach(34, L), L <= 14, L>= 3, not reach(41, L+1), not reach(48, L+1), not reach(56, L+1).
reach(41,L + 1) :- reach(34, L), L <= 14, L>= 3, not reach(26, L+1), not reach(48, L+1), not reach(56, L+1).
reach(48,L + 1) :- reach(34, L), L <= 14, L>= 3, not reach(26, L+1), not reach(41, L+1), not reach(56, L+1).
reach(56,L + 1) :- reach(34, L), L <= 14, L>= 3, not reach(26, L+1), not reach(41, L+1), not reach(48, L+1).
reach(26,L + 1) :- reach(35, L), L <= 14, L>= 3, not reach(47, L+1), not reach(82, L+1).
reach(47,L + 1) :- reach(35, L), L <= 14, L>= 3, not reach(26, L+1), not reach(82, L+1).
reach(82,L + 1) :- reach(35, L), L <= 14, L>= 3, not reach(26, L+1), not reach(47, L+1).
reach(44,L + 1) :- reach(29, L), L <= 14, L>= 2, not reach(6, L+1), not reach(15, L+1).
reach(6,L + 1) :- reach(29, L), L <= 14, L>= 2, not reach(44, L+1), not reach(15, L+1).
reach(15,L + 1) :- reach(29, L), L <= 14, L>= 2, not reach(44, L+1), not reach(6, L+1).
reach(29,L + 1) :- reach(44, L), L <= 15, L>= 3, not reach(91, L+1), not reach(88, L+1).
reach(91,L + 1) :- reach(44, L), L <= 15, L>= 3, not reach(29, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(44, L), L <= 15, L>= 3, not reach(29, L+1), not reach(91, L+1).
reach(43,L + 1) :- reach(30, L), L <= 13, L>= 5, not reach(41, L+1), not reach(91, L+1).
reach(41,L + 1) :- reach(30, L), L <= 13, L>= 5, not reach(43, L+1), not reach(91, L+1).
reach(91,L + 1) :- reach(30, L), L <= 13, L>= 5, not reach(43, L+1), not reach(41, L+1).
reach(30,L + 1) :- reach(43, L), L <= 13, L>= 4, not reach(36, L+1), not reach(56, L+1), not reach(80, L+1), not reach(93, L+1).
reach(36,L + 1) :- reach(43, L), L <= 13, L>= 4, not reach(30, L+1), not reach(56, L+1), not reach(80, L+1), not reach(93, L+1).
reach(56,L + 1) :- reach(43, L), L <= 13, L>= 4, not reach(30, L+1), not reach(36, L+1), not reach(80, L+1), not reach(93, L+1).
reach(80,L + 1) :- reach(43, L), L <= 13, L>= 4, not reach(30, L+1), not reach(36, L+1), not reach(56, L+1), not reach(93, L+1).
reach(93,L + 1) :- reach(43, L), L <= 13, L>= 4, not reach(30, L+1), not reach(36, L+1), not reach(56, L+1), not reach(80, L+1).
reach(30,L + 1) :- reach(41, L), L <= 13, L>= 4, not reach(34, L+1), not reach(80, L+1).
reach(34,L + 1) :- reach(41, L), L <= 13, L>= 4, not reach(30, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(41, L), L <= 13, L>= 4, not reach(30, L+1), not reach(34, L+1).
reach(31,L + 1) :- reach(39, L), L <= 13, L>= 2, not reach(71, L+1), not reach(25, L+1), not reach(66, L+1).
reach(71,L + 1) :- reach(39, L), L <= 13, L>= 2, not reach(31, L+1), not reach(25, L+1), not reach(66, L+1).
reach(25,L + 1) :- reach(39, L), L <= 13, L>= 2, not reach(31, L+1), not reach(71, L+1), not reach(66, L+1).
reach(66,L + 1) :- reach(39, L), L <= 13, L>= 2, not reach(31, L+1), not reach(71, L+1), not reach(25, L+1).
reach(31,L + 1) :- reach(49, L), L <= 13, L>= 2, not reach(51, L+1), not reach(80, L+1), not reach(26, L+1).
reach(51,L + 1) :- reach(49, L), L <= 13, L>= 2, not reach(31, L+1), not reach(80, L+1), not reach(26, L+1).
reach(80,L + 1) :- reach(49, L), L <= 13, L>= 2, not reach(31, L+1), not reach(51, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(49, L), L <= 13, L>= 2, not reach(31, L+1), not reach(51, L+1), not reach(80, L+1).
reach(34,L + 1) :- reach(48, L), L <= 15, L>= 3, not reach(38, L+1), not reach(66, L+1), not reach(68, L+1), not reach(80, L+1), not reach(4, L+1).
reach(38,L + 1) :- reach(48, L), L <= 15, L>= 3, not reach(34, L+1), not reach(66, L+1), not reach(68, L+1), not reach(80, L+1), not reach(4, L+1).
reach(66,L + 1) :- reach(48, L), L <= 15, L>= 3, not reach(34, L+1), not reach(38, L+1), not reach(68, L+1), not reach(80, L+1), not reach(4, L+1).
reach(68,L + 1) :- reach(48, L), L <= 15, L>= 3, not reach(34, L+1), not reach(38, L+1), not reach(66, L+1), not reach(80, L+1), not reach(4, L+1).
reach(80,L + 1) :- reach(48, L), L <= 15, L>= 3, not reach(34, L+1), not reach(38, L+1), not reach(66, L+1), not reach(68, L+1), not reach(4, L+1).
reach(4,L + 1) :- reach(48, L), L <= 15, L>= 3, not reach(34, L+1), not reach(38, L+1), not reach(66, L+1), not reach(68, L+1), not reach(80, L+1).
reach(35,L + 1) :- reach(47, L), L <= 15, L>= 4, not reach(46, L+1), not reach(77, L+1), not reach(80, L+1), not reach(38, L+1).
reach(46,L + 1) :- reach(47, L), L <= 15, L>= 4, not reach(35, L+1), not reach(77, L+1), not reach(80, L+1), not reach(38, L+1).
reach(77,L + 1) :- reach(47, L), L <= 15, L>= 4, not reach(35, L+1), not reach(46, L+1), not reach(80, L+1), not reach(38, L+1).
reach(80,L + 1) :- reach(47, L), L <= 15, L>= 4, not reach(35, L+1), not reach(46, L+1), not reach(77, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(47, L), L <= 15, L>= 4, not reach(35, L+1), not reach(46, L+1), not reach(77, L+1), not reach(80, L+1).
reach(43,L + 1) :- reach(36, L), L <= 13, L>= 3, not reach(81, L+1), not reach(77, L+1).
reach(81,L + 1) :- reach(36, L), L <= 13, L>= 3, not reach(43, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(36, L), L <= 13, L>= 3, not reach(43, L+1), not reach(81, L+1).
reach(45,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(57, L+1), not reach(44, L+1).
reach(57,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(45, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(88, L), L <= 16, L>= 4, not reach(45, L+1), not reach(57, L+1).
reach(45,L + 1) :- reach(38, L), L <= 16, L>= 4, not reach(48, L+1), not reach(47, L+1).
reach(48,L + 1) :- reach(38, L), L <= 16, L>= 4, not reach(45, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(38, L), L <= 16, L>= 4, not reach(45, L+1), not reach(48, L+1).
reach(40,L + 1) :- reach(57, L), L <= 15, L>= 4, not reach(46, L+1), not reach(80, L+1), not reach(88, L+1), not reach(82, L+1).
reach(46,L + 1) :- reach(57, L), L <= 15, L>= 4, not reach(40, L+1), not reach(80, L+1), not reach(88, L+1), not reach(82, L+1).
reach(80,L + 1) :- reach(57, L), L <= 15, L>= 4, not reach(40, L+1), not reach(46, L+1), not reach(88, L+1), not reach(82, L+1).
reach(88,L + 1) :- reach(57, L), L <= 15, L>= 4, not reach(40, L+1), not reach(46, L+1), not reach(80, L+1), not reach(82, L+1).
reach(82,L + 1) :- reach(57, L), L <= 15, L>= 4, not reach(40, L+1), not reach(46, L+1), not reach(80, L+1), not reach(88, L+1).
reach(43,L + 1) :- reach(93, L), L <= 13, L>= 5, not reach(53, L+1), not reach(68, L+1).
reach(53,L + 1) :- reach(93, L), L <= 13, L>= 5, not reach(43, L+1), not reach(68, L+1).
reach(68,L + 1) :- reach(93, L), L <= 13, L>= 5, not reach(43, L+1), not reach(53, L+1).
reach(57,L + 1) :- reach(46, L), L <= 14, L>= 4, not reach(73, L+1), not reach(47, L+1).
reach(73,L + 1) :- reach(46, L), L <= 14, L>= 4, not reach(57, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(46, L), L <= 14, L>= 4, not reach(57, L+1), not reach(73, L+1).
reach(48,L + 1) :- reach(68, L), L <= 14, L>= 4, not reach(62, L+1), not reach(93, L+1).
reach(62,L + 1) :- reach(68, L), L <= 14, L>= 4, not reach(48, L+1), not reach(93, L+1).
reach(93,L + 1) :- reach(68, L), L <= 14, L>= 4, not reach(48, L+1), not reach(62, L+1).
reach(51,L + 1) :- reach(87, L), L <= 13, L>= 3, not reach(72, L+1), not reach(91, L+1).
reach(72,L + 1) :- reach(87, L), L <= 13, L>= 3, not reach(51, L+1), not reach(91, L+1).
reach(91,L + 1) :- reach(87, L), L <= 13, L>= 3, not reach(51, L+1), not reach(72, L+1).
reach(52,L + 1) :- reach(89, L), L <= 12, L>= 3, not reach(73, L+1), not reach(81, L+1).
reach(73,L + 1) :- reach(89, L), L <= 12, L>= 3, not reach(52, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(89, L), L <= 12, L>= 3, not reach(52, L+1), not reach(73, L+1).
reach(67,L + 1) :- reach(62, L), L <= 13, L>= 4, not reach(68, L+1), not reach(77, L+1).
reach(68,L + 1) :- reach(62, L), L <= 13, L>= 4, not reach(67, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(62, L), L <= 13, L>= 4, not reach(67, L+1), not reach(68, L+1).
reach(62,L + 1) :- reach(67, L), L <= 13, L>= 3, not reach(77, L+1), not reach(52, L+1).
reach(77,L + 1) :- reach(67, L), L <= 13, L>= 3, not reach(62, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(67, L), L <= 13, L>= 3, not reach(62, L+1), not reach(77, L+1).
reach(72,L + 1) :- reach(69, L), L <= 12, L>= 2, not reach(81, L+1), not reach(13, L+1).
reach(81,L + 1) :- reach(69, L), L <= 12, L>= 2, not reach(72, L+1), not reach(13, L+1).
reach(13,L + 1) :- reach(69, L), L <= 12, L>= 2, not reach(72, L+1), not reach(81, L+1).

