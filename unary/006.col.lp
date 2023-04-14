reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17), not reach(X, 18), not reach(X, 19).
start(43).
goal(241).
:- reach(X, L), reach(X, L'), L != L'.
reach(6,L + 1) :- reach(49, L), L <= 7, L>= 7, not reach(50, L+1), not reach(65, L+1).
reach(50,L + 1) :- reach(49, L), L <= 7, L>= 7, not reach(6, L+1), not reach(65, L+1).
reach(65,L + 1) :- reach(49, L), L <= 7, L>= 7, not reach(6, L+1), not reach(50, L+1).
reach(9,L + 1) :- reach(50, L), L <= 6, L>= 6, not reach(49, L+1), not reach(51, L+1), not reach(66, L+1).
reach(49,L + 1) :- reach(50, L), L <= 6, L>= 6, not reach(9, L+1), not reach(51, L+1), not reach(66, L+1).
reach(51,L + 1) :- reach(50, L), L <= 6, L>= 6, not reach(9, L+1), not reach(49, L+1), not reach(66, L+1).
reach(66,L + 1) :- reach(50, L), L <= 6, L>= 6, not reach(9, L+1), not reach(49, L+1), not reach(51, L+1).
reach(10,L + 1) :- reach(13, L), L <= 11, L>= 10, not reach(16, L+1), not reach(17, L+1), not reach(248, L+1).
reach(16,L + 1) :- reach(13, L), L <= 11, L>= 10, not reach(10, L+1), not reach(17, L+1), not reach(248, L+1).
reach(17,L + 1) :- reach(13, L), L <= 11, L>= 10, not reach(10, L+1), not reach(16, L+1), not reach(248, L+1).
reach(248,L + 1) :- reach(13, L), L <= 11, L>= 10, not reach(10, L+1), not reach(16, L+1), not reach(17, L+1).
reach(13,L + 1) :- reach(16, L), L <= 10, L>= 9, not reach(19, L+1), not reach(20, L+1).
reach(19,L + 1) :- reach(16, L), L <= 10, L>= 9, not reach(13, L+1), not reach(20, L+1).
reach(20,L + 1) :- reach(16, L), L <= 10, L>= 9, not reach(13, L+1), not reach(19, L+1).
reach(13,L + 1) :- reach(248, L), L <= 12, L>= 11, not reach(232, L+1), not reach(247, L+1), not reach(249, L+1).
reach(232,L + 1) :- reach(248, L), L <= 12, L>= 11, not reach(13, L+1), not reach(247, L+1), not reach(249, L+1).
reach(247,L + 1) :- reach(248, L), L <= 12, L>= 11, not reach(13, L+1), not reach(232, L+1), not reach(249, L+1).
reach(249,L + 1) :- reach(248, L), L <= 12, L>= 11, not reach(13, L+1), not reach(232, L+1), not reach(247, L+1).
reach(16,L + 1) :- reach(19, L), L <= 9, L>= 8, not reach(22, L+1), not reach(23, L+1).
reach(22,L + 1) :- reach(19, L), L <= 9, L>= 8, not reach(16, L+1), not reach(23, L+1).
reach(23,L + 1) :- reach(19, L), L <= 9, L>= 8, not reach(16, L+1), not reach(22, L+1).
reach(19,L + 1) :- reach(22, L), L <= 8, L>= 7, not reach(25, L+1), not reach(26, L+1).
reach(25,L + 1) :- reach(22, L), L <= 8, L>= 7, not reach(19, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(22, L), L <= 8, L>= 7, not reach(19, L+1), not reach(25, L+1).
reach(22,L + 1) :- reach(25, L), L <= 7, L>= 6, not reach(28, L+1), not reach(29, L+1).
reach(28,L + 1) :- reach(25, L), L <= 7, L>= 6, not reach(22, L+1), not reach(29, L+1).
reach(29,L + 1) :- reach(25, L), L <= 7, L>= 6, not reach(22, L+1), not reach(28, L+1).
reach(25,L + 1) :- reach(28, L), L <= 6, L>= 5, not reach(31, L+1), not reach(32, L+1).
reach(31,L + 1) :- reach(28, L), L <= 6, L>= 5, not reach(25, L+1), not reach(32, L+1).
reach(32,L + 1) :- reach(28, L), L <= 6, L>= 5, not reach(25, L+1), not reach(31, L+1).
reach(28,L + 1) :- reach(31, L), L <= 5, L>= 4, not reach(34, L+1), not reach(35, L+1).
reach(34,L + 1) :- reach(31, L), L <= 5, L>= 4, not reach(28, L+1), not reach(35, L+1).
reach(35,L + 1) :- reach(31, L), L <= 5, L>= 4, not reach(28, L+1), not reach(34, L+1).
reach(31,L + 1) :- reach(34, L), L <= 4, L>= 3, not reach(37, L+1), not reach(38, L+1).
reach(37,L + 1) :- reach(34, L), L <= 4, L>= 3, not reach(31, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(34, L), L <= 4, L>= 3, not reach(31, L+1), not reach(37, L+1).
reach(32,L + 1) :- reach(36, L), L <= 5, L>= 5, not reach(33, L+1), not reach(39, L+1), not reach(59, L+1).
reach(33,L + 1) :- reach(36, L), L <= 5, L>= 5, not reach(32, L+1), not reach(39, L+1), not reach(59, L+1).
reach(39,L + 1) :- reach(36, L), L <= 5, L>= 5, not reach(32, L+1), not reach(33, L+1), not reach(59, L+1).
reach(59,L + 1) :- reach(36, L), L <= 5, L>= 5, not reach(32, L+1), not reach(33, L+1), not reach(39, L+1).
reach(34,L + 1) :- reach(37, L), L <= 3, L>= 2, not reach(40, L+1), not reach(41, L+1).
reach(40,L + 1) :- reach(37, L), L <= 3, L>= 2, not reach(34, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(37, L), L <= 3, L>= 2, not reach(34, L+1), not reach(40, L+1).
reach(35,L + 1) :- reach(39, L), L <= 4, L>= 4, not reach(36, L+1), not reach(42, L+1), not reach(60, L+1), not reach(41, L+1).
reach(36,L + 1) :- reach(39, L), L <= 4, L>= 4, not reach(35, L+1), not reach(42, L+1), not reach(60, L+1), not reach(41, L+1).
reach(42,L + 1) :- reach(39, L), L <= 4, L>= 4, not reach(35, L+1), not reach(36, L+1), not reach(60, L+1), not reach(41, L+1).
reach(60,L + 1) :- reach(39, L), L <= 4, L>= 4, not reach(35, L+1), not reach(36, L+1), not reach(42, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(39, L), L <= 4, L>= 4, not reach(35, L+1), not reach(36, L+1), not reach(42, L+1), not reach(60, L+1).
reach(36,L + 1) :- reach(59, L), L <= 6, L>= 6, not reach(58, L+1), not reach(60, L+1), not reach(75, L+1), not reach(96, L+1).
reach(58,L + 1) :- reach(59, L), L <= 6, L>= 6, not reach(36, L+1), not reach(60, L+1), not reach(75, L+1), not reach(96, L+1).
reach(60,L + 1) :- reach(59, L), L <= 6, L>= 6, not reach(36, L+1), not reach(58, L+1), not reach(75, L+1), not reach(96, L+1).
reach(75,L + 1) :- reach(59, L), L <= 6, L>= 6, not reach(36, L+1), not reach(58, L+1), not reach(60, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(59, L), L <= 6, L>= 6, not reach(36, L+1), not reach(58, L+1), not reach(60, L+1), not reach(75, L+1).
reach(37,L + 1) :- reach(40, L), L <= 3, L>= 1, not reach(43, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(40, L), L <= 3, L>= 1, not reach(37, L+1), not reach(43, L+1).
reach(37,L + 1) :- reach(41, L), L <= 3, L>= 3, not reach(38, L+1), not reach(39, L+1), not reach(44, L+1), not reach(45, L+1).
reach(38,L + 1) :- reach(41, L), L <= 3, L>= 3, not reach(37, L+1), not reach(39, L+1), not reach(44, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(41, L), L <= 3, L>= 3, not reach(37, L+1), not reach(38, L+1), not reach(44, L+1), not reach(45, L+1).
reach(44,L + 1) :- reach(41, L), L <= 3, L>= 3, not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(41, L), L <= 3, L>= 3, not reach(37, L+1), not reach(38, L+1), not reach(39, L+1), not reach(44, L+1).
reach(39,L + 1) :- reach(60, L), L <= 5, L>= 5, not reach(59, L+1), not reach(61, L+1), not reach(76, L+1).
reach(59,L + 1) :- reach(60, L), L <= 5, L>= 5, not reach(39, L+1), not reach(61, L+1), not reach(76, L+1).
reach(61,L + 1) :- reach(60, L), L <= 5, L>= 5, not reach(39, L+1), not reach(59, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(60, L), L <= 5, L>= 5, not reach(39, L+1), not reach(59, L+1), not reach(61, L+1).
reach(40,L + 1) :- reach(43, L), L <= 2, L>= 0, not reach(46, L+1).
reach(46,L + 1) :- reach(43, L), L <= 2, L>= 0, not reach(40, L+1).
reach(40,L + 1) :- reach(44, L), L <= 4, L>= 2, not reach(41, L+1), not reach(47, L+1).
reach(41,L + 1) :- reach(44, L), L <= 4, L>= 2, not reach(40, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(44, L), L <= 4, L>= 2, not reach(40, L+1), not reach(41, L+1).
reach(41,L + 1) :- reach(45, L), L <= 4, L>= 4, not reach(42, L+1), not reach(47, L+1), not reach(62, L+1).
reach(42,L + 1) :- reach(45, L), L <= 4, L>= 4, not reach(41, L+1), not reach(47, L+1), not reach(62, L+1).
reach(47,L + 1) :- reach(45, L), L <= 4, L>= 4, not reach(41, L+1), not reach(42, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(45, L), L <= 4, L>= 4, not reach(41, L+1), not reach(42, L+1), not reach(47, L+1).
reach(48,L + 1) :- reach(46, L), L <= 3, L>= 1, not reach(43, L+1).
reach(44,L + 1) :- reach(47, L), L <= 5, L>= 3, not reach(45, L+1), not reach(48, L+1), not reach(63, L+1).
reach(45,L + 1) :- reach(47, L), L <= 5, L>= 3, not reach(44, L+1), not reach(48, L+1), not reach(63, L+1).
reach(48,L + 1) :- reach(47, L), L <= 5, L>= 3, not reach(44, L+1), not reach(45, L+1), not reach(63, L+1).
reach(63,L + 1) :- reach(47, L), L <= 5, L>= 3, not reach(44, L+1), not reach(45, L+1), not reach(48, L+1).
reach(45,L + 1) :- reach(62, L), L <= 5, L>= 5, not reach(61, L+1), not reach(78, L+1).
reach(61,L + 1) :- reach(62, L), L <= 5, L>= 5, not reach(45, L+1), not reach(78, L+1).
reach(78,L + 1) :- reach(62, L), L <= 5, L>= 5, not reach(45, L+1), not reach(61, L+1).
reach(46,L + 1) :- reach(48, L), L <= 4, L>= 2, not reach(47, L+1), not reach(64, L+1).
reach(47,L + 1) :- reach(48, L), L <= 4, L>= 2, not reach(46, L+1), not reach(64, L+1).
reach(64,L + 1) :- reach(48, L), L <= 4, L>= 2, not reach(46, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(63, L), L <= 6, L>= 4, not reach(64, L+1), not reach(79, L+1), not reach(66, L+1).
reach(64,L + 1) :- reach(63, L), L <= 6, L>= 4, not reach(47, L+1), not reach(79, L+1), not reach(66, L+1).
reach(79,L + 1) :- reach(63, L), L <= 6, L>= 4, not reach(47, L+1), not reach(64, L+1), not reach(66, L+1).
reach(66,L + 1) :- reach(63, L), L <= 6, L>= 4, not reach(47, L+1), not reach(64, L+1), not reach(79, L+1).
reach(48,L + 1) :- reach(64, L), L <= 5, L>= 3, not reach(63, L+1), not reach(80, L+1).
reach(63,L + 1) :- reach(64, L), L <= 5, L>= 3, not reach(48, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(64, L), L <= 5, L>= 3, not reach(48, L+1), not reach(63, L+1).
reach(49,L + 1) :- reach(65, L), L <= 8, L>= 6, not reach(66, L+1), not reach(81, L+1).
reach(66,L + 1) :- reach(65, L), L <= 8, L>= 6, not reach(49, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(65, L), L <= 8, L>= 6, not reach(49, L+1), not reach(66, L+1).
reach(50,L + 1) :- reach(66, L), L <= 7, L>= 5, not reach(63, L+1), not reach(65, L+1), not reach(67, L+1), not reach(82, L+1).
reach(63,L + 1) :- reach(66, L), L <= 7, L>= 5, not reach(50, L+1), not reach(65, L+1), not reach(67, L+1), not reach(82, L+1).
reach(65,L + 1) :- reach(66, L), L <= 7, L>= 5, not reach(50, L+1), not reach(63, L+1), not reach(67, L+1), not reach(82, L+1).
reach(67,L + 1) :- reach(66, L), L <= 7, L>= 5, not reach(50, L+1), not reach(63, L+1), not reach(65, L+1), not reach(82, L+1).
reach(82,L + 1) :- reach(66, L), L <= 7, L>= 5, not reach(50, L+1), not reach(63, L+1), not reach(65, L+1), not reach(67, L+1).
reach(51,L + 1) :- reach(67, L), L <= 6, L>= 6, not reach(66, L+1), not reach(68, L+1), not reach(83, L+1).
reach(66,L + 1) :- reach(67, L), L <= 6, L>= 6, not reach(51, L+1), not reach(68, L+1), not reach(83, L+1).
reach(68,L + 1) :- reach(67, L), L <= 6, L>= 6, not reach(51, L+1), not reach(66, L+1), not reach(83, L+1).
reach(83,L + 1) :- reach(67, L), L <= 6, L>= 6, not reach(51, L+1), not reach(66, L+1), not reach(68, L+1).
reach(59,L + 1) :- reach(96, L), L <= 7, L>= 5, not reach(80, L+1), not reach(95, L+1), not reach(112, L+1).
reach(80,L + 1) :- reach(96, L), L <= 7, L>= 5, not reach(59, L+1), not reach(95, L+1), not reach(112, L+1).
reach(95,L + 1) :- reach(96, L), L <= 7, L>= 5, not reach(59, L+1), not reach(80, L+1), not reach(112, L+1).
reach(112,L + 1) :- reach(96, L), L <= 7, L>= 5, not reach(59, L+1), not reach(80, L+1), not reach(95, L+1).
reach(62,L + 1) :- reach(78, L), L <= 6, L>= 6, not reach(77, L+1), not reach(79, L+1), not reach(94, L+1).
reach(77,L + 1) :- reach(78, L), L <= 6, L>= 6, not reach(62, L+1), not reach(79, L+1), not reach(94, L+1).
reach(79,L + 1) :- reach(78, L), L <= 6, L>= 6, not reach(62, L+1), not reach(77, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(78, L), L <= 6, L>= 6, not reach(62, L+1), not reach(77, L+1), not reach(79, L+1).
reach(63,L + 1) :- reach(79, L), L <= 7, L>= 5, not reach(78, L+1), not reach(80, L+1), not reach(95, L+1).
reach(78,L + 1) :- reach(79, L), L <= 7, L>= 5, not reach(63, L+1), not reach(80, L+1), not reach(95, L+1).
reach(80,L + 1) :- reach(79, L), L <= 7, L>= 5, not reach(63, L+1), not reach(78, L+1), not reach(95, L+1).
reach(95,L + 1) :- reach(79, L), L <= 7, L>= 5, not reach(63, L+1), not reach(78, L+1), not reach(80, L+1).
reach(64,L + 1) :- reach(80, L), L <= 6, L>= 4, not reach(79, L+1), not reach(96, L+1).
reach(79,L + 1) :- reach(80, L), L <= 6, L>= 4, not reach(64, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(80, L), L <= 6, L>= 4, not reach(64, L+1), not reach(79, L+1).
reach(65,L + 1) :- reach(81, L), L <= 9, L>= 7, not reach(82, L+1), not reach(97, L+1).
reach(82,L + 1) :- reach(81, L), L <= 9, L>= 7, not reach(65, L+1), not reach(97, L+1).
reach(97,L + 1) :- reach(81, L), L <= 9, L>= 7, not reach(65, L+1), not reach(82, L+1).
reach(66,L + 1) :- reach(82, L), L <= 8, L>= 6, not reach(81, L+1), not reach(83, L+1), not reach(98, L+1).
reach(81,L + 1) :- reach(82, L), L <= 8, L>= 6, not reach(66, L+1), not reach(83, L+1), not reach(98, L+1).
reach(83,L + 1) :- reach(82, L), L <= 8, L>= 6, not reach(66, L+1), not reach(81, L+1), not reach(98, L+1).
reach(98,L + 1) :- reach(82, L), L <= 8, L>= 6, not reach(66, L+1), not reach(81, L+1), not reach(83, L+1).
reach(67,L + 1) :- reach(83, L), L <= 7, L>= 7, not reach(82, L+1), not reach(84, L+1), not reach(99, L+1).
reach(82,L + 1) :- reach(83, L), L <= 7, L>= 7, not reach(67, L+1), not reach(84, L+1), not reach(99, L+1).
reach(84,L + 1) :- reach(83, L), L <= 7, L>= 7, not reach(67, L+1), not reach(82, L+1), not reach(99, L+1).
reach(99,L + 1) :- reach(83, L), L <= 7, L>= 7, not reach(67, L+1), not reach(82, L+1), not reach(84, L+1).
reach(78,L + 1) :- reach(94, L), L <= 7, L>= 7, not reach(93, L+1), not reach(95, L+1), not reach(110, L+1).
reach(93,L + 1) :- reach(94, L), L <= 7, L>= 7, not reach(78, L+1), not reach(95, L+1), not reach(110, L+1).
reach(95,L + 1) :- reach(94, L), L <= 7, L>= 7, not reach(78, L+1), not reach(93, L+1), not reach(110, L+1).
reach(110,L + 1) :- reach(94, L), L <= 7, L>= 7, not reach(78, L+1), not reach(93, L+1), not reach(95, L+1).
reach(79,L + 1) :- reach(95, L), L <= 8, L>= 6, not reach(94, L+1), not reach(96, L+1), not reach(111, L+1).
reach(94,L + 1) :- reach(95, L), L <= 8, L>= 6, not reach(79, L+1), not reach(96, L+1), not reach(111, L+1).
reach(96,L + 1) :- reach(95, L), L <= 8, L>= 6, not reach(79, L+1), not reach(94, L+1), not reach(111, L+1).
reach(111,L + 1) :- reach(95, L), L <= 8, L>= 6, not reach(79, L+1), not reach(94, L+1), not reach(96, L+1).
reach(81,L + 1) :- reach(97, L), L <= 10, L>= 8, not reach(98, L+1), not reach(113, L+1).
reach(98,L + 1) :- reach(97, L), L <= 10, L>= 8, not reach(81, L+1), not reach(113, L+1).
reach(113,L + 1) :- reach(97, L), L <= 10, L>= 8, not reach(81, L+1), not reach(98, L+1).
reach(82,L + 1) :- reach(98, L), L <= 9, L>= 7, not reach(97, L+1), not reach(99, L+1), not reach(114, L+1).
reach(97,L + 1) :- reach(98, L), L <= 9, L>= 7, not reach(82, L+1), not reach(99, L+1), not reach(114, L+1).
reach(99,L + 1) :- reach(98, L), L <= 9, L>= 7, not reach(82, L+1), not reach(97, L+1), not reach(114, L+1).
reach(114,L + 1) :- reach(98, L), L <= 9, L>= 7, not reach(82, L+1), not reach(97, L+1), not reach(99, L+1).
reach(83,L + 1) :- reach(99, L), L <= 8, L>= 8, not reach(98, L+1), not reach(100, L+1), not reach(115, L+1).
reach(98,L + 1) :- reach(99, L), L <= 8, L>= 8, not reach(83, L+1), not reach(100, L+1), not reach(115, L+1).
reach(100,L + 1) :- reach(99, L), L <= 8, L>= 8, not reach(83, L+1), not reach(98, L+1), not reach(115, L+1).
reach(115,L + 1) :- reach(99, L), L <= 8, L>= 8, not reach(83, L+1), not reach(98, L+1), not reach(100, L+1).
reach(88,L + 1) :- reach(147, L), L <= 11, L>= 11, not reach(131, L+1), not reach(146, L+1), not reach(148, L+1), not reach(162, L+1).
reach(131,L + 1) :- reach(147, L), L <= 11, L>= 11, not reach(88, L+1), not reach(146, L+1), not reach(148, L+1), not reach(162, L+1).
reach(146,L + 1) :- reach(147, L), L <= 11, L>= 11, not reach(88, L+1), not reach(131, L+1), not reach(148, L+1), not reach(162, L+1).
reach(148,L + 1) :- reach(147, L), L <= 11, L>= 11, not reach(88, L+1), not reach(131, L+1), not reach(146, L+1), not reach(162, L+1).
reach(162,L + 1) :- reach(147, L), L <= 11, L>= 11, not reach(88, L+1), not reach(131, L+1), not reach(146, L+1), not reach(148, L+1).
reach(94,L + 1) :- reach(110, L), L <= 8, L>= 8, not reach(109, L+1), not reach(111, L+1), not reach(126, L+1).
reach(109,L + 1) :- reach(110, L), L <= 8, L>= 8, not reach(94, L+1), not reach(111, L+1), not reach(126, L+1).
reach(111,L + 1) :- reach(110, L), L <= 8, L>= 8, not reach(94, L+1), not reach(109, L+1), not reach(126, L+1).
reach(126,L + 1) :- reach(110, L), L <= 8, L>= 8, not reach(94, L+1), not reach(109, L+1), not reach(111, L+1).
reach(95,L + 1) :- reach(111, L), L <= 9, L>= 7, not reach(110, L+1), not reach(112, L+1), not reach(127, L+1), not reach(199, L+1).
reach(110,L + 1) :- reach(111, L), L <= 9, L>= 7, not reach(95, L+1), not reach(112, L+1), not reach(127, L+1), not reach(199, L+1).
reach(112,L + 1) :- reach(111, L), L <= 9, L>= 7, not reach(95, L+1), not reach(110, L+1), not reach(127, L+1), not reach(199, L+1).
reach(127,L + 1) :- reach(111, L), L <= 9, L>= 7, not reach(95, L+1), not reach(110, L+1), not reach(112, L+1), not reach(199, L+1).
reach(199,L + 1) :- reach(111, L), L <= 9, L>= 7, not reach(95, L+1), not reach(110, L+1), not reach(112, L+1), not reach(127, L+1).
reach(96,L + 1) :- reach(112, L), L <= 8, L>= 6, not reach(111, L+1), not reach(128, L+1).
reach(111,L + 1) :- reach(112, L), L <= 8, L>= 6, not reach(96, L+1), not reach(128, L+1).
reach(128,L + 1) :- reach(112, L), L <= 8, L>= 6, not reach(96, L+1), not reach(111, L+1).
reach(97,L + 1) :- reach(113, L), L <= 11, L>= 9, not reach(114, L+1), not reach(129, L+1).
reach(114,L + 1) :- reach(113, L), L <= 11, L>= 9, not reach(97, L+1), not reach(129, L+1).
reach(129,L + 1) :- reach(113, L), L <= 11, L>= 9, not reach(97, L+1), not reach(114, L+1).
reach(98,L + 1) :- reach(114, L), L <= 10, L>= 8, not reach(113, L+1), not reach(115, L+1), not reach(130, L+1).
reach(113,L + 1) :- reach(114, L), L <= 10, L>= 8, not reach(98, L+1), not reach(115, L+1), not reach(130, L+1).
reach(115,L + 1) :- reach(114, L), L <= 10, L>= 8, not reach(98, L+1), not reach(113, L+1), not reach(130, L+1).
reach(130,L + 1) :- reach(114, L), L <= 10, L>= 8, not reach(98, L+1), not reach(113, L+1), not reach(115, L+1).
reach(99,L + 1) :- reach(115, L), L <= 9, L>= 9, not reach(114, L+1), not reach(116, L+1), not reach(131, L+1).
reach(114,L + 1) :- reach(115, L), L <= 9, L>= 9, not reach(99, L+1), not reach(116, L+1), not reach(131, L+1).
reach(116,L + 1) :- reach(115, L), L <= 9, L>= 9, not reach(99, L+1), not reach(114, L+1), not reach(131, L+1).
reach(131,L + 1) :- reach(115, L), L <= 9, L>= 9, not reach(99, L+1), not reach(114, L+1), not reach(116, L+1).
reach(111,L + 1) :- reach(127, L), L <= 8, L>= 8, not reach(126, L+1), not reach(128, L+1), not reach(143, L+1).
reach(126,L + 1) :- reach(127, L), L <= 8, L>= 8, not reach(111, L+1), not reach(128, L+1), not reach(143, L+1).
reach(128,L + 1) :- reach(127, L), L <= 8, L>= 8, not reach(111, L+1), not reach(126, L+1), not reach(143, L+1).
reach(143,L + 1) :- reach(127, L), L <= 8, L>= 8, not reach(111, L+1), not reach(126, L+1), not reach(128, L+1).
reach(111,L + 1) :- reach(199, L), L <= 10, L>= 8, not reach(183, L+1), not reach(198, L+1), not reach(200, L+1), not reach(215, L+1).
reach(183,L + 1) :- reach(199, L), L <= 10, L>= 8, not reach(111, L+1), not reach(198, L+1), not reach(200, L+1), not reach(215, L+1).
reach(198,L + 1) :- reach(199, L), L <= 10, L>= 8, not reach(111, L+1), not reach(183, L+1), not reach(200, L+1), not reach(215, L+1).
reach(200,L + 1) :- reach(199, L), L <= 10, L>= 8, not reach(111, L+1), not reach(183, L+1), not reach(198, L+1), not reach(215, L+1).
reach(215,L + 1) :- reach(199, L), L <= 10, L>= 8, not reach(111, L+1), not reach(183, L+1), not reach(198, L+1), not reach(200, L+1).
reach(112,L + 1) :- reach(128, L), L <= 7, L>= 7, not reach(127, L+1), not reach(144, L+1).
reach(127,L + 1) :- reach(128, L), L <= 7, L>= 7, not reach(112, L+1), not reach(144, L+1).
reach(144,L + 1) :- reach(128, L), L <= 7, L>= 7, not reach(112, L+1), not reach(127, L+1).
reach(113,L + 1) :- reach(129, L), L <= 12, L>= 10, not reach(130, L+1), not reach(145, L+1).
reach(130,L + 1) :- reach(129, L), L <= 12, L>= 10, not reach(113, L+1), not reach(145, L+1).
reach(145,L + 1) :- reach(129, L), L <= 12, L>= 10, not reach(113, L+1), not reach(130, L+1).
reach(114,L + 1) :- reach(130, L), L <= 11, L>= 9, not reach(129, L+1), not reach(131, L+1), not reach(146, L+1).
reach(129,L + 1) :- reach(130, L), L <= 11, L>= 9, not reach(114, L+1), not reach(131, L+1), not reach(146, L+1).
reach(131,L + 1) :- reach(130, L), L <= 11, L>= 9, not reach(114, L+1), not reach(129, L+1), not reach(146, L+1).
reach(146,L + 1) :- reach(130, L), L <= 11, L>= 9, not reach(114, L+1), not reach(129, L+1), not reach(131, L+1).
reach(115,L + 1) :- reach(131, L), L <= 10, L>= 10, not reach(130, L+1), not reach(132, L+1), not reach(147, L+1).
reach(130,L + 1) :- reach(131, L), L <= 10, L>= 10, not reach(115, L+1), not reach(132, L+1), not reach(147, L+1).
reach(132,L + 1) :- reach(131, L), L <= 10, L>= 10, not reach(115, L+1), not reach(130, L+1), not reach(147, L+1).
reach(147,L + 1) :- reach(131, L), L <= 10, L>= 10, not reach(115, L+1), not reach(130, L+1), not reach(132, L+1).
reach(129,L + 1) :- reach(145, L), L <= 13, L>= 11, not reach(146, L+1), not reach(160, L+1).
reach(146,L + 1) :- reach(145, L), L <= 13, L>= 11, not reach(129, L+1), not reach(160, L+1).
reach(160,L + 1) :- reach(145, L), L <= 13, L>= 11, not reach(129, L+1), not reach(146, L+1).
reach(130,L + 1) :- reach(146, L), L <= 12, L>= 10, not reach(145, L+1), not reach(147, L+1), not reach(161, L+1).
reach(145,L + 1) :- reach(146, L), L <= 12, L>= 10, not reach(130, L+1), not reach(147, L+1), not reach(161, L+1).
reach(147,L + 1) :- reach(146, L), L <= 12, L>= 10, not reach(130, L+1), not reach(145, L+1), not reach(161, L+1).
reach(161,L + 1) :- reach(146, L), L <= 12, L>= 10, not reach(130, L+1), not reach(145, L+1), not reach(147, L+1).
reach(145,L + 1) :- reach(160, L), L <= 14, L>= 12, not reach(161, L+1), not reach(177, L+1).
reach(161,L + 1) :- reach(160, L), L <= 14, L>= 12, not reach(145, L+1), not reach(177, L+1).
reach(177,L + 1) :- reach(160, L), L <= 14, L>= 12, not reach(145, L+1), not reach(161, L+1).
reach(146,L + 1) :- reach(161, L), L <= 13, L>= 11, not reach(160, L+1), not reach(162, L+1), not reach(178, L+1).
reach(160,L + 1) :- reach(161, L), L <= 13, L>= 11, not reach(146, L+1), not reach(162, L+1), not reach(178, L+1).
reach(162,L + 1) :- reach(161, L), L <= 13, L>= 11, not reach(146, L+1), not reach(160, L+1), not reach(178, L+1).
reach(178,L + 1) :- reach(161, L), L <= 13, L>= 11, not reach(146, L+1), not reach(160, L+1), not reach(162, L+1).
reach(147,L + 1) :- reach(162, L), L <= 12, L>= 12, not reach(161, L+1), not reach(163, L+1), not reach(179, L+1).
reach(161,L + 1) :- reach(162, L), L <= 12, L>= 12, not reach(147, L+1), not reach(163, L+1), not reach(179, L+1).
reach(163,L + 1) :- reach(162, L), L <= 12, L>= 12, not reach(147, L+1), not reach(161, L+1), not reach(179, L+1).
reach(179,L + 1) :- reach(162, L), L <= 12, L>= 12, not reach(147, L+1), not reach(161, L+1), not reach(163, L+1).
reach(156,L + 1) :- reach(183, L), L <= 9, L>= 9, not reach(166, L+1), not reach(182, L+1), not reach(184, L+1), not reach(199, L+1).
reach(166,L + 1) :- reach(183, L), L <= 9, L>= 9, not reach(156, L+1), not reach(182, L+1), not reach(184, L+1), not reach(199, L+1).
reach(182,L + 1) :- reach(183, L), L <= 9, L>= 9, not reach(156, L+1), not reach(166, L+1), not reach(184, L+1), not reach(199, L+1).
reach(184,L + 1) :- reach(183, L), L <= 9, L>= 9, not reach(156, L+1), not reach(166, L+1), not reach(182, L+1), not reach(199, L+1).
reach(199,L + 1) :- reach(183, L), L <= 9, L>= 9, not reach(156, L+1), not reach(166, L+1), not reach(182, L+1), not reach(184, L+1).
reach(160,L + 1) :- reach(177, L), L <= 15, L>= 13, not reach(193, L+1).
reach(193,L + 1) :- reach(177, L), L <= 15, L>= 13, not reach(160, L+1).
reach(161,L + 1) :- reach(178, L), L <= 14, L>= 12, not reach(179, L+1), not reach(194, L+1).
reach(179,L + 1) :- reach(178, L), L <= 14, L>= 12, not reach(161, L+1), not reach(194, L+1).
reach(194,L + 1) :- reach(178, L), L <= 14, L>= 12, not reach(161, L+1), not reach(179, L+1).
reach(162,L + 1) :- reach(179, L), L <= 13, L>= 13, not reach(178, L+1), not reach(180, L+1).
reach(178,L + 1) :- reach(179, L), L <= 13, L>= 13, not reach(162, L+1), not reach(180, L+1).
reach(180,L + 1) :- reach(179, L), L <= 13, L>= 13, not reach(162, L+1), not reach(178, L+1).
reach(163,L + 1) :- reach(180, L), L <= 12, L>= 12, not reach(179, L+1), not reach(181, L+1), not reach(196, L+1).
reach(179,L + 1) :- reach(180, L), L <= 12, L>= 12, not reach(163, L+1), not reach(181, L+1), not reach(196, L+1).
reach(181,L + 1) :- reach(180, L), L <= 12, L>= 12, not reach(163, L+1), not reach(179, L+1), not reach(196, L+1).
reach(196,L + 1) :- reach(180, L), L <= 12, L>= 12, not reach(163, L+1), not reach(179, L+1), not reach(181, L+1).
reach(164,L + 1) :- reach(181, L), L <= 11, L>= 11, not reach(180, L+1), not reach(182, L+1), not reach(197, L+1).
reach(180,L + 1) :- reach(181, L), L <= 11, L>= 11, not reach(164, L+1), not reach(182, L+1), not reach(197, L+1).
reach(182,L + 1) :- reach(181, L), L <= 11, L>= 11, not reach(164, L+1), not reach(180, L+1), not reach(197, L+1).
reach(197,L + 1) :- reach(181, L), L <= 11, L>= 11, not reach(164, L+1), not reach(180, L+1), not reach(182, L+1).
reach(165,L + 1) :- reach(182, L), L <= 10, L>= 10, not reach(181, L+1), not reach(183, L+1), not reach(198, L+1).
reach(181,L + 1) :- reach(182, L), L <= 10, L>= 10, not reach(165, L+1), not reach(183, L+1), not reach(198, L+1).
reach(183,L + 1) :- reach(182, L), L <= 10, L>= 10, not reach(165, L+1), not reach(181, L+1), not reach(198, L+1).
reach(198,L + 1) :- reach(182, L), L <= 10, L>= 10, not reach(165, L+1), not reach(181, L+1), not reach(183, L+1).
reach(177,L + 1) :- reach(193, L), L <= 16, L>= 14, not reach(194, L+1), not reach(209, L+1).
reach(194,L + 1) :- reach(193, L), L <= 16, L>= 14, not reach(177, L+1), not reach(209, L+1).
reach(209,L + 1) :- reach(193, L), L <= 16, L>= 14, not reach(177, L+1), not reach(194, L+1).
reach(178,L + 1) :- reach(194, L), L <= 15, L>= 13, not reach(193, L+1), not reach(195, L+1), not reach(210, L+1).
reach(193,L + 1) :- reach(194, L), L <= 15, L>= 13, not reach(178, L+1), not reach(195, L+1), not reach(210, L+1).
reach(195,L + 1) :- reach(194, L), L <= 15, L>= 13, not reach(178, L+1), not reach(193, L+1), not reach(210, L+1).
reach(210,L + 1) :- reach(194, L), L <= 15, L>= 13, not reach(178, L+1), not reach(193, L+1), not reach(195, L+1).
reach(180,L + 1) :- reach(196, L), L <= 13, L>= 11, not reach(195, L+1), not reach(197, L+1), not reach(212, L+1).
reach(195,L + 1) :- reach(196, L), L <= 13, L>= 11, not reach(180, L+1), not reach(197, L+1), not reach(212, L+1).
reach(197,L + 1) :- reach(196, L), L <= 13, L>= 11, not reach(180, L+1), not reach(195, L+1), not reach(212, L+1).
reach(212,L + 1) :- reach(196, L), L <= 13, L>= 11, not reach(180, L+1), not reach(195, L+1), not reach(197, L+1).
reach(181,L + 1) :- reach(197, L), L <= 12, L>= 10, not reach(196, L+1), not reach(198, L+1), not reach(213, L+1).
reach(196,L + 1) :- reach(197, L), L <= 12, L>= 10, not reach(181, L+1), not reach(198, L+1), not reach(213, L+1).
reach(198,L + 1) :- reach(197, L), L <= 12, L>= 10, not reach(181, L+1), not reach(196, L+1), not reach(213, L+1).
reach(213,L + 1) :- reach(197, L), L <= 12, L>= 10, not reach(181, L+1), not reach(196, L+1), not reach(198, L+1).
reach(182,L + 1) :- reach(198, L), L <= 11, L>= 9, not reach(197, L+1), not reach(199, L+1), not reach(214, L+1).
reach(197,L + 1) :- reach(198, L), L <= 11, L>= 9, not reach(182, L+1), not reach(199, L+1), not reach(214, L+1).
reach(199,L + 1) :- reach(198, L), L <= 11, L>= 9, not reach(182, L+1), not reach(197, L+1), not reach(214, L+1).
reach(214,L + 1) :- reach(198, L), L <= 11, L>= 9, not reach(182, L+1), not reach(197, L+1), not reach(199, L+1).
reach(184,L + 1) :- reach(200, L), L <= 9, L>= 9, not reach(199, L+1), not reach(201, L+1), not reach(216, L+1).
reach(199,L + 1) :- reach(200, L), L <= 9, L>= 9, not reach(184, L+1), not reach(201, L+1), not reach(216, L+1).
reach(201,L + 1) :- reach(200, L), L <= 9, L>= 9, not reach(184, L+1), not reach(199, L+1), not reach(216, L+1).
reach(216,L + 1) :- reach(200, L), L <= 9, L>= 9, not reach(184, L+1), not reach(199, L+1), not reach(201, L+1).
reach(193,L + 1) :- reach(209, L), L <= 17, L>= 15, not reach(210, L+1), not reach(225, L+1).
reach(210,L + 1) :- reach(209, L), L <= 17, L>= 15, not reach(193, L+1), not reach(225, L+1).
reach(225,L + 1) :- reach(209, L), L <= 17, L>= 15, not reach(193, L+1), not reach(210, L+1).
reach(194,L + 1) :- reach(195, L), L <= 14, L>= 12, not reach(196, L+1), not reach(211, L+1).
reach(196,L + 1) :- reach(195, L), L <= 14, L>= 12, not reach(194, L+1), not reach(211, L+1).
reach(211,L + 1) :- reach(195, L), L <= 14, L>= 12, not reach(194, L+1), not reach(196, L+1).
reach(194,L + 1) :- reach(210, L), L <= 16, L>= 14, not reach(209, L+1), not reach(211, L+1), not reach(226, L+1).
reach(209,L + 1) :- reach(210, L), L <= 16, L>= 14, not reach(194, L+1), not reach(211, L+1), not reach(226, L+1).
reach(211,L + 1) :- reach(210, L), L <= 16, L>= 14, not reach(194, L+1), not reach(209, L+1), not reach(226, L+1).
reach(226,L + 1) :- reach(210, L), L <= 16, L>= 14, not reach(194, L+1), not reach(209, L+1), not reach(211, L+1).
reach(195,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(210, L+1), not reach(212, L+1), not reach(227, L+1).
reach(210,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(195, L+1), not reach(212, L+1), not reach(227, L+1).
reach(212,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(195, L+1), not reach(210, L+1), not reach(227, L+1).
reach(227,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(195, L+1), not reach(210, L+1), not reach(212, L+1).
reach(196,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(211, L+1), not reach(213, L+1), not reach(228, L+1).
reach(211,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(196, L+1), not reach(213, L+1), not reach(228, L+1).
reach(213,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(196, L+1), not reach(211, L+1), not reach(228, L+1).
reach(228,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(196, L+1), not reach(211, L+1), not reach(213, L+1).
reach(197,L + 1) :- reach(213, L), L <= 13, L>= 11, not reach(212, L+1), not reach(214, L+1), not reach(229, L+1).
reach(212,L + 1) :- reach(213, L), L <= 13, L>= 11, not reach(197, L+1), not reach(214, L+1), not reach(229, L+1).
reach(214,L + 1) :- reach(213, L), L <= 13, L>= 11, not reach(197, L+1), not reach(212, L+1), not reach(229, L+1).
reach(229,L + 1) :- reach(213, L), L <= 13, L>= 11, not reach(197, L+1), not reach(212, L+1), not reach(214, L+1).
reach(198,L + 1) :- reach(214, L), L <= 12, L>= 10, not reach(213, L+1), not reach(215, L+1), not reach(230, L+1).
reach(213,L + 1) :- reach(214, L), L <= 12, L>= 10, not reach(198, L+1), not reach(215, L+1), not reach(230, L+1).
reach(215,L + 1) :- reach(214, L), L <= 12, L>= 10, not reach(198, L+1), not reach(213, L+1), not reach(230, L+1).
reach(230,L + 1) :- reach(214, L), L <= 12, L>= 10, not reach(198, L+1), not reach(213, L+1), not reach(215, L+1).
reach(199,L + 1) :- reach(215, L), L <= 11, L>= 9, not reach(214, L+1), not reach(216, L+1), not reach(231, L+1).
reach(214,L + 1) :- reach(215, L), L <= 11, L>= 9, not reach(199, L+1), not reach(216, L+1), not reach(231, L+1).
reach(216,L + 1) :- reach(215, L), L <= 11, L>= 9, not reach(199, L+1), not reach(214, L+1), not reach(231, L+1).
reach(231,L + 1) :- reach(215, L), L <= 11, L>= 9, not reach(199, L+1), not reach(214, L+1), not reach(216, L+1).
reach(200,L + 1) :- reach(216, L), L <= 10, L>= 10, not reach(215, L+1), not reach(217, L+1), not reach(232, L+1).
reach(215,L + 1) :- reach(216, L), L <= 10, L>= 10, not reach(200, L+1), not reach(217, L+1), not reach(232, L+1).
reach(217,L + 1) :- reach(216, L), L <= 10, L>= 10, not reach(200, L+1), not reach(215, L+1), not reach(232, L+1).
reach(232,L + 1) :- reach(216, L), L <= 10, L>= 10, not reach(200, L+1), not reach(215, L+1), not reach(217, L+1).
reach(209,L + 1) :- reach(225, L), L <= 18, L>= 16, not reach(226, L+1), not reach(241, L+1).
reach(226,L + 1) :- reach(225, L), L <= 18, L>= 16, not reach(209, L+1), not reach(241, L+1).
reach(241,L + 1) :- reach(225, L), L <= 18, L>= 16, not reach(209, L+1), not reach(226, L+1).
reach(210,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(225, L+1), not reach(227, L+1), not reach(242, L+1).
reach(225,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(210, L+1), not reach(227, L+1), not reach(242, L+1).
reach(227,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(210, L+1), not reach(225, L+1), not reach(242, L+1).
reach(242,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(210, L+1), not reach(225, L+1), not reach(227, L+1).
reach(211,L + 1) :- reach(227, L), L <= 16, L>= 14, not reach(226, L+1), not reach(228, L+1), not reach(243, L+1).
reach(226,L + 1) :- reach(227, L), L <= 16, L>= 14, not reach(211, L+1), not reach(228, L+1), not reach(243, L+1).
reach(228,L + 1) :- reach(227, L), L <= 16, L>= 14, not reach(211, L+1), not reach(226, L+1), not reach(243, L+1).
reach(243,L + 1) :- reach(227, L), L <= 16, L>= 14, not reach(211, L+1), not reach(226, L+1), not reach(228, L+1).
reach(212,L + 1) :- reach(228, L), L <= 15, L>= 13, not reach(227, L+1), not reach(244, L+1).
reach(227,L + 1) :- reach(228, L), L <= 15, L>= 13, not reach(212, L+1), not reach(244, L+1).
reach(244,L + 1) :- reach(228, L), L <= 15, L>= 13, not reach(212, L+1), not reach(227, L+1).
reach(213,L + 1) :- reach(229, L), L <= 14, L>= 12, not reach(230, L+1), not reach(245, L+1).
reach(230,L + 1) :- reach(229, L), L <= 14, L>= 12, not reach(213, L+1), not reach(245, L+1).
reach(245,L + 1) :- reach(229, L), L <= 14, L>= 12, not reach(213, L+1), not reach(230, L+1).
reach(214,L + 1) :- reach(230, L), L <= 13, L>= 11, not reach(229, L+1), not reach(231, L+1), not reach(246, L+1).
reach(229,L + 1) :- reach(230, L), L <= 13, L>= 11, not reach(214, L+1), not reach(231, L+1), not reach(246, L+1).
reach(231,L + 1) :- reach(230, L), L <= 13, L>= 11, not reach(214, L+1), not reach(229, L+1), not reach(246, L+1).
reach(246,L + 1) :- reach(230, L), L <= 13, L>= 11, not reach(214, L+1), not reach(229, L+1), not reach(231, L+1).
reach(215,L + 1) :- reach(231, L), L <= 12, L>= 10, not reach(230, L+1), not reach(232, L+1), not reach(247, L+1).
reach(230,L + 1) :- reach(231, L), L <= 12, L>= 10, not reach(215, L+1), not reach(232, L+1), not reach(247, L+1).
reach(232,L + 1) :- reach(231, L), L <= 12, L>= 10, not reach(215, L+1), not reach(230, L+1), not reach(247, L+1).
reach(247,L + 1) :- reach(231, L), L <= 12, L>= 10, not reach(215, L+1), not reach(230, L+1), not reach(232, L+1).
reach(216,L + 1) :- reach(232, L), L <= 11, L>= 11, not reach(231, L+1), not reach(233, L+1), not reach(248, L+1).
reach(231,L + 1) :- reach(232, L), L <= 11, L>= 11, not reach(216, L+1), not reach(233, L+1), not reach(248, L+1).
reach(233,L + 1) :- reach(232, L), L <= 11, L>= 11, not reach(216, L+1), not reach(231, L+1), not reach(248, L+1).
reach(248,L + 1) :- reach(232, L), L <= 11, L>= 11, not reach(216, L+1), not reach(231, L+1), not reach(233, L+1).
reach(226,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(241, L+1), not reach(243, L+1).
reach(241,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(226, L+1), not reach(243, L+1).
reach(243,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(226, L+1), not reach(241, L+1).
reach(227,L + 1) :- reach(243, L), L <= 17, L>= 15, not reach(242, L+1), not reach(244, L+1).
reach(242,L + 1) :- reach(243, L), L <= 17, L>= 15, not reach(227, L+1), not reach(244, L+1).
reach(244,L + 1) :- reach(243, L), L <= 17, L>= 15, not reach(227, L+1), not reach(242, L+1).
reach(228,L + 1) :- reach(244, L), L <= 16, L>= 14, not reach(243, L+1), not reach(245, L+1).
reach(243,L + 1) :- reach(244, L), L <= 16, L>= 14, not reach(228, L+1), not reach(245, L+1).
reach(245,L + 1) :- reach(244, L), L <= 16, L>= 14, not reach(228, L+1), not reach(243, L+1).
reach(229,L + 1) :- reach(245, L), L <= 15, L>= 13, not reach(244, L+1), not reach(246, L+1).
reach(244,L + 1) :- reach(245, L), L <= 15, L>= 13, not reach(229, L+1), not reach(246, L+1).
reach(246,L + 1) :- reach(245, L), L <= 15, L>= 13, not reach(229, L+1), not reach(244, L+1).
reach(230,L + 1) :- reach(246, L), L <= 14, L>= 12, not reach(245, L+1), not reach(247, L+1).
reach(245,L + 1) :- reach(246, L), L <= 14, L>= 12, not reach(230, L+1), not reach(247, L+1).
reach(247,L + 1) :- reach(246, L), L <= 14, L>= 12, not reach(230, L+1), not reach(245, L+1).
reach(231,L + 1) :- reach(247, L), L <= 13, L>= 11, not reach(246, L+1), not reach(248, L+1).
reach(246,L + 1) :- reach(247, L), L <= 13, L>= 11, not reach(231, L+1), not reach(248, L+1).
reach(248,L + 1) :- reach(247, L), L <= 13, L>= 11, not reach(231, L+1), not reach(246, L+1).

