reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17), not reach(X, 18), not reach(X, 19), not reach(X, 20), not reach(X, 21), not reach(X, 22), not reach(X, 23), not reach(X, 24).
start(22).
goal(273).
:- reach(X, L), reach(X, L'), L != L'.
reach(4,L + 1) :- reach(2, L), L <= 10, L>= 7, not reach(5, L+1), not reach(3, L+1).
reach(5,L + 1) :- reach(2, L), L <= 10, L>= 7, not reach(4, L+1), not reach(3, L+1).
reach(3,L + 1) :- reach(2, L), L <= 10, L>= 7, not reach(4, L+1), not reach(5, L+1).
reach(5,L + 1) :- reach(3, L), L <= 10, L>= 8, not reach(6, L+1), not reach(2, L+1).
reach(6,L + 1) :- reach(3, L), L <= 10, L>= 8, not reach(5, L+1), not reach(2, L+1).
reach(2,L + 1) :- reach(3, L), L <= 10, L>= 8, not reach(5, L+1), not reach(6, L+1).
reach(2,L + 1) :- reach(4, L), L <= 10, L>= 6, not reach(7, L+1), not reach(8, L+1).
reach(7,L + 1) :- reach(4, L), L <= 10, L>= 6, not reach(2, L+1), not reach(8, L+1).
reach(8,L + 1) :- reach(4, L), L <= 10, L>= 6, not reach(2, L+1), not reach(7, L+1).
reach(2,L + 1) :- reach(5, L), L <= 11, L>= 8, not reach(3, L+1), not reach(8, L+1), not reach(9, L+1), not reach(12, L+1).
reach(3,L + 1) :- reach(5, L), L <= 11, L>= 8, not reach(2, L+1), not reach(8, L+1), not reach(9, L+1), not reach(12, L+1).
reach(8,L + 1) :- reach(5, L), L <= 11, L>= 8, not reach(2, L+1), not reach(3, L+1), not reach(9, L+1), not reach(12, L+1).
reach(9,L + 1) :- reach(5, L), L <= 11, L>= 8, not reach(2, L+1), not reach(3, L+1), not reach(8, L+1), not reach(12, L+1).
reach(12,L + 1) :- reach(5, L), L <= 11, L>= 8, not reach(2, L+1), not reach(3, L+1), not reach(8, L+1), not reach(9, L+1).
reach(3,L + 1) :- reach(6, L), L <= 10, L>= 9, not reach(9, L+1), not reach(52, L+1).
reach(9,L + 1) :- reach(6, L), L <= 10, L>= 9, not reach(3, L+1), not reach(52, L+1).
reach(52,L + 1) :- reach(6, L), L <= 10, L>= 9, not reach(3, L+1), not reach(9, L+1).
reach(4,L + 1) :- reach(7, L), L <= 9, L>= 5, not reach(10, L+1), not reach(11, L+1).
reach(10,L + 1) :- reach(7, L), L <= 9, L>= 5, not reach(4, L+1), not reach(11, L+1).
reach(11,L + 1) :- reach(7, L), L <= 9, L>= 5, not reach(4, L+1), not reach(10, L+1).
reach(4,L + 1) :- reach(8, L), L <= 11, L>= 7, not reach(5, L+1), not reach(11, L+1), not reach(12, L+1).
reach(5,L + 1) :- reach(8, L), L <= 11, L>= 7, not reach(4, L+1), not reach(11, L+1), not reach(12, L+1).
reach(11,L + 1) :- reach(8, L), L <= 11, L>= 7, not reach(4, L+1), not reach(5, L+1), not reach(12, L+1).
reach(12,L + 1) :- reach(8, L), L <= 11, L>= 7, not reach(4, L+1), not reach(5, L+1), not reach(11, L+1).
reach(5,L + 1) :- reach(9, L), L <= 11, L>= 9, not reach(6, L+1), not reach(12, L+1), not reach(53, L+1).
reach(6,L + 1) :- reach(9, L), L <= 11, L>= 9, not reach(5, L+1), not reach(12, L+1), not reach(53, L+1).
reach(12,L + 1) :- reach(9, L), L <= 11, L>= 9, not reach(5, L+1), not reach(6, L+1), not reach(53, L+1).
reach(53,L + 1) :- reach(9, L), L <= 11, L>= 9, not reach(5, L+1), not reach(6, L+1), not reach(12, L+1).
reach(5,L + 1) :- reach(12, L), L <= 12, L>= 8, not reach(8, L+1), not reach(9, L+1), not reach(15, L+1), not reach(54, L+1), not reach(140, L+1).
reach(8,L + 1) :- reach(12, L), L <= 12, L>= 8, not reach(5, L+1), not reach(9, L+1), not reach(15, L+1), not reach(54, L+1), not reach(140, L+1).
reach(9,L + 1) :- reach(12, L), L <= 12, L>= 8, not reach(5, L+1), not reach(8, L+1), not reach(15, L+1), not reach(54, L+1), not reach(140, L+1).
reach(15,L + 1) :- reach(12, L), L <= 12, L>= 8, not reach(5, L+1), not reach(8, L+1), not reach(9, L+1), not reach(54, L+1), not reach(140, L+1).
reach(54,L + 1) :- reach(12, L), L <= 12, L>= 8, not reach(5, L+1), not reach(8, L+1), not reach(9, L+1), not reach(15, L+1), not reach(140, L+1).
reach(140,L + 1) :- reach(12, L), L <= 12, L>= 8, not reach(5, L+1), not reach(8, L+1), not reach(9, L+1), not reach(15, L+1), not reach(54, L+1).
reach(6,L + 1) :- reach(52, L), L <= 11, L>= 10, not reach(53, L+1), not reach(69, L+1).
reach(53,L + 1) :- reach(52, L), L <= 11, L>= 10, not reach(6, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(52, L), L <= 11, L>= 10, not reach(6, L+1), not reach(53, L+1).
reach(7,L + 1) :- reach(10, L), L <= 8, L>= 4, not reach(13, L+1), not reach(14, L+1).
reach(13,L + 1) :- reach(10, L), L <= 8, L>= 4, not reach(7, L+1), not reach(14, L+1).
reach(14,L + 1) :- reach(10, L), L <= 8, L>= 4, not reach(7, L+1), not reach(13, L+1).
reach(7,L + 1) :- reach(11, L), L <= 10, L>= 6, not reach(8, L+1), not reach(14, L+1), not reach(15, L+1).
reach(8,L + 1) :- reach(11, L), L <= 10, L>= 6, not reach(7, L+1), not reach(14, L+1), not reach(15, L+1).
reach(14,L + 1) :- reach(11, L), L <= 10, L>= 6, not reach(7, L+1), not reach(8, L+1), not reach(15, L+1).
reach(15,L + 1) :- reach(11, L), L <= 10, L>= 6, not reach(7, L+1), not reach(8, L+1), not reach(14, L+1).
reach(9,L + 1) :- reach(53, L), L <= 10, L>= 10, not reach(52, L+1), not reach(54, L+1), not reach(70, L+1).
reach(52,L + 1) :- reach(53, L), L <= 10, L>= 10, not reach(9, L+1), not reach(54, L+1), not reach(70, L+1).
reach(54,L + 1) :- reach(53, L), L <= 10, L>= 10, not reach(9, L+1), not reach(52, L+1), not reach(70, L+1).
reach(70,L + 1) :- reach(53, L), L <= 10, L>= 10, not reach(9, L+1), not reach(52, L+1), not reach(54, L+1).
reach(10,L + 1) :- reach(13, L), L <= 7, L>= 3, not reach(16, L+1), not reach(17, L+1).
reach(16,L + 1) :- reach(13, L), L <= 7, L>= 3, not reach(10, L+1), not reach(17, L+1).
reach(17,L + 1) :- reach(13, L), L <= 7, L>= 3, not reach(10, L+1), not reach(16, L+1).
reach(10,L + 1) :- reach(14, L), L <= 9, L>= 5, not reach(11, L+1), not reach(17, L+1), not reach(18, L+1).
reach(11,L + 1) :- reach(14, L), L <= 9, L>= 5, not reach(10, L+1), not reach(17, L+1), not reach(18, L+1).
reach(17,L + 1) :- reach(14, L), L <= 9, L>= 5, not reach(10, L+1), not reach(11, L+1), not reach(18, L+1).
reach(18,L + 1) :- reach(14, L), L <= 9, L>= 5, not reach(10, L+1), not reach(11, L+1), not reach(17, L+1).
reach(11,L + 1) :- reach(15, L), L <= 11, L>= 7, not reach(12, L+1), not reach(18, L+1), not reach(55, L+1).
reach(12,L + 1) :- reach(15, L), L <= 11, L>= 7, not reach(11, L+1), not reach(18, L+1), not reach(55, L+1).
reach(18,L + 1) :- reach(15, L), L <= 11, L>= 7, not reach(11, L+1), not reach(12, L+1), not reach(55, L+1).
reach(55,L + 1) :- reach(15, L), L <= 11, L>= 7, not reach(11, L+1), not reach(12, L+1), not reach(18, L+1).
reach(12,L + 1) :- reach(54, L), L <= 11, L>= 9, not reach(53, L+1), not reach(55, L+1), not reach(71, L+1).
reach(53,L + 1) :- reach(54, L), L <= 11, L>= 9, not reach(12, L+1), not reach(55, L+1), not reach(71, L+1).
reach(55,L + 1) :- reach(54, L), L <= 11, L>= 9, not reach(12, L+1), not reach(53, L+1), not reach(71, L+1).
reach(71,L + 1) :- reach(54, L), L <= 11, L>= 9, not reach(12, L+1), not reach(53, L+1), not reach(55, L+1).
reach(12,L + 1) :- reach(140, L), L <= 13, L>= 9, not reach(123, L+1), not reach(139, L+1), not reach(141, L+1), not reach(157, L+1).
reach(123,L + 1) :- reach(140, L), L <= 13, L>= 9, not reach(12, L+1), not reach(139, L+1), not reach(141, L+1), not reach(157, L+1).
reach(139,L + 1) :- reach(140, L), L <= 13, L>= 9, not reach(12, L+1), not reach(123, L+1), not reach(141, L+1), not reach(157, L+1).
reach(141,L + 1) :- reach(140, L), L <= 13, L>= 9, not reach(12, L+1), not reach(123, L+1), not reach(139, L+1), not reach(157, L+1).
reach(157,L + 1) :- reach(140, L), L <= 13, L>= 9, not reach(12, L+1), not reach(123, L+1), not reach(139, L+1), not reach(141, L+1).
reach(13,L + 1) :- reach(16, L), L <= 6, L>= 2, not reach(19, L+1), not reach(20, L+1).
reach(19,L + 1) :- reach(16, L), L <= 6, L>= 2, not reach(13, L+1), not reach(20, L+1).
reach(20,L + 1) :- reach(16, L), L <= 6, L>= 2, not reach(13, L+1), not reach(19, L+1).
reach(13,L + 1) :- reach(17, L), L <= 8, L>= 4, not reach(14, L+1), not reach(20, L+1), not reach(21, L+1).
reach(14,L + 1) :- reach(17, L), L <= 8, L>= 4, not reach(13, L+1), not reach(20, L+1), not reach(21, L+1).
reach(20,L + 1) :- reach(17, L), L <= 8, L>= 4, not reach(13, L+1), not reach(14, L+1), not reach(21, L+1).
reach(21,L + 1) :- reach(17, L), L <= 8, L>= 4, not reach(13, L+1), not reach(14, L+1), not reach(20, L+1).
reach(14,L + 1) :- reach(18, L), L <= 10, L>= 6, not reach(15, L+1), not reach(21, L+1), not reach(56, L+1).
reach(15,L + 1) :- reach(18, L), L <= 10, L>= 6, not reach(14, L+1), not reach(21, L+1), not reach(56, L+1).
reach(21,L + 1) :- reach(18, L), L <= 10, L>= 6, not reach(14, L+1), not reach(15, L+1), not reach(56, L+1).
reach(56,L + 1) :- reach(18, L), L <= 10, L>= 6, not reach(14, L+1), not reach(15, L+1), not reach(21, L+1).
reach(15,L + 1) :- reach(55, L), L <= 10, L>= 8, not reach(54, L+1), not reach(72, L+1).
reach(54,L + 1) :- reach(55, L), L <= 10, L>= 8, not reach(15, L+1), not reach(72, L+1).
reach(72,L + 1) :- reach(55, L), L <= 10, L>= 8, not reach(15, L+1), not reach(54, L+1).
reach(16,L + 1) :- reach(19, L), L <= 5, L>= 1, not reach(22, L+1), not reach(23, L+1).
reach(23,L + 1) :- reach(19, L), L <= 5, L>= 1, not reach(16, L+1), not reach(22, L+1).
reach(16,L + 1) :- reach(20, L), L <= 7, L>= 3, not reach(17, L+1), not reach(23, L+1), not reach(24, L+1).
reach(17,L + 1) :- reach(20, L), L <= 7, L>= 3, not reach(16, L+1), not reach(23, L+1), not reach(24, L+1).
reach(23,L + 1) :- reach(20, L), L <= 7, L>= 3, not reach(16, L+1), not reach(17, L+1), not reach(24, L+1).
reach(24,L + 1) :- reach(20, L), L <= 7, L>= 3, not reach(16, L+1), not reach(17, L+1), not reach(23, L+1).
reach(17,L + 1) :- reach(21, L), L <= 9, L>= 5, not reach(18, L+1), not reach(24, L+1), not reach(57, L+1).
reach(18,L + 1) :- reach(21, L), L <= 9, L>= 5, not reach(17, L+1), not reach(24, L+1), not reach(57, L+1).
reach(24,L + 1) :- reach(21, L), L <= 9, L>= 5, not reach(17, L+1), not reach(18, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(21, L), L <= 9, L>= 5, not reach(17, L+1), not reach(18, L+1), not reach(24, L+1).
reach(18,L + 1) :- reach(56, L), L <= 9, L>= 7, not reach(57, L+1), not reach(73, L+1).
reach(57,L + 1) :- reach(56, L), L <= 9, L>= 7, not reach(18, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(56, L), L <= 9, L>= 7, not reach(18, L+1), not reach(57, L+1).
reach(19,L + 1) :- reach(22, L), L <= 5, L>= 0, not reach(25, L+1), not reach(26, L+1).
reach(25,L + 1) :- reach(22, L), L <= 5, L>= 0, not reach(19, L+1), not reach(26, L+1).
reach(26,L + 1) :- reach(22, L), L <= 5, L>= 0, not reach(19, L+1), not reach(25, L+1).
reach(19,L + 1) :- reach(23, L), L <= 6, L>= 2, not reach(20, L+1), not reach(26, L+1), not reach(27, L+1).
reach(20,L + 1) :- reach(23, L), L <= 6, L>= 2, not reach(19, L+1), not reach(26, L+1), not reach(27, L+1).
reach(26,L + 1) :- reach(23, L), L <= 6, L>= 2, not reach(19, L+1), not reach(20, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(23, L), L <= 6, L>= 2, not reach(19, L+1), not reach(20, L+1), not reach(26, L+1).
reach(20,L + 1) :- reach(24, L), L <= 8, L>= 4, not reach(21, L+1), not reach(27, L+1), not reach(58, L+1).
reach(21,L + 1) :- reach(24, L), L <= 8, L>= 4, not reach(20, L+1), not reach(27, L+1), not reach(58, L+1).
reach(27,L + 1) :- reach(24, L), L <= 8, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(24, L), L <= 8, L>= 4, not reach(20, L+1), not reach(21, L+1), not reach(27, L+1).
reach(21,L + 1) :- reach(57, L), L <= 10, L>= 6, not reach(56, L+1), not reach(58, L+1), not reach(74, L+1).
reach(56,L + 1) :- reach(57, L), L <= 10, L>= 6, not reach(21, L+1), not reach(58, L+1), not reach(74, L+1).
reach(58,L + 1) :- reach(57, L), L <= 10, L>= 6, not reach(21, L+1), not reach(56, L+1), not reach(74, L+1).
reach(74,L + 1) :- reach(57, L), L <= 10, L>= 6, not reach(21, L+1), not reach(56, L+1), not reach(58, L+1).
reach(28,L + 1) :- reach(25, L), L <= 6, L>= 1, not reach(22, L+1), not reach(29, L+1).
reach(29,L + 1) :- reach(25, L), L <= 6, L>= 1, not reach(22, L+1), not reach(28, L+1).
reach(23,L + 1) :- reach(26, L), L <= 5, L>= 1, not reach(22, L+1), not reach(29, L+1), not reach(30, L+1).
reach(29,L + 1) :- reach(26, L), L <= 5, L>= 1, not reach(22, L+1), not reach(23, L+1), not reach(30, L+1).
reach(30,L + 1) :- reach(26, L), L <= 5, L>= 1, not reach(22, L+1), not reach(23, L+1), not reach(29, L+1).
reach(23,L + 1) :- reach(27, L), L <= 7, L>= 3, not reach(24, L+1), not reach(30, L+1), not reach(59, L+1).
reach(24,L + 1) :- reach(27, L), L <= 7, L>= 3, not reach(23, L+1), not reach(30, L+1), not reach(59, L+1).
reach(30,L + 1) :- reach(27, L), L <= 7, L>= 3, not reach(23, L+1), not reach(24, L+1), not reach(59, L+1).
reach(59,L + 1) :- reach(27, L), L <= 7, L>= 3, not reach(23, L+1), not reach(24, L+1), not reach(30, L+1).
reach(24,L + 1) :- reach(58, L), L <= 9, L>= 5, not reach(57, L+1), not reach(59, L+1), not reach(75, L+1).
reach(57,L + 1) :- reach(58, L), L <= 9, L>= 5, not reach(24, L+1), not reach(59, L+1), not reach(75, L+1).
reach(59,L + 1) :- reach(58, L), L <= 9, L>= 5, not reach(24, L+1), not reach(57, L+1), not reach(75, L+1).
reach(75,L + 1) :- reach(58, L), L <= 9, L>= 5, not reach(24, L+1), not reach(57, L+1), not reach(59, L+1).
reach(25,L + 1) :- reach(28, L), L <= 7, L>= 2, not reach(31, L+1), not reach(32, L+1).
reach(31,L + 1) :- reach(28, L), L <= 7, L>= 2, not reach(25, L+1), not reach(32, L+1).
reach(32,L + 1) :- reach(28, L), L <= 7, L>= 2, not reach(25, L+1), not reach(31, L+1).
reach(25,L + 1) :- reach(29, L), L <= 6, L>= 2, not reach(26, L+1), not reach(32, L+1), not reach(33, L+1).
reach(26,L + 1) :- reach(29, L), L <= 6, L>= 2, not reach(25, L+1), not reach(32, L+1), not reach(33, L+1).
reach(32,L + 1) :- reach(29, L), L <= 6, L>= 2, not reach(25, L+1), not reach(26, L+1), not reach(33, L+1).
reach(33,L + 1) :- reach(29, L), L <= 6, L>= 2, not reach(25, L+1), not reach(26, L+1), not reach(32, L+1).
reach(26,L + 1) :- reach(30, L), L <= 6, L>= 2, not reach(27, L+1), not reach(33, L+1), not reach(60, L+1).
reach(27,L + 1) :- reach(30, L), L <= 6, L>= 2, not reach(26, L+1), not reach(33, L+1), not reach(60, L+1).
reach(33,L + 1) :- reach(30, L), L <= 6, L>= 2, not reach(26, L+1), not reach(27, L+1), not reach(60, L+1).
reach(60,L + 1) :- reach(30, L), L <= 6, L>= 2, not reach(26, L+1), not reach(27, L+1), not reach(33, L+1).
reach(27,L + 1) :- reach(59, L), L <= 8, L>= 4, not reach(58, L+1), not reach(60, L+1), not reach(76, L+1).
reach(58,L + 1) :- reach(59, L), L <= 8, L>= 4, not reach(27, L+1), not reach(60, L+1), not reach(76, L+1).
reach(60,L + 1) :- reach(59, L), L <= 8, L>= 4, not reach(27, L+1), not reach(58, L+1), not reach(76, L+1).
reach(76,L + 1) :- reach(59, L), L <= 8, L>= 4, not reach(27, L+1), not reach(58, L+1), not reach(60, L+1).
reach(28,L + 1) :- reach(31, L), L <= 8, L>= 3, not reach(35, L+1), not reach(38, L+1).
reach(35,L + 1) :- reach(31, L), L <= 8, L>= 3, not reach(28, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(31, L), L <= 8, L>= 3, not reach(28, L+1), not reach(35, L+1).
reach(28,L + 1) :- reach(32, L), L <= 7, L>= 3, not reach(29, L+1), not reach(35, L+1), not reach(36, L+1).
reach(29,L + 1) :- reach(32, L), L <= 7, L>= 3, not reach(28, L+1), not reach(35, L+1), not reach(36, L+1).
reach(35,L + 1) :- reach(32, L), L <= 7, L>= 3, not reach(28, L+1), not reach(29, L+1), not reach(36, L+1).
reach(36,L + 1) :- reach(32, L), L <= 7, L>= 3, not reach(28, L+1), not reach(29, L+1), not reach(35, L+1).
reach(29,L + 1) :- reach(33, L), L <= 7, L>= 3, not reach(30, L+1), not reach(36, L+1), not reach(61, L+1).
reach(30,L + 1) :- reach(33, L), L <= 7, L>= 3, not reach(29, L+1), not reach(36, L+1), not reach(61, L+1).
reach(36,L + 1) :- reach(33, L), L <= 7, L>= 3, not reach(29, L+1), not reach(30, L+1), not reach(61, L+1).
reach(61,L + 1) :- reach(33, L), L <= 7, L>= 3, not reach(29, L+1), not reach(30, L+1), not reach(36, L+1).
reach(30,L + 1) :- reach(60, L), L <= 7, L>= 3, not reach(59, L+1), not reach(61, L+1), not reach(77, L+1).
reach(59,L + 1) :- reach(60, L), L <= 7, L>= 3, not reach(30, L+1), not reach(61, L+1), not reach(77, L+1).
reach(61,L + 1) :- reach(60, L), L <= 7, L>= 3, not reach(30, L+1), not reach(59, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(60, L), L <= 7, L>= 3, not reach(30, L+1), not reach(59, L+1), not reach(61, L+1).
reach(31,L + 1) :- reach(35, L), L <= 8, L>= 4, not reach(32, L+1), not reach(38, L+1), not reach(39, L+1).
reach(32,L + 1) :- reach(35, L), L <= 8, L>= 4, not reach(31, L+1), not reach(38, L+1), not reach(39, L+1).
reach(38,L + 1) :- reach(35, L), L <= 8, L>= 4, not reach(31, L+1), not reach(32, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(35, L), L <= 8, L>= 4, not reach(31, L+1), not reach(32, L+1), not reach(38, L+1).
reach(32,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(33, L+1), not reach(39, L+1), not reach(62, L+1).
reach(33,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(32, L+1), not reach(39, L+1), not reach(62, L+1).
reach(39,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(32, L+1), not reach(33, L+1), not reach(62, L+1).
reach(62,L + 1) :- reach(36, L), L <= 8, L>= 4, not reach(32, L+1), not reach(33, L+1), not reach(39, L+1).
reach(33,L + 1) :- reach(61, L), L <= 8, L>= 4, not reach(60, L+1), not reach(62, L+1), not reach(78, L+1).
reach(60,L + 1) :- reach(61, L), L <= 8, L>= 4, not reach(33, L+1), not reach(62, L+1), not reach(78, L+1).
reach(62,L + 1) :- reach(61, L), L <= 8, L>= 4, not reach(33, L+1), not reach(60, L+1), not reach(78, L+1).
reach(78,L + 1) :- reach(61, L), L <= 8, L>= 4, not reach(33, L+1), not reach(60, L+1), not reach(62, L+1).
reach(35,L + 1) :- reach(38, L), L <= 9, L>= 4, not reach(41, L+1), not reach(42, L+1), not reach(31, L+1).
reach(41,L + 1) :- reach(38, L), L <= 9, L>= 4, not reach(35, L+1), not reach(42, L+1), not reach(31, L+1).
reach(42,L + 1) :- reach(38, L), L <= 9, L>= 4, not reach(35, L+1), not reach(41, L+1), not reach(31, L+1).
reach(31,L + 1) :- reach(38, L), L <= 9, L>= 4, not reach(35, L+1), not reach(41, L+1), not reach(42, L+1).
reach(35,L + 1) :- reach(39, L), L <= 7, L>= 5, not reach(36, L+1), not reach(42, L+1), not reach(63, L+1).
reach(36,L + 1) :- reach(39, L), L <= 7, L>= 5, not reach(35, L+1), not reach(42, L+1), not reach(63, L+1).
reach(42,L + 1) :- reach(39, L), L <= 7, L>= 5, not reach(35, L+1), not reach(36, L+1), not reach(63, L+1).
reach(63,L + 1) :- reach(39, L), L <= 7, L>= 5, not reach(35, L+1), not reach(36, L+1), not reach(42, L+1).
reach(36,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(61, L+1), not reach(63, L+1), not reach(79, L+1).
reach(61,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(36, L+1), not reach(63, L+1), not reach(79, L+1).
reach(63,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(36, L+1), not reach(61, L+1), not reach(79, L+1).
reach(79,L + 1) :- reach(62, L), L <= 9, L>= 5, not reach(36, L+1), not reach(61, L+1), not reach(63, L+1).
reach(40,L + 1) :- reach(37, L), L <= 11, L>= 6, not reach(41, L+1), not reach(163, L+1).
reach(41,L + 1) :- reach(37, L), L <= 11, L>= 6, not reach(40, L+1), not reach(163, L+1).
reach(163,L + 1) :- reach(37, L), L <= 11, L>= 6, not reach(40, L+1), not reach(41, L+1).
reach(37,L + 1) :- reach(40, L), L <= 10, L>= 7, not reach(43, L+1), not reach(44, L+1).
reach(43,L + 1) :- reach(40, L), L <= 10, L>= 7, not reach(37, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(40, L), L <= 10, L>= 7, not reach(37, L+1), not reach(43, L+1).
reach(37,L + 1) :- reach(41, L), L <= 10, L>= 5, not reach(38, L+1), not reach(44, L+1), not reach(45, L+1).
reach(38,L + 1) :- reach(41, L), L <= 10, L>= 5, not reach(37, L+1), not reach(44, L+1), not reach(45, L+1).
reach(44,L + 1) :- reach(41, L), L <= 10, L>= 5, not reach(37, L+1), not reach(38, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(41, L), L <= 10, L>= 5, not reach(37, L+1), not reach(38, L+1), not reach(44, L+1).
reach(37,L + 1) :- reach(163, L), L <= 12, L>= 7, not reach(146, L+1), not reach(162, L+1), not reach(164, L+1), not reach(180, L+1).
reach(146,L + 1) :- reach(163, L), L <= 12, L>= 7, not reach(37, L+1), not reach(162, L+1), not reach(164, L+1), not reach(180, L+1).
reach(162,L + 1) :- reach(163, L), L <= 12, L>= 7, not reach(37, L+1), not reach(146, L+1), not reach(164, L+1), not reach(180, L+1).
reach(164,L + 1) :- reach(163, L), L <= 12, L>= 7, not reach(37, L+1), not reach(146, L+1), not reach(162, L+1), not reach(180, L+1).
reach(180,L + 1) :- reach(163, L), L <= 12, L>= 7, not reach(37, L+1), not reach(146, L+1), not reach(162, L+1), not reach(164, L+1).
reach(38,L + 1) :- reach(42, L), L <= 8, L>= 5, not reach(39, L+1), not reach(45, L+1), not reach(64, L+1).
reach(39,L + 1) :- reach(42, L), L <= 8, L>= 5, not reach(38, L+1), not reach(45, L+1), not reach(64, L+1).
reach(45,L + 1) :- reach(42, L), L <= 8, L>= 5, not reach(38, L+1), not reach(39, L+1), not reach(64, L+1).
reach(64,L + 1) :- reach(42, L), L <= 8, L>= 5, not reach(38, L+1), not reach(39, L+1), not reach(45, L+1).
reach(39,L + 1) :- reach(63, L), L <= 8, L>= 6, not reach(62, L+1), not reach(64, L+1), not reach(80, L+1).
reach(62,L + 1) :- reach(63, L), L <= 8, L>= 6, not reach(39, L+1), not reach(64, L+1), not reach(80, L+1).
reach(64,L + 1) :- reach(63, L), L <= 8, L>= 6, not reach(39, L+1), not reach(62, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(63, L), L <= 8, L>= 6, not reach(39, L+1), not reach(62, L+1), not reach(64, L+1).
reach(40,L + 1) :- reach(43, L), L <= 9, L>= 8, not reach(47, L+1), not reach(49, L+1).
reach(47,L + 1) :- reach(43, L), L <= 9, L>= 8, not reach(40, L+1), not reach(49, L+1).
reach(49,L + 1) :- reach(43, L), L <= 9, L>= 8, not reach(40, L+1), not reach(47, L+1).
reach(40,L + 1) :- reach(44, L), L <= 9, L>= 6, not reach(41, L+1), not reach(47, L+1), not reach(48, L+1).
reach(41,L + 1) :- reach(44, L), L <= 9, L>= 6, not reach(40, L+1), not reach(47, L+1), not reach(48, L+1).
reach(47,L + 1) :- reach(44, L), L <= 9, L>= 6, not reach(40, L+1), not reach(41, L+1), not reach(48, L+1).
reach(48,L + 1) :- reach(44, L), L <= 9, L>= 6, not reach(40, L+1), not reach(41, L+1), not reach(47, L+1).
reach(41,L + 1) :- reach(45, L), L <= 9, L>= 6, not reach(42, L+1), not reach(48, L+1), not reach(65, L+1).
reach(42,L + 1) :- reach(45, L), L <= 9, L>= 6, not reach(41, L+1), not reach(48, L+1), not reach(65, L+1).
reach(48,L + 1) :- reach(45, L), L <= 9, L>= 6, not reach(41, L+1), not reach(42, L+1), not reach(65, L+1).
reach(65,L + 1) :- reach(45, L), L <= 9, L>= 6, not reach(41, L+1), not reach(42, L+1), not reach(48, L+1).
reach(42,L + 1) :- reach(64, L), L <= 9, L>= 6, not reach(63, L+1), not reach(65, L+1), not reach(81, L+1).
reach(63,L + 1) :- reach(64, L), L <= 9, L>= 6, not reach(42, L+1), not reach(65, L+1), not reach(81, L+1).
reach(65,L + 1) :- reach(64, L), L <= 9, L>= 6, not reach(42, L+1), not reach(63, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(64, L), L <= 9, L>= 6, not reach(42, L+1), not reach(63, L+1), not reach(65, L+1).
reach(43,L + 1) :- reach(47, L), L <= 8, L>= 7, not reach(44, L+1), not reach(49, L+1), not reach(50, L+1).
reach(44,L + 1) :- reach(47, L), L <= 8, L>= 7, not reach(43, L+1), not reach(49, L+1), not reach(50, L+1).
reach(49,L + 1) :- reach(47, L), L <= 8, L>= 7, not reach(43, L+1), not reach(44, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(47, L), L <= 8, L>= 7, not reach(43, L+1), not reach(44, L+1), not reach(49, L+1).
reach(44,L + 1) :- reach(48, L), L <= 8, L>= 7, not reach(45, L+1), not reach(50, L+1).
reach(45,L + 1) :- reach(48, L), L <= 8, L>= 7, not reach(44, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(48, L), L <= 8, L>= 7, not reach(44, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(65, L), L <= 10, L>= 7, not reach(64, L+1), not reach(66, L+1), not reach(82, L+1).
reach(64,L + 1) :- reach(65, L), L <= 10, L>= 7, not reach(45, L+1), not reach(66, L+1), not reach(82, L+1).
reach(66,L + 1) :- reach(65, L), L <= 10, L>= 7, not reach(45, L+1), not reach(64, L+1), not reach(82, L+1).
reach(82,L + 1) :- reach(65, L), L <= 10, L>= 7, not reach(45, L+1), not reach(64, L+1), not reach(66, L+1).
reach(47,L + 1) :- reach(49, L), L <= 8, L>= 8, not reach(51, L+1), not reach(43, L+1).
reach(51,L + 1) :- reach(49, L), L <= 8, L>= 8, not reach(47, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(49, L), L <= 8, L>= 8, not reach(47, L+1), not reach(51, L+1).
reach(52,L + 1) :- reach(69, L), L <= 12, L>= 11, not reach(70, L+1), not reach(86, L+1).
reach(70,L + 1) :- reach(69, L), L <= 12, L>= 11, not reach(52, L+1), not reach(86, L+1).
reach(86,L + 1) :- reach(69, L), L <= 12, L>= 11, not reach(52, L+1), not reach(70, L+1).
reach(53,L + 1) :- reach(70, L), L <= 11, L>= 11, not reach(69, L+1), not reach(71, L+1), not reach(87, L+1).
reach(69,L + 1) :- reach(70, L), L <= 11, L>= 11, not reach(53, L+1), not reach(71, L+1), not reach(87, L+1).
reach(71,L + 1) :- reach(70, L), L <= 11, L>= 11, not reach(53, L+1), not reach(69, L+1), not reach(87, L+1).
reach(87,L + 1) :- reach(70, L), L <= 11, L>= 11, not reach(53, L+1), not reach(69, L+1), not reach(71, L+1).
reach(54,L + 1) :- reach(71, L), L <= 10, L>= 10, not reach(70, L+1), not reach(72, L+1), not reach(88, L+1).
reach(70,L + 1) :- reach(71, L), L <= 10, L>= 10, not reach(54, L+1), not reach(72, L+1), not reach(88, L+1).
reach(72,L + 1) :- reach(71, L), L <= 10, L>= 10, not reach(54, L+1), not reach(70, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(71, L), L <= 10, L>= 10, not reach(54, L+1), not reach(70, L+1), not reach(72, L+1).
reach(55,L + 1) :- reach(72, L), L <= 9, L>= 9, not reach(71, L+1), not reach(73, L+1).
reach(71,L + 1) :- reach(72, L), L <= 9, L>= 9, not reach(55, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(72, L), L <= 9, L>= 9, not reach(55, L+1), not reach(71, L+1).
reach(56,L + 1) :- reach(73, L), L <= 10, L>= 8, not reach(72, L+1), not reach(74, L+1), not reach(89, L+1).
reach(72,L + 1) :- reach(73, L), L <= 10, L>= 8, not reach(56, L+1), not reach(74, L+1), not reach(89, L+1).
reach(74,L + 1) :- reach(73, L), L <= 10, L>= 8, not reach(56, L+1), not reach(72, L+1), not reach(89, L+1).
reach(89,L + 1) :- reach(73, L), L <= 10, L>= 8, not reach(56, L+1), not reach(72, L+1), not reach(74, L+1).
reach(57,L + 1) :- reach(74, L), L <= 11, L>= 7, not reach(73, L+1), not reach(75, L+1), not reach(90, L+1).
reach(73,L + 1) :- reach(74, L), L <= 11, L>= 7, not reach(57, L+1), not reach(75, L+1), not reach(90, L+1).
reach(75,L + 1) :- reach(74, L), L <= 11, L>= 7, not reach(57, L+1), not reach(73, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(74, L), L <= 11, L>= 7, not reach(57, L+1), not reach(73, L+1), not reach(75, L+1).
reach(58,L + 1) :- reach(75, L), L <= 10, L>= 6, not reach(74, L+1), not reach(76, L+1), not reach(91, L+1).
reach(74,L + 1) :- reach(75, L), L <= 10, L>= 6, not reach(58, L+1), not reach(76, L+1), not reach(91, L+1).
reach(76,L + 1) :- reach(75, L), L <= 10, L>= 6, not reach(58, L+1), not reach(74, L+1), not reach(91, L+1).
reach(91,L + 1) :- reach(75, L), L <= 10, L>= 6, not reach(58, L+1), not reach(74, L+1), not reach(76, L+1).
reach(59,L + 1) :- reach(76, L), L <= 9, L>= 5, not reach(75, L+1), not reach(77, L+1), not reach(92, L+1).
reach(75,L + 1) :- reach(76, L), L <= 9, L>= 5, not reach(59, L+1), not reach(77, L+1), not reach(92, L+1).
reach(77,L + 1) :- reach(76, L), L <= 9, L>= 5, not reach(59, L+1), not reach(75, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(76, L), L <= 9, L>= 5, not reach(59, L+1), not reach(75, L+1), not reach(77, L+1).
reach(60,L + 1) :- reach(77, L), L <= 8, L>= 4, not reach(76, L+1), not reach(78, L+1), not reach(93, L+1).
reach(76,L + 1) :- reach(77, L), L <= 8, L>= 4, not reach(60, L+1), not reach(78, L+1), not reach(93, L+1).
reach(78,L + 1) :- reach(77, L), L <= 8, L>= 4, not reach(60, L+1), not reach(76, L+1), not reach(93, L+1).
reach(93,L + 1) :- reach(77, L), L <= 8, L>= 4, not reach(60, L+1), not reach(76, L+1), not reach(78, L+1).
reach(61,L + 1) :- reach(78, L), L <= 9, L>= 5, not reach(77, L+1), not reach(79, L+1), not reach(94, L+1).
reach(77,L + 1) :- reach(78, L), L <= 9, L>= 5, not reach(61, L+1), not reach(79, L+1), not reach(94, L+1).
reach(79,L + 1) :- reach(78, L), L <= 9, L>= 5, not reach(61, L+1), not reach(77, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(78, L), L <= 9, L>= 5, not reach(61, L+1), not reach(77, L+1), not reach(79, L+1).
reach(62,L + 1) :- reach(79, L), L <= 10, L>= 6, not reach(78, L+1), not reach(80, L+1), not reach(95, L+1).
reach(78,L + 1) :- reach(79, L), L <= 10, L>= 6, not reach(62, L+1), not reach(80, L+1), not reach(95, L+1).
reach(80,L + 1) :- reach(79, L), L <= 10, L>= 6, not reach(62, L+1), not reach(78, L+1), not reach(95, L+1).
reach(95,L + 1) :- reach(79, L), L <= 10, L>= 6, not reach(62, L+1), not reach(78, L+1), not reach(80, L+1).
reach(63,L + 1) :- reach(80, L), L <= 9, L>= 7, not reach(79, L+1), not reach(81, L+1), not reach(96, L+1).
reach(79,L + 1) :- reach(80, L), L <= 9, L>= 7, not reach(63, L+1), not reach(81, L+1), not reach(96, L+1).
reach(81,L + 1) :- reach(80, L), L <= 9, L>= 7, not reach(63, L+1), not reach(79, L+1), not reach(96, L+1).
reach(96,L + 1) :- reach(80, L), L <= 9, L>= 7, not reach(63, L+1), not reach(79, L+1), not reach(81, L+1).
reach(64,L + 1) :- reach(81, L), L <= 10, L>= 7, not reach(80, L+1), not reach(82, L+1), not reach(97, L+1).
reach(80,L + 1) :- reach(81, L), L <= 10, L>= 7, not reach(64, L+1), not reach(82, L+1), not reach(97, L+1).
reach(82,L + 1) :- reach(81, L), L <= 10, L>= 7, not reach(64, L+1), not reach(80, L+1), not reach(97, L+1).
reach(97,L + 1) :- reach(81, L), L <= 10, L>= 7, not reach(64, L+1), not reach(80, L+1), not reach(82, L+1).
reach(65,L + 1) :- reach(66, L), L <= 9, L>= 8, not reach(67, L+1), not reach(83, L+1).
reach(67,L + 1) :- reach(66, L), L <= 9, L>= 8, not reach(65, L+1), not reach(83, L+1).
reach(83,L + 1) :- reach(66, L), L <= 9, L>= 8, not reach(65, L+1), not reach(67, L+1).
reach(65,L + 1) :- reach(82, L), L <= 11, L>= 8, not reach(81, L+1), not reach(83, L+1), not reach(98, L+1), not reach(113, L+1).
reach(81,L + 1) :- reach(82, L), L <= 11, L>= 8, not reach(65, L+1), not reach(83, L+1), not reach(98, L+1), not reach(113, L+1).
reach(83,L + 1) :- reach(82, L), L <= 11, L>= 8, not reach(65, L+1), not reach(81, L+1), not reach(98, L+1), not reach(113, L+1).
reach(98,L + 1) :- reach(82, L), L <= 11, L>= 8, not reach(65, L+1), not reach(81, L+1), not reach(83, L+1), not reach(113, L+1).
reach(113,L + 1) :- reach(82, L), L <= 11, L>= 8, not reach(65, L+1), not reach(81, L+1), not reach(83, L+1), not reach(98, L+1).
reach(66,L + 1) :- reach(83, L), L <= 10, L>= 9, not reach(82, L+1), not reach(84, L+1), not reach(99, L+1).
reach(82,L + 1) :- reach(83, L), L <= 10, L>= 9, not reach(66, L+1), not reach(84, L+1), not reach(99, L+1).
reach(84,L + 1) :- reach(83, L), L <= 10, L>= 9, not reach(66, L+1), not reach(82, L+1), not reach(99, L+1).
reach(99,L + 1) :- reach(83, L), L <= 10, L>= 9, not reach(66, L+1), not reach(82, L+1), not reach(84, L+1).
reach(69,L + 1) :- reach(86, L), L <= 13, L>= 12, not reach(87, L+1), not reach(102, L+1), not reach(135, L+1).
reach(87,L + 1) :- reach(86, L), L <= 13, L>= 12, not reach(69, L+1), not reach(102, L+1), not reach(135, L+1).
reach(102,L + 1) :- reach(86, L), L <= 13, L>= 12, not reach(69, L+1), not reach(87, L+1), not reach(135, L+1).
reach(135,L + 1) :- reach(86, L), L <= 13, L>= 12, not reach(69, L+1), not reach(87, L+1), not reach(102, L+1).
reach(70,L + 1) :- reach(87, L), L <= 12, L>= 12, not reach(86, L+1), not reach(88, L+1), not reach(103, L+1).
reach(86,L + 1) :- reach(87, L), L <= 12, L>= 12, not reach(70, L+1), not reach(88, L+1), not reach(103, L+1).
reach(88,L + 1) :- reach(87, L), L <= 12, L>= 12, not reach(70, L+1), not reach(86, L+1), not reach(103, L+1).
reach(103,L + 1) :- reach(87, L), L <= 12, L>= 12, not reach(70, L+1), not reach(86, L+1), not reach(88, L+1).
reach(71,L + 1) :- reach(88, L), L <= 11, L>= 11, not reach(87, L+1), not reach(104, L+1), not reach(105, L+1).
reach(87,L + 1) :- reach(88, L), L <= 11, L>= 11, not reach(71, L+1), not reach(104, L+1), not reach(105, L+1).
reach(104,L + 1) :- reach(88, L), L <= 11, L>= 11, not reach(71, L+1), not reach(87, L+1), not reach(105, L+1).
reach(105,L + 1) :- reach(88, L), L <= 11, L>= 11, not reach(71, L+1), not reach(87, L+1), not reach(104, L+1).
reach(73,L + 1) :- reach(89, L), L <= 11, L>= 9, not reach(104, L+1), not reach(90, L+1), not reach(107, L+1).
reach(104,L + 1) :- reach(89, L), L <= 11, L>= 9, not reach(73, L+1), not reach(90, L+1), not reach(107, L+1).
reach(90,L + 1) :- reach(89, L), L <= 11, L>= 9, not reach(73, L+1), not reach(104, L+1), not reach(107, L+1).
reach(107,L + 1) :- reach(89, L), L <= 11, L>= 9, not reach(73, L+1), not reach(104, L+1), not reach(90, L+1).
reach(74,L + 1) :- reach(90, L), L <= 12, L>= 8, not reach(89, L+1), not reach(91, L+1), not reach(108, L+1).
reach(89,L + 1) :- reach(90, L), L <= 12, L>= 8, not reach(74, L+1), not reach(91, L+1), not reach(108, L+1).
reach(91,L + 1) :- reach(90, L), L <= 12, L>= 8, not reach(74, L+1), not reach(89, L+1), not reach(108, L+1).
reach(108,L + 1) :- reach(90, L), L <= 12, L>= 8, not reach(74, L+1), not reach(89, L+1), not reach(91, L+1).
reach(75,L + 1) :- reach(91, L), L <= 11, L>= 7, not reach(90, L+1), not reach(92, L+1), not reach(109, L+1).
reach(90,L + 1) :- reach(91, L), L <= 11, L>= 7, not reach(75, L+1), not reach(92, L+1), not reach(109, L+1).
reach(92,L + 1) :- reach(91, L), L <= 11, L>= 7, not reach(75, L+1), not reach(90, L+1), not reach(109, L+1).
reach(109,L + 1) :- reach(91, L), L <= 11, L>= 7, not reach(75, L+1), not reach(90, L+1), not reach(92, L+1).
reach(76,L + 1) :- reach(92, L), L <= 10, L>= 6, not reach(91, L+1), not reach(93, L+1), not reach(110, L+1).
reach(91,L + 1) :- reach(92, L), L <= 10, L>= 6, not reach(76, L+1), not reach(93, L+1), not reach(110, L+1).
reach(93,L + 1) :- reach(92, L), L <= 10, L>= 6, not reach(76, L+1), not reach(91, L+1), not reach(110, L+1).
reach(110,L + 1) :- reach(92, L), L <= 10, L>= 6, not reach(76, L+1), not reach(91, L+1), not reach(93, L+1).
reach(77,L + 1) :- reach(93, L), L <= 9, L>= 5, not reach(92, L+1), not reach(94, L+1), not reach(111, L+1).
reach(92,L + 1) :- reach(93, L), L <= 9, L>= 5, not reach(77, L+1), not reach(94, L+1), not reach(111, L+1).
reach(94,L + 1) :- reach(93, L), L <= 9, L>= 5, not reach(77, L+1), not reach(92, L+1), not reach(111, L+1).
reach(111,L + 1) :- reach(93, L), L <= 9, L>= 5, not reach(77, L+1), not reach(92, L+1), not reach(94, L+1).
reach(78,L + 1) :- reach(94, L), L <= 10, L>= 6, not reach(93, L+1), not reach(95, L+1), not reach(112, L+1).
reach(93,L + 1) :- reach(94, L), L <= 10, L>= 6, not reach(78, L+1), not reach(95, L+1), not reach(112, L+1).
reach(95,L + 1) :- reach(94, L), L <= 10, L>= 6, not reach(78, L+1), not reach(93, L+1), not reach(112, L+1).
reach(112,L + 1) :- reach(94, L), L <= 10, L>= 6, not reach(78, L+1), not reach(93, L+1), not reach(95, L+1).
reach(79,L + 1) :- reach(95, L), L <= 11, L>= 7, not reach(94, L+1), not reach(96, L+1), not reach(113, L+1).
reach(94,L + 1) :- reach(95, L), L <= 11, L>= 7, not reach(79, L+1), not reach(96, L+1), not reach(113, L+1).
reach(96,L + 1) :- reach(95, L), L <= 11, L>= 7, not reach(79, L+1), not reach(94, L+1), not reach(113, L+1).
reach(113,L + 1) :- reach(95, L), L <= 11, L>= 7, not reach(79, L+1), not reach(94, L+1), not reach(96, L+1).
reach(80,L + 1) :- reach(96, L), L <= 10, L>= 8, not reach(95, L+1), not reach(97, L+1), not reach(114, L+1).
reach(95,L + 1) :- reach(96, L), L <= 10, L>= 8, not reach(80, L+1), not reach(97, L+1), not reach(114, L+1).
reach(97,L + 1) :- reach(96, L), L <= 10, L>= 8, not reach(80, L+1), not reach(95, L+1), not reach(114, L+1).
reach(114,L + 1) :- reach(96, L), L <= 10, L>= 8, not reach(80, L+1), not reach(95, L+1), not reach(97, L+1).
reach(81,L + 1) :- reach(97, L), L <= 9, L>= 8, not reach(96, L+1), not reach(98, L+1), not reach(115, L+1).
reach(96,L + 1) :- reach(97, L), L <= 9, L>= 8, not reach(81, L+1), not reach(98, L+1), not reach(115, L+1).
reach(98,L + 1) :- reach(97, L), L <= 9, L>= 8, not reach(81, L+1), not reach(96, L+1), not reach(115, L+1).
reach(115,L + 1) :- reach(97, L), L <= 9, L>= 8, not reach(81, L+1), not reach(96, L+1), not reach(98, L+1).
reach(82,L + 1) :- reach(98, L), L <= 10, L>= 9, not reach(97, L+1), not reach(99, L+1), not reach(116, L+1).
reach(97,L + 1) :- reach(98, L), L <= 10, L>= 9, not reach(82, L+1), not reach(99, L+1), not reach(116, L+1).
reach(99,L + 1) :- reach(98, L), L <= 10, L>= 9, not reach(82, L+1), not reach(97, L+1), not reach(116, L+1).
reach(116,L + 1) :- reach(98, L), L <= 10, L>= 9, not reach(82, L+1), not reach(97, L+1), not reach(99, L+1).
reach(82,L + 1) :- reach(113, L), L <= 12, L>= 8, not reach(95, L+1), not reach(112, L+1), not reach(114, L+1), not reach(130, L+1).
reach(95,L + 1) :- reach(113, L), L <= 12, L>= 8, not reach(82, L+1), not reach(112, L+1), not reach(114, L+1), not reach(130, L+1).
reach(112,L + 1) :- reach(113, L), L <= 12, L>= 8, not reach(82, L+1), not reach(95, L+1), not reach(114, L+1), not reach(130, L+1).
reach(114,L + 1) :- reach(113, L), L <= 12, L>= 8, not reach(82, L+1), not reach(95, L+1), not reach(112, L+1), not reach(130, L+1).
reach(130,L + 1) :- reach(113, L), L <= 12, L>= 8, not reach(82, L+1), not reach(95, L+1), not reach(112, L+1), not reach(114, L+1).
reach(86,L + 1) :- reach(102, L), L <= 14, L>= 13, not reach(103, L+1), not reach(120, L+1).
reach(103,L + 1) :- reach(102, L), L <= 14, L>= 13, not reach(86, L+1), not reach(120, L+1).
reach(120,L + 1) :- reach(102, L), L <= 14, L>= 13, not reach(86, L+1), not reach(103, L+1).
reach(87,L + 1) :- reach(103, L), L <= 13, L>= 13, not reach(102, L+1), not reach(105, L+1), not reach(121, L+1).
reach(102,L + 1) :- reach(103, L), L <= 13, L>= 13, not reach(87, L+1), not reach(105, L+1), not reach(121, L+1).
reach(105,L + 1) :- reach(103, L), L <= 13, L>= 13, not reach(87, L+1), not reach(102, L+1), not reach(121, L+1).
reach(121,L + 1) :- reach(103, L), L <= 13, L>= 13, not reach(87, L+1), not reach(102, L+1), not reach(105, L+1).
reach(88,L + 1) :- reach(104, L), L <= 10, L>= 10, not reach(89, L+1), not reach(106, L+1).
reach(89,L + 1) :- reach(104, L), L <= 10, L>= 10, not reach(88, L+1), not reach(106, L+1).
reach(106,L + 1) :- reach(104, L), L <= 10, L>= 10, not reach(88, L+1), not reach(89, L+1).
reach(88,L + 1) :- reach(105, L), L <= 12, L>= 12, not reach(103, L+1), not reach(106, L+1), not reach(122, L+1).
reach(103,L + 1) :- reach(105, L), L <= 12, L>= 12, not reach(88, L+1), not reach(106, L+1), not reach(122, L+1).
reach(106,L + 1) :- reach(105, L), L <= 12, L>= 12, not reach(88, L+1), not reach(103, L+1), not reach(122, L+1).
reach(122,L + 1) :- reach(105, L), L <= 12, L>= 12, not reach(88, L+1), not reach(103, L+1), not reach(106, L+1).
reach(89,L + 1) :- reach(107, L), L <= 12, L>= 10, not reach(106, L+1), not reach(108, L+1), not reach(124, L+1).
reach(106,L + 1) :- reach(107, L), L <= 12, L>= 10, not reach(89, L+1), not reach(108, L+1), not reach(124, L+1).
reach(108,L + 1) :- reach(107, L), L <= 12, L>= 10, not reach(89, L+1), not reach(106, L+1), not reach(124, L+1).
reach(124,L + 1) :- reach(107, L), L <= 12, L>= 10, not reach(89, L+1), not reach(106, L+1), not reach(108, L+1).
reach(90,L + 1) :- reach(108, L), L <= 13, L>= 9, not reach(107, L+1), not reach(109, L+1), not reach(125, L+1).
reach(107,L + 1) :- reach(108, L), L <= 13, L>= 9, not reach(90, L+1), not reach(109, L+1), not reach(125, L+1).
reach(109,L + 1) :- reach(108, L), L <= 13, L>= 9, not reach(90, L+1), not reach(107, L+1), not reach(125, L+1).
reach(125,L + 1) :- reach(108, L), L <= 13, L>= 9, not reach(90, L+1), not reach(107, L+1), not reach(109, L+1).
reach(91,L + 1) :- reach(109, L), L <= 12, L>= 8, not reach(108, L+1), not reach(110, L+1), not reach(126, L+1).
reach(108,L + 1) :- reach(109, L), L <= 12, L>= 8, not reach(91, L+1), not reach(110, L+1), not reach(126, L+1).
reach(110,L + 1) :- reach(109, L), L <= 12, L>= 8, not reach(91, L+1), not reach(108, L+1), not reach(126, L+1).
reach(126,L + 1) :- reach(109, L), L <= 12, L>= 8, not reach(91, L+1), not reach(108, L+1), not reach(110, L+1).
reach(92,L + 1) :- reach(110, L), L <= 11, L>= 7, not reach(109, L+1), not reach(111, L+1), not reach(127, L+1).
reach(109,L + 1) :- reach(110, L), L <= 11, L>= 7, not reach(92, L+1), not reach(111, L+1), not reach(127, L+1).
reach(111,L + 1) :- reach(110, L), L <= 11, L>= 7, not reach(92, L+1), not reach(109, L+1), not reach(127, L+1).
reach(127,L + 1) :- reach(110, L), L <= 11, L>= 7, not reach(92, L+1), not reach(109, L+1), not reach(111, L+1).
reach(93,L + 1) :- reach(111, L), L <= 10, L>= 6, not reach(110, L+1), not reach(112, L+1), not reach(128, L+1).
reach(110,L + 1) :- reach(111, L), L <= 10, L>= 6, not reach(93, L+1), not reach(112, L+1), not reach(128, L+1).
reach(112,L + 1) :- reach(111, L), L <= 10, L>= 6, not reach(93, L+1), not reach(110, L+1), not reach(128, L+1).
reach(128,L + 1) :- reach(111, L), L <= 10, L>= 6, not reach(93, L+1), not reach(110, L+1), not reach(112, L+1).
reach(94,L + 1) :- reach(112, L), L <= 11, L>= 7, not reach(111, L+1), not reach(113, L+1), not reach(129, L+1).
reach(111,L + 1) :- reach(112, L), L <= 11, L>= 7, not reach(94, L+1), not reach(113, L+1), not reach(129, L+1).
reach(113,L + 1) :- reach(112, L), L <= 11, L>= 7, not reach(94, L+1), not reach(111, L+1), not reach(129, L+1).
reach(129,L + 1) :- reach(112, L), L <= 11, L>= 7, not reach(94, L+1), not reach(111, L+1), not reach(113, L+1).
reach(96,L + 1) :- reach(114, L), L <= 11, L>= 9, not reach(113, L+1), not reach(115, L+1), not reach(131, L+1).
reach(113,L + 1) :- reach(114, L), L <= 11, L>= 9, not reach(96, L+1), not reach(115, L+1), not reach(131, L+1).
reach(115,L + 1) :- reach(114, L), L <= 11, L>= 9, not reach(96, L+1), not reach(113, L+1), not reach(131, L+1).
reach(131,L + 1) :- reach(114, L), L <= 11, L>= 9, not reach(96, L+1), not reach(113, L+1), not reach(115, L+1).
reach(97,L + 1) :- reach(115, L), L <= 10, L>= 9, not reach(114, L+1), not reach(116, L+1), not reach(132, L+1).
reach(114,L + 1) :- reach(115, L), L <= 10, L>= 9, not reach(97, L+1), not reach(116, L+1), not reach(132, L+1).
reach(116,L + 1) :- reach(115, L), L <= 10, L>= 9, not reach(97, L+1), not reach(114, L+1), not reach(132, L+1).
reach(132,L + 1) :- reach(115, L), L <= 10, L>= 9, not reach(97, L+1), not reach(114, L+1), not reach(116, L+1).
reach(102,L + 1) :- reach(120, L), L <= 15, L>= 13, not reach(121, L+1), not reach(137, L+1).
reach(121,L + 1) :- reach(120, L), L <= 15, L>= 13, not reach(102, L+1), not reach(137, L+1).
reach(137,L + 1) :- reach(120, L), L <= 15, L>= 13, not reach(102, L+1), not reach(121, L+1).
reach(103,L + 1) :- reach(121, L), L <= 14, L>= 12, not reach(120, L+1), not reach(122, L+1), not reach(138, L+1).
reach(120,L + 1) :- reach(121, L), L <= 14, L>= 12, not reach(103, L+1), not reach(122, L+1), not reach(138, L+1).
reach(122,L + 1) :- reach(121, L), L <= 14, L>= 12, not reach(103, L+1), not reach(120, L+1), not reach(138, L+1).
reach(138,L + 1) :- reach(121, L), L <= 14, L>= 12, not reach(103, L+1), not reach(120, L+1), not reach(122, L+1).
reach(104,L + 1) :- reach(106, L), L <= 11, L>= 11, not reach(105, L+1), not reach(107, L+1), not reach(123, L+1).
reach(105,L + 1) :- reach(106, L), L <= 11, L>= 11, not reach(104, L+1), not reach(107, L+1), not reach(123, L+1).
reach(107,L + 1) :- reach(106, L), L <= 11, L>= 11, not reach(104, L+1), not reach(105, L+1), not reach(123, L+1).
reach(123,L + 1) :- reach(106, L), L <= 11, L>= 11, not reach(104, L+1), not reach(105, L+1), not reach(107, L+1).
reach(105,L + 1) :- reach(122, L), L <= 13, L>= 11, not reach(121, L+1), not reach(123, L+1), not reach(139, L+1).
reach(121,L + 1) :- reach(122, L), L <= 13, L>= 11, not reach(105, L+1), not reach(123, L+1), not reach(139, L+1).
reach(123,L + 1) :- reach(122, L), L <= 13, L>= 11, not reach(105, L+1), not reach(121, L+1), not reach(139, L+1).
reach(139,L + 1) :- reach(122, L), L <= 13, L>= 11, not reach(105, L+1), not reach(121, L+1), not reach(123, L+1).
reach(106,L + 1) :- reach(123, L), L <= 12, L>= 10, not reach(122, L+1), not reach(124, L+1), not reach(140, L+1).
reach(122,L + 1) :- reach(123, L), L <= 12, L>= 10, not reach(106, L+1), not reach(124, L+1), not reach(140, L+1).
reach(124,L + 1) :- reach(123, L), L <= 12, L>= 10, not reach(106, L+1), not reach(122, L+1), not reach(140, L+1).
reach(140,L + 1) :- reach(123, L), L <= 12, L>= 10, not reach(106, L+1), not reach(122, L+1), not reach(124, L+1).
reach(107,L + 1) :- reach(124, L), L <= 13, L>= 11, not reach(123, L+1), not reach(125, L+1), not reach(141, L+1).
reach(123,L + 1) :- reach(124, L), L <= 13, L>= 11, not reach(107, L+1), not reach(125, L+1), not reach(141, L+1).
reach(125,L + 1) :- reach(124, L), L <= 13, L>= 11, not reach(107, L+1), not reach(123, L+1), not reach(141, L+1).
reach(141,L + 1) :- reach(124, L), L <= 13, L>= 11, not reach(107, L+1), not reach(123, L+1), not reach(125, L+1).
reach(108,L + 1) :- reach(125, L), L <= 14, L>= 10, not reach(124, L+1), not reach(126, L+1), not reach(142, L+1), not reach(251, L+1).
reach(124,L + 1) :- reach(125, L), L <= 14, L>= 10, not reach(108, L+1), not reach(126, L+1), not reach(142, L+1), not reach(251, L+1).
reach(126,L + 1) :- reach(125, L), L <= 14, L>= 10, not reach(108, L+1), not reach(124, L+1), not reach(142, L+1), not reach(251, L+1).
reach(142,L + 1) :- reach(125, L), L <= 14, L>= 10, not reach(108, L+1), not reach(124, L+1), not reach(126, L+1), not reach(251, L+1).
reach(251,L + 1) :- reach(125, L), L <= 14, L>= 10, not reach(108, L+1), not reach(124, L+1), not reach(126, L+1), not reach(142, L+1).
reach(109,L + 1) :- reach(126, L), L <= 13, L>= 9, not reach(125, L+1), not reach(127, L+1), not reach(143, L+1).
reach(125,L + 1) :- reach(126, L), L <= 13, L>= 9, not reach(109, L+1), not reach(127, L+1), not reach(143, L+1).
reach(127,L + 1) :- reach(126, L), L <= 13, L>= 9, not reach(109, L+1), not reach(125, L+1), not reach(143, L+1).
reach(143,L + 1) :- reach(126, L), L <= 13, L>= 9, not reach(109, L+1), not reach(125, L+1), not reach(127, L+1).
reach(110,L + 1) :- reach(127, L), L <= 12, L>= 8, not reach(126, L+1), not reach(128, L+1), not reach(144, L+1).
reach(126,L + 1) :- reach(127, L), L <= 12, L>= 8, not reach(110, L+1), not reach(128, L+1), not reach(144, L+1).
reach(128,L + 1) :- reach(127, L), L <= 12, L>= 8, not reach(110, L+1), not reach(126, L+1), not reach(144, L+1).
reach(144,L + 1) :- reach(127, L), L <= 12, L>= 8, not reach(110, L+1), not reach(126, L+1), not reach(128, L+1).
reach(111,L + 1) :- reach(128, L), L <= 11, L>= 7, not reach(127, L+1), not reach(129, L+1), not reach(145, L+1).
reach(127,L + 1) :- reach(128, L), L <= 11, L>= 7, not reach(111, L+1), not reach(129, L+1), not reach(145, L+1).
reach(129,L + 1) :- reach(128, L), L <= 11, L>= 7, not reach(111, L+1), not reach(127, L+1), not reach(145, L+1).
reach(145,L + 1) :- reach(128, L), L <= 11, L>= 7, not reach(111, L+1), not reach(127, L+1), not reach(129, L+1).
reach(112,L + 1) :- reach(129, L), L <= 12, L>= 8, not reach(128, L+1), not reach(130, L+1), not reach(146, L+1).
reach(128,L + 1) :- reach(129, L), L <= 12, L>= 8, not reach(112, L+1), not reach(130, L+1), not reach(146, L+1).
reach(130,L + 1) :- reach(129, L), L <= 12, L>= 8, not reach(112, L+1), not reach(128, L+1), not reach(146, L+1).
reach(146,L + 1) :- reach(129, L), L <= 12, L>= 8, not reach(112, L+1), not reach(128, L+1), not reach(130, L+1).
reach(113,L + 1) :- reach(130, L), L <= 13, L>= 9, not reach(129, L+1), not reach(131, L+1), not reach(147, L+1).
reach(129,L + 1) :- reach(130, L), L <= 13, L>= 9, not reach(113, L+1), not reach(131, L+1), not reach(147, L+1).
reach(131,L + 1) :- reach(130, L), L <= 13, L>= 9, not reach(113, L+1), not reach(129, L+1), not reach(147, L+1).
reach(147,L + 1) :- reach(130, L), L <= 13, L>= 9, not reach(113, L+1), not reach(129, L+1), not reach(131, L+1).
reach(114,L + 1) :- reach(131, L), L <= 12, L>= 10, not reach(130, L+1), not reach(132, L+1), not reach(148, L+1).
reach(130,L + 1) :- reach(131, L), L <= 12, L>= 10, not reach(114, L+1), not reach(132, L+1), not reach(148, L+1).
reach(132,L + 1) :- reach(131, L), L <= 12, L>= 10, not reach(114, L+1), not reach(130, L+1), not reach(148, L+1).
reach(148,L + 1) :- reach(131, L), L <= 12, L>= 10, not reach(114, L+1), not reach(130, L+1), not reach(132, L+1).
reach(115,L + 1) :- reach(132, L), L <= 11, L>= 10, not reach(131, L+1), not reach(133, L+1), not reach(149, L+1).
reach(131,L + 1) :- reach(132, L), L <= 11, L>= 10, not reach(115, L+1), not reach(133, L+1), not reach(149, L+1).
reach(133,L + 1) :- reach(132, L), L <= 11, L>= 10, not reach(115, L+1), not reach(131, L+1), not reach(149, L+1).
reach(149,L + 1) :- reach(132, L), L <= 11, L>= 10, not reach(115, L+1), not reach(131, L+1), not reach(133, L+1).
reach(120,L + 1) :- reach(137, L), L <= 16, L>= 12, not reach(138, L+1), not reach(154, L+1).
reach(138,L + 1) :- reach(137, L), L <= 16, L>= 12, not reach(120, L+1), not reach(154, L+1).
reach(154,L + 1) :- reach(137, L), L <= 16, L>= 12, not reach(120, L+1), not reach(138, L+1).
reach(121,L + 1) :- reach(138, L), L <= 15, L>= 11, not reach(137, L+1), not reach(139, L+1), not reach(155, L+1).
reach(137,L + 1) :- reach(138, L), L <= 15, L>= 11, not reach(121, L+1), not reach(139, L+1), not reach(155, L+1).
reach(139,L + 1) :- reach(138, L), L <= 15, L>= 11, not reach(121, L+1), not reach(137, L+1), not reach(155, L+1).
reach(155,L + 1) :- reach(138, L), L <= 15, L>= 11, not reach(121, L+1), not reach(137, L+1), not reach(139, L+1).
reach(122,L + 1) :- reach(139, L), L <= 14, L>= 10, not reach(138, L+1), not reach(140, L+1), not reach(156, L+1).
reach(138,L + 1) :- reach(139, L), L <= 14, L>= 10, not reach(122, L+1), not reach(140, L+1), not reach(156, L+1).
reach(140,L + 1) :- reach(139, L), L <= 14, L>= 10, not reach(122, L+1), not reach(138, L+1), not reach(156, L+1).
reach(156,L + 1) :- reach(139, L), L <= 14, L>= 10, not reach(122, L+1), not reach(138, L+1), not reach(140, L+1).
reach(124,L + 1) :- reach(141, L), L <= 12, L>= 10, not reach(140, L+1), not reach(142, L+1), not reach(158, L+1).
reach(140,L + 1) :- reach(141, L), L <= 12, L>= 10, not reach(124, L+1), not reach(142, L+1), not reach(158, L+1).
reach(142,L + 1) :- reach(141, L), L <= 12, L>= 10, not reach(124, L+1), not reach(140, L+1), not reach(158, L+1).
reach(158,L + 1) :- reach(141, L), L <= 12, L>= 10, not reach(124, L+1), not reach(140, L+1), not reach(142, L+1).
reach(125,L + 1) :- reach(142, L), L <= 13, L>= 11, not reach(141, L+1), not reach(143, L+1), not reach(159, L+1).
reach(141,L + 1) :- reach(142, L), L <= 13, L>= 11, not reach(125, L+1), not reach(143, L+1), not reach(159, L+1).
reach(143,L + 1) :- reach(142, L), L <= 13, L>= 11, not reach(125, L+1), not reach(141, L+1), not reach(159, L+1).
reach(159,L + 1) :- reach(142, L), L <= 13, L>= 11, not reach(125, L+1), not reach(141, L+1), not reach(143, L+1).
reach(125,L + 1) :- reach(251, L), L <= 15, L>= 11, not reach(233, L+1), not reach(252, L+1), not reach(268, L+1).
reach(233,L + 1) :- reach(251, L), L <= 15, L>= 11, not reach(125, L+1), not reach(252, L+1), not reach(268, L+1).
reach(252,L + 1) :- reach(251, L), L <= 15, L>= 11, not reach(125, L+1), not reach(233, L+1), not reach(268, L+1).
reach(268,L + 1) :- reach(251, L), L <= 15, L>= 11, not reach(125, L+1), not reach(233, L+1), not reach(252, L+1).
reach(126,L + 1) :- reach(143, L), L <= 12, L>= 10, not reach(142, L+1), not reach(160, L+1).
reach(142,L + 1) :- reach(143, L), L <= 12, L>= 10, not reach(126, L+1), not reach(160, L+1).
reach(160,L + 1) :- reach(143, L), L <= 12, L>= 10, not reach(126, L+1), not reach(142, L+1).
reach(127,L + 1) :- reach(144, L), L <= 11, L>= 9, not reach(145, L+1), not reach(161, L+1).
reach(145,L + 1) :- reach(144, L), L <= 11, L>= 9, not reach(127, L+1), not reach(161, L+1).
reach(161,L + 1) :- reach(144, L), L <= 11, L>= 9, not reach(127, L+1), not reach(145, L+1).
reach(128,L + 1) :- reach(145, L), L <= 12, L>= 8, not reach(144, L+1), not reach(146, L+1), not reach(162, L+1).
reach(144,L + 1) :- reach(145, L), L <= 12, L>= 8, not reach(128, L+1), not reach(146, L+1), not reach(162, L+1).
reach(146,L + 1) :- reach(145, L), L <= 12, L>= 8, not reach(128, L+1), not reach(144, L+1), not reach(162, L+1).
reach(162,L + 1) :- reach(145, L), L <= 12, L>= 8, not reach(128, L+1), not reach(144, L+1), not reach(146, L+1).
reach(129,L + 1) :- reach(146, L), L <= 13, L>= 8, not reach(145, L+1), not reach(147, L+1), not reach(163, L+1).
reach(145,L + 1) :- reach(146, L), L <= 13, L>= 8, not reach(129, L+1), not reach(147, L+1), not reach(163, L+1).
reach(147,L + 1) :- reach(146, L), L <= 13, L>= 8, not reach(129, L+1), not reach(145, L+1), not reach(163, L+1).
reach(163,L + 1) :- reach(146, L), L <= 13, L>= 8, not reach(129, L+1), not reach(145, L+1), not reach(147, L+1).
reach(130,L + 1) :- reach(147, L), L <= 14, L>= 9, not reach(146, L+1), not reach(148, L+1), not reach(164, L+1), not reach(236, L+1).
reach(146,L + 1) :- reach(147, L), L <= 14, L>= 9, not reach(130, L+1), not reach(148, L+1), not reach(164, L+1), not reach(236, L+1).
reach(148,L + 1) :- reach(147, L), L <= 14, L>= 9, not reach(130, L+1), not reach(146, L+1), not reach(164, L+1), not reach(236, L+1).
reach(164,L + 1) :- reach(147, L), L <= 14, L>= 9, not reach(130, L+1), not reach(146, L+1), not reach(148, L+1), not reach(236, L+1).
reach(236,L + 1) :- reach(147, L), L <= 14, L>= 9, not reach(130, L+1), not reach(146, L+1), not reach(148, L+1), not reach(164, L+1).
reach(131,L + 1) :- reach(148, L), L <= 13, L>= 10, not reach(147, L+1), not reach(149, L+1), not reach(165, L+1).
reach(147,L + 1) :- reach(148, L), L <= 13, L>= 10, not reach(131, L+1), not reach(149, L+1), not reach(165, L+1).
reach(149,L + 1) :- reach(148, L), L <= 13, L>= 10, not reach(131, L+1), not reach(147, L+1), not reach(165, L+1).
reach(165,L + 1) :- reach(148, L), L <= 13, L>= 10, not reach(131, L+1), not reach(147, L+1), not reach(149, L+1).
reach(132,L + 1) :- reach(149, L), L <= 12, L>= 11, not reach(148, L+1), not reach(150, L+1), not reach(166, L+1).
reach(148,L + 1) :- reach(149, L), L <= 12, L>= 11, not reach(132, L+1), not reach(150, L+1), not reach(166, L+1).
reach(150,L + 1) :- reach(149, L), L <= 12, L>= 11, not reach(132, L+1), not reach(148, L+1), not reach(166, L+1).
reach(166,L + 1) :- reach(149, L), L <= 12, L>= 11, not reach(132, L+1), not reach(148, L+1), not reach(150, L+1).
reach(137,L + 1) :- reach(154, L), L <= 17, L>= 13, not reach(155, L+1), not reach(171, L+1).
reach(155,L + 1) :- reach(154, L), L <= 17, L>= 13, not reach(137, L+1), not reach(171, L+1).
reach(171,L + 1) :- reach(154, L), L <= 17, L>= 13, not reach(137, L+1), not reach(155, L+1).
reach(138,L + 1) :- reach(155, L), L <= 16, L>= 12, not reach(154, L+1), not reach(156, L+1), not reach(172, L+1).
reach(154,L + 1) :- reach(155, L), L <= 16, L>= 12, not reach(138, L+1), not reach(156, L+1), not reach(172, L+1).
reach(156,L + 1) :- reach(155, L), L <= 16, L>= 12, not reach(138, L+1), not reach(154, L+1), not reach(172, L+1).
reach(172,L + 1) :- reach(155, L), L <= 16, L>= 12, not reach(138, L+1), not reach(154, L+1), not reach(156, L+1).
reach(139,L + 1) :- reach(156, L), L <= 15, L>= 11, not reach(155, L+1), not reach(157, L+1), not reach(173, L+1).
reach(155,L + 1) :- reach(156, L), L <= 15, L>= 11, not reach(139, L+1), not reach(157, L+1), not reach(173, L+1).
reach(157,L + 1) :- reach(156, L), L <= 15, L>= 11, not reach(139, L+1), not reach(155, L+1), not reach(173, L+1).
reach(173,L + 1) :- reach(156, L), L <= 15, L>= 11, not reach(139, L+1), not reach(155, L+1), not reach(157, L+1).
reach(140,L + 1) :- reach(157, L), L <= 14, L>= 10, not reach(156, L+1), not reach(158, L+1), not reach(174, L+1).
reach(156,L + 1) :- reach(157, L), L <= 14, L>= 10, not reach(140, L+1), not reach(158, L+1), not reach(174, L+1).
reach(158,L + 1) :- reach(157, L), L <= 14, L>= 10, not reach(140, L+1), not reach(156, L+1), not reach(174, L+1).
reach(174,L + 1) :- reach(157, L), L <= 14, L>= 10, not reach(140, L+1), not reach(156, L+1), not reach(158, L+1).
reach(141,L + 1) :- reach(158, L), L <= 13, L>= 11, not reach(157, L+1), not reach(159, L+1), not reach(175, L+1).
reach(157,L + 1) :- reach(158, L), L <= 13, L>= 11, not reach(141, L+1), not reach(159, L+1), not reach(175, L+1).
reach(159,L + 1) :- reach(158, L), L <= 13, L>= 11, not reach(141, L+1), not reach(157, L+1), not reach(175, L+1).
reach(175,L + 1) :- reach(158, L), L <= 13, L>= 11, not reach(141, L+1), not reach(157, L+1), not reach(159, L+1).
reach(142,L + 1) :- reach(159, L), L <= 12, L>= 11, not reach(158, L+1), not reach(160, L+1), not reach(176, L+1).
reach(158,L + 1) :- reach(159, L), L <= 12, L>= 11, not reach(142, L+1), not reach(160, L+1), not reach(176, L+1).
reach(160,L + 1) :- reach(159, L), L <= 12, L>= 11, not reach(142, L+1), not reach(158, L+1), not reach(176, L+1).
reach(176,L + 1) :- reach(159, L), L <= 12, L>= 11, not reach(142, L+1), not reach(158, L+1), not reach(160, L+1).
reach(143,L + 1) :- reach(160, L), L <= 11, L>= 10, not reach(159, L+1), not reach(161, L+1).
reach(159,L + 1) :- reach(160, L), L <= 11, L>= 10, not reach(143, L+1), not reach(161, L+1).
reach(161,L + 1) :- reach(160, L), L <= 11, L>= 10, not reach(143, L+1), not reach(159, L+1).
reach(144,L + 1) :- reach(161, L), L <= 11, L>= 9, not reach(160, L+1), not reach(162, L+1), not reach(178, L+1).
reach(160,L + 1) :- reach(161, L), L <= 11, L>= 9, not reach(144, L+1), not reach(162, L+1), not reach(178, L+1).
reach(162,L + 1) :- reach(161, L), L <= 11, L>= 9, not reach(144, L+1), not reach(160, L+1), not reach(178, L+1).
reach(178,L + 1) :- reach(161, L), L <= 11, L>= 9, not reach(144, L+1), not reach(160, L+1), not reach(162, L+1).
reach(145,L + 1) :- reach(162, L), L <= 12, L>= 8, not reach(161, L+1), not reach(163, L+1), not reach(179, L+1).
reach(161,L + 1) :- reach(162, L), L <= 12, L>= 8, not reach(145, L+1), not reach(163, L+1), not reach(179, L+1).
reach(163,L + 1) :- reach(162, L), L <= 12, L>= 8, not reach(145, L+1), not reach(161, L+1), not reach(179, L+1).
reach(179,L + 1) :- reach(162, L), L <= 12, L>= 8, not reach(145, L+1), not reach(161, L+1), not reach(163, L+1).
reach(147,L + 1) :- reach(164, L), L <= 13, L>= 8, not reach(163, L+1), not reach(165, L+1), not reach(181, L+1).
reach(163,L + 1) :- reach(164, L), L <= 13, L>= 8, not reach(147, L+1), not reach(165, L+1), not reach(181, L+1).
reach(165,L + 1) :- reach(164, L), L <= 13, L>= 8, not reach(147, L+1), not reach(163, L+1), not reach(181, L+1).
reach(181,L + 1) :- reach(164, L), L <= 13, L>= 8, not reach(147, L+1), not reach(163, L+1), not reach(165, L+1).
reach(147,L + 1) :- reach(236, L), L <= 15, L>= 10, not reach(220, L+1), not reach(235, L+1), not reach(237, L+1), not reach(254, L+1).
reach(220,L + 1) :- reach(236, L), L <= 15, L>= 10, not reach(147, L+1), not reach(235, L+1), not reach(237, L+1), not reach(254, L+1).
reach(235,L + 1) :- reach(236, L), L <= 15, L>= 10, not reach(147, L+1), not reach(220, L+1), not reach(237, L+1), not reach(254, L+1).
reach(237,L + 1) :- reach(236, L), L <= 15, L>= 10, not reach(147, L+1), not reach(220, L+1), not reach(235, L+1), not reach(254, L+1).
reach(254,L + 1) :- reach(236, L), L <= 15, L>= 10, not reach(147, L+1), not reach(220, L+1), not reach(235, L+1), not reach(237, L+1).
reach(148,L + 1) :- reach(165, L), L <= 12, L>= 9, not reach(164, L+1), not reach(166, L+1), not reach(182, L+1).
reach(164,L + 1) :- reach(165, L), L <= 12, L>= 9, not reach(148, L+1), not reach(166, L+1), not reach(182, L+1).
reach(166,L + 1) :- reach(165, L), L <= 12, L>= 9, not reach(148, L+1), not reach(164, L+1), not reach(182, L+1).
reach(182,L + 1) :- reach(165, L), L <= 12, L>= 9, not reach(148, L+1), not reach(164, L+1), not reach(166, L+1).
reach(149,L + 1) :- reach(166, L), L <= 13, L>= 10, not reach(165, L+1), not reach(167, L+1), not reach(183, L+1).
reach(165,L + 1) :- reach(166, L), L <= 13, L>= 10, not reach(149, L+1), not reach(167, L+1), not reach(183, L+1).
reach(167,L + 1) :- reach(166, L), L <= 13, L>= 10, not reach(149, L+1), not reach(165, L+1), not reach(183, L+1).
reach(183,L + 1) :- reach(166, L), L <= 13, L>= 10, not reach(149, L+1), not reach(165, L+1), not reach(167, L+1).
reach(150,L + 1) :- reach(167, L), L <= 12, L>= 11, not reach(166, L+1), not reach(168, L+1), not reach(184, L+1).
reach(166,L + 1) :- reach(167, L), L <= 12, L>= 11, not reach(150, L+1), not reach(168, L+1), not reach(184, L+1).
reach(168,L + 1) :- reach(167, L), L <= 12, L>= 11, not reach(150, L+1), not reach(166, L+1), not reach(184, L+1).
reach(184,L + 1) :- reach(167, L), L <= 12, L>= 11, not reach(150, L+1), not reach(166, L+1), not reach(168, L+1).
reach(151,L + 1) :- reach(168, L), L <= 12, L>= 12, not reach(167, L+1), not reach(169, L+1), not reach(185, L+1).
reach(167,L + 1) :- reach(168, L), L <= 12, L>= 12, not reach(151, L+1), not reach(169, L+1), not reach(185, L+1).
reach(169,L + 1) :- reach(168, L), L <= 12, L>= 12, not reach(151, L+1), not reach(167, L+1), not reach(185, L+1).
reach(185,L + 1) :- reach(168, L), L <= 12, L>= 12, not reach(151, L+1), not reach(167, L+1), not reach(169, L+1).
reach(154,L + 1) :- reach(171, L), L <= 18, L>= 14, not reach(172, L+1), not reach(188, L+1).
reach(172,L + 1) :- reach(171, L), L <= 18, L>= 14, not reach(154, L+1), not reach(188, L+1).
reach(188,L + 1) :- reach(171, L), L <= 18, L>= 14, not reach(154, L+1), not reach(172, L+1).
reach(155,L + 1) :- reach(172, L), L <= 17, L>= 13, not reach(171, L+1), not reach(173, L+1), not reach(189, L+1).
reach(171,L + 1) :- reach(172, L), L <= 17, L>= 13, not reach(155, L+1), not reach(173, L+1), not reach(189, L+1).
reach(173,L + 1) :- reach(172, L), L <= 17, L>= 13, not reach(155, L+1), not reach(171, L+1), not reach(189, L+1).
reach(189,L + 1) :- reach(172, L), L <= 17, L>= 13, not reach(155, L+1), not reach(171, L+1), not reach(173, L+1).
reach(156,L + 1) :- reach(173, L), L <= 16, L>= 12, not reach(172, L+1), not reach(174, L+1), not reach(190, L+1).
reach(172,L + 1) :- reach(173, L), L <= 16, L>= 12, not reach(156, L+1), not reach(174, L+1), not reach(190, L+1).
reach(174,L + 1) :- reach(173, L), L <= 16, L>= 12, not reach(156, L+1), not reach(172, L+1), not reach(190, L+1).
reach(190,L + 1) :- reach(173, L), L <= 16, L>= 12, not reach(156, L+1), not reach(172, L+1), not reach(174, L+1).
reach(157,L + 1) :- reach(174, L), L <= 15, L>= 11, not reach(173, L+1), not reach(175, L+1), not reach(191, L+1).
reach(173,L + 1) :- reach(174, L), L <= 15, L>= 11, not reach(157, L+1), not reach(175, L+1), not reach(191, L+1).
reach(175,L + 1) :- reach(174, L), L <= 15, L>= 11, not reach(157, L+1), not reach(173, L+1), not reach(191, L+1).
reach(191,L + 1) :- reach(174, L), L <= 15, L>= 11, not reach(157, L+1), not reach(173, L+1), not reach(175, L+1).
reach(158,L + 1) :- reach(175, L), L <= 14, L>= 12, not reach(174, L+1), not reach(176, L+1), not reach(192, L+1).
reach(174,L + 1) :- reach(175, L), L <= 14, L>= 12, not reach(158, L+1), not reach(176, L+1), not reach(192, L+1).
reach(176,L + 1) :- reach(175, L), L <= 14, L>= 12, not reach(158, L+1), not reach(174, L+1), not reach(192, L+1).
reach(192,L + 1) :- reach(175, L), L <= 14, L>= 12, not reach(158, L+1), not reach(174, L+1), not reach(176, L+1).
reach(159,L + 1) :- reach(176, L), L <= 13, L>= 12, not reach(175, L+1), not reach(177, L+1), not reach(193, L+1).
reach(175,L + 1) :- reach(176, L), L <= 13, L>= 12, not reach(159, L+1), not reach(177, L+1), not reach(193, L+1).
reach(177,L + 1) :- reach(176, L), L <= 13, L>= 12, not reach(159, L+1), not reach(175, L+1), not reach(193, L+1).
reach(193,L + 1) :- reach(176, L), L <= 13, L>= 12, not reach(159, L+1), not reach(175, L+1), not reach(177, L+1).
reach(161,L + 1) :- reach(178, L), L <= 12, L>= 10, not reach(177, L+1), not reach(179, L+1), not reach(195, L+1).
reach(177,L + 1) :- reach(178, L), L <= 12, L>= 10, not reach(161, L+1), not reach(179, L+1), not reach(195, L+1).
reach(179,L + 1) :- reach(178, L), L <= 12, L>= 10, not reach(161, L+1), not reach(177, L+1), not reach(195, L+1).
reach(195,L + 1) :- reach(178, L), L <= 12, L>= 10, not reach(161, L+1), not reach(177, L+1), not reach(179, L+1).
reach(162,L + 1) :- reach(179, L), L <= 13, L>= 9, not reach(178, L+1), not reach(180, L+1), not reach(196, L+1).
reach(178,L + 1) :- reach(179, L), L <= 13, L>= 9, not reach(162, L+1), not reach(180, L+1), not reach(196, L+1).
reach(180,L + 1) :- reach(179, L), L <= 13, L>= 9, not reach(162, L+1), not reach(178, L+1), not reach(196, L+1).
reach(196,L + 1) :- reach(179, L), L <= 13, L>= 9, not reach(162, L+1), not reach(178, L+1), not reach(180, L+1).
reach(163,L + 1) :- reach(180, L), L <= 12, L>= 8, not reach(179, L+1), not reach(181, L+1), not reach(197, L+1).
reach(179,L + 1) :- reach(180, L), L <= 12, L>= 8, not reach(163, L+1), not reach(181, L+1), not reach(197, L+1).
reach(181,L + 1) :- reach(180, L), L <= 12, L>= 8, not reach(163, L+1), not reach(179, L+1), not reach(197, L+1).
reach(197,L + 1) :- reach(180, L), L <= 12, L>= 8, not reach(163, L+1), not reach(179, L+1), not reach(181, L+1).
reach(164,L + 1) :- reach(181, L), L <= 12, L>= 9, not reach(180, L+1), not reach(182, L+1), not reach(198, L+1).
reach(180,L + 1) :- reach(181, L), L <= 12, L>= 9, not reach(164, L+1), not reach(182, L+1), not reach(198, L+1).
reach(182,L + 1) :- reach(181, L), L <= 12, L>= 9, not reach(164, L+1), not reach(180, L+1), not reach(198, L+1).
reach(198,L + 1) :- reach(181, L), L <= 12, L>= 9, not reach(164, L+1), not reach(180, L+1), not reach(182, L+1).
reach(165,L + 1) :- reach(182, L), L <= 13, L>= 10, not reach(181, L+1), not reach(183, L+1), not reach(199, L+1).
reach(181,L + 1) :- reach(182, L), L <= 13, L>= 10, not reach(165, L+1), not reach(183, L+1), not reach(199, L+1).
reach(183,L + 1) :- reach(182, L), L <= 13, L>= 10, not reach(165, L+1), not reach(181, L+1), not reach(199, L+1).
reach(199,L + 1) :- reach(182, L), L <= 13, L>= 10, not reach(165, L+1), not reach(181, L+1), not reach(183, L+1).
reach(166,L + 1) :- reach(183, L), L <= 14, L>= 11, not reach(182, L+1), not reach(184, L+1), not reach(200, L+1).
reach(182,L + 1) :- reach(183, L), L <= 14, L>= 11, not reach(166, L+1), not reach(184, L+1), not reach(200, L+1).
reach(184,L + 1) :- reach(183, L), L <= 14, L>= 11, not reach(166, L+1), not reach(182, L+1), not reach(200, L+1).
reach(200,L + 1) :- reach(183, L), L <= 14, L>= 11, not reach(166, L+1), not reach(182, L+1), not reach(184, L+1).
reach(167,L + 1) :- reach(184, L), L <= 13, L>= 12, not reach(183, L+1), not reach(185, L+1), not reach(201, L+1).
reach(183,L + 1) :- reach(184, L), L <= 13, L>= 12, not reach(167, L+1), not reach(185, L+1), not reach(201, L+1).
reach(185,L + 1) :- reach(184, L), L <= 13, L>= 12, not reach(167, L+1), not reach(183, L+1), not reach(201, L+1).
reach(201,L + 1) :- reach(184, L), L <= 13, L>= 12, not reach(167, L+1), not reach(183, L+1), not reach(185, L+1).
reach(168,L + 1) :- reach(185, L), L <= 13, L>= 13, not reach(184, L+1), not reach(186, L+1), not reach(202, L+1).
reach(184,L + 1) :- reach(185, L), L <= 13, L>= 13, not reach(168, L+1), not reach(186, L+1), not reach(202, L+1).
reach(186,L + 1) :- reach(185, L), L <= 13, L>= 13, not reach(168, L+1), not reach(184, L+1), not reach(202, L+1).
reach(202,L + 1) :- reach(185, L), L <= 13, L>= 13, not reach(168, L+1), not reach(184, L+1), not reach(186, L+1).
reach(171,L + 1) :- reach(188, L), L <= 19, L>= 14, not reach(189, L+1), not reach(205, L+1), not reach(270, L+1).
reach(189,L + 1) :- reach(188, L), L <= 19, L>= 14, not reach(171, L+1), not reach(205, L+1), not reach(270, L+1).
reach(205,L + 1) :- reach(188, L), L <= 19, L>= 14, not reach(171, L+1), not reach(189, L+1), not reach(270, L+1).
reach(270,L + 1) :- reach(188, L), L <= 19, L>= 14, not reach(171, L+1), not reach(189, L+1), not reach(205, L+1).
reach(172,L + 1) :- reach(189, L), L <= 18, L>= 14, not reach(188, L+1), not reach(190, L+1), not reach(206, L+1).
reach(188,L + 1) :- reach(189, L), L <= 18, L>= 14, not reach(172, L+1), not reach(190, L+1), not reach(206, L+1).
reach(190,L + 1) :- reach(189, L), L <= 18, L>= 14, not reach(172, L+1), not reach(188, L+1), not reach(206, L+1).
reach(206,L + 1) :- reach(189, L), L <= 18, L>= 14, not reach(172, L+1), not reach(188, L+1), not reach(190, L+1).
reach(173,L + 1) :- reach(190, L), L <= 17, L>= 13, not reach(189, L+1), not reach(191, L+1), not reach(207, L+1).
reach(189,L + 1) :- reach(190, L), L <= 17, L>= 13, not reach(173, L+1), not reach(191, L+1), not reach(207, L+1).
reach(191,L + 1) :- reach(190, L), L <= 17, L>= 13, not reach(173, L+1), not reach(189, L+1), not reach(207, L+1).
reach(207,L + 1) :- reach(190, L), L <= 17, L>= 13, not reach(173, L+1), not reach(189, L+1), not reach(191, L+1).
reach(174,L + 1) :- reach(191, L), L <= 16, L>= 12, not reach(190, L+1), not reach(192, L+1), not reach(208, L+1).
reach(190,L + 1) :- reach(191, L), L <= 16, L>= 12, not reach(174, L+1), not reach(192, L+1), not reach(208, L+1).
reach(192,L + 1) :- reach(191, L), L <= 16, L>= 12, not reach(174, L+1), not reach(190, L+1), not reach(208, L+1).
reach(208,L + 1) :- reach(191, L), L <= 16, L>= 12, not reach(174, L+1), not reach(190, L+1), not reach(192, L+1).
reach(175,L + 1) :- reach(192, L), L <= 15, L>= 13, not reach(191, L+1), not reach(193, L+1), not reach(209, L+1).
reach(191,L + 1) :- reach(192, L), L <= 15, L>= 13, not reach(175, L+1), not reach(193, L+1), not reach(209, L+1).
reach(193,L + 1) :- reach(192, L), L <= 15, L>= 13, not reach(175, L+1), not reach(191, L+1), not reach(209, L+1).
reach(209,L + 1) :- reach(192, L), L <= 15, L>= 13, not reach(175, L+1), not reach(191, L+1), not reach(193, L+1).
reach(176,L + 1) :- reach(177, L), L <= 13, L>= 11, not reach(178, L+1), not reach(194, L+1).
reach(178,L + 1) :- reach(177, L), L <= 13, L>= 11, not reach(176, L+1), not reach(194, L+1).
reach(194,L + 1) :- reach(177, L), L <= 13, L>= 11, not reach(176, L+1), not reach(178, L+1).
reach(176,L + 1) :- reach(193, L), L <= 14, L>= 13, not reach(192, L+1), not reach(194, L+1), not reach(210, L+1).
reach(192,L + 1) :- reach(193, L), L <= 14, L>= 13, not reach(176, L+1), not reach(194, L+1), not reach(210, L+1).
reach(194,L + 1) :- reach(193, L), L <= 14, L>= 13, not reach(176, L+1), not reach(192, L+1), not reach(210, L+1).
reach(210,L + 1) :- reach(193, L), L <= 14, L>= 13, not reach(176, L+1), not reach(192, L+1), not reach(194, L+1).
reach(177,L + 1) :- reach(194, L), L <= 14, L>= 12, not reach(193, L+1), not reach(195, L+1), not reach(211, L+1).
reach(193,L + 1) :- reach(194, L), L <= 14, L>= 12, not reach(177, L+1), not reach(195, L+1), not reach(211, L+1).
reach(195,L + 1) :- reach(194, L), L <= 14, L>= 12, not reach(177, L+1), not reach(193, L+1), not reach(211, L+1).
reach(211,L + 1) :- reach(194, L), L <= 14, L>= 12, not reach(177, L+1), not reach(193, L+1), not reach(195, L+1).
reach(178,L + 1) :- reach(195, L), L <= 13, L>= 11, not reach(194, L+1), not reach(196, L+1), not reach(212, L+1).
reach(194,L + 1) :- reach(195, L), L <= 13, L>= 11, not reach(178, L+1), not reach(196, L+1), not reach(212, L+1).
reach(196,L + 1) :- reach(195, L), L <= 13, L>= 11, not reach(178, L+1), not reach(194, L+1), not reach(212, L+1).
reach(212,L + 1) :- reach(195, L), L <= 13, L>= 11, not reach(178, L+1), not reach(194, L+1), not reach(196, L+1).
reach(179,L + 1) :- reach(196, L), L <= 14, L>= 10, not reach(195, L+1), not reach(197, L+1), not reach(213, L+1).
reach(195,L + 1) :- reach(196, L), L <= 14, L>= 10, not reach(179, L+1), not reach(197, L+1), not reach(213, L+1).
reach(197,L + 1) :- reach(196, L), L <= 14, L>= 10, not reach(179, L+1), not reach(195, L+1), not reach(213, L+1).
reach(213,L + 1) :- reach(196, L), L <= 14, L>= 10, not reach(179, L+1), not reach(195, L+1), not reach(197, L+1).
reach(180,L + 1) :- reach(197, L), L <= 13, L>= 9, not reach(196, L+1), not reach(198, L+1), not reach(214, L+1).
reach(196,L + 1) :- reach(197, L), L <= 13, L>= 9, not reach(180, L+1), not reach(198, L+1), not reach(214, L+1).
reach(198,L + 1) :- reach(197, L), L <= 13, L>= 9, not reach(180, L+1), not reach(196, L+1), not reach(214, L+1).
reach(214,L + 1) :- reach(197, L), L <= 13, L>= 9, not reach(180, L+1), not reach(196, L+1), not reach(198, L+1).
reach(181,L + 1) :- reach(198, L), L <= 13, L>= 10, not reach(197, L+1), not reach(199, L+1), not reach(215, L+1).
reach(197,L + 1) :- reach(198, L), L <= 13, L>= 10, not reach(181, L+1), not reach(199, L+1), not reach(215, L+1).
reach(199,L + 1) :- reach(198, L), L <= 13, L>= 10, not reach(181, L+1), not reach(197, L+1), not reach(215, L+1).
reach(215,L + 1) :- reach(198, L), L <= 13, L>= 10, not reach(181, L+1), not reach(197, L+1), not reach(199, L+1).
reach(182,L + 1) :- reach(199, L), L <= 14, L>= 11, not reach(198, L+1), not reach(200, L+1), not reach(216, L+1).
reach(198,L + 1) :- reach(199, L), L <= 14, L>= 11, not reach(182, L+1), not reach(200, L+1), not reach(216, L+1).
reach(200,L + 1) :- reach(199, L), L <= 14, L>= 11, not reach(182, L+1), not reach(198, L+1), not reach(216, L+1).
reach(216,L + 1) :- reach(199, L), L <= 14, L>= 11, not reach(182, L+1), not reach(198, L+1), not reach(200, L+1).
reach(183,L + 1) :- reach(200, L), L <= 15, L>= 12, not reach(199, L+1), not reach(201, L+1), not reach(217, L+1).
reach(199,L + 1) :- reach(200, L), L <= 15, L>= 12, not reach(183, L+1), not reach(201, L+1), not reach(217, L+1).
reach(201,L + 1) :- reach(200, L), L <= 15, L>= 12, not reach(183, L+1), not reach(199, L+1), not reach(217, L+1).
reach(217,L + 1) :- reach(200, L), L <= 15, L>= 12, not reach(183, L+1), not reach(199, L+1), not reach(201, L+1).
reach(184,L + 1) :- reach(201, L), L <= 14, L>= 13, not reach(200, L+1), not reach(202, L+1), not reach(218, L+1).
reach(200,L + 1) :- reach(201, L), L <= 14, L>= 13, not reach(184, L+1), not reach(202, L+1), not reach(218, L+1).
reach(202,L + 1) :- reach(201, L), L <= 14, L>= 13, not reach(184, L+1), not reach(200, L+1), not reach(218, L+1).
reach(218,L + 1) :- reach(201, L), L <= 14, L>= 13, not reach(184, L+1), not reach(200, L+1), not reach(202, L+1).
reach(185,L + 1) :- reach(202, L), L <= 14, L>= 13, not reach(201, L+1), not reach(203, L+1), not reach(219, L+1).
reach(201,L + 1) :- reach(202, L), L <= 14, L>= 13, not reach(185, L+1), not reach(203, L+1), not reach(219, L+1).
reach(203,L + 1) :- reach(202, L), L <= 14, L>= 13, not reach(185, L+1), not reach(201, L+1), not reach(219, L+1).
reach(219,L + 1) :- reach(202, L), L <= 14, L>= 13, not reach(185, L+1), not reach(201, L+1), not reach(203, L+1).
reach(186,L + 1) :- reach(203, L), L <= 13, L>= 12, not reach(202, L+1), not reach(220, L+1), not reach(221, L+1).
reach(202,L + 1) :- reach(203, L), L <= 13, L>= 12, not reach(186, L+1), not reach(220, L+1), not reach(221, L+1).
reach(220,L + 1) :- reach(203, L), L <= 13, L>= 12, not reach(186, L+1), not reach(202, L+1), not reach(221, L+1).
reach(221,L + 1) :- reach(203, L), L <= 13, L>= 12, not reach(186, L+1), not reach(202, L+1), not reach(220, L+1).
reach(188,L + 1) :- reach(205, L), L <= 20, L>= 15, not reach(206, L+1), not reach(222, L+1).
reach(206,L + 1) :- reach(205, L), L <= 20, L>= 15, not reach(188, L+1), not reach(222, L+1).
reach(222,L + 1) :- reach(205, L), L <= 20, L>= 15, not reach(188, L+1), not reach(206, L+1).
reach(188,L + 1) :- reach(270, L), L <= 18, L>= 13, not reach(253, L+1), not reach(269, L+1), not reach(271, L+1), not reach(287, L+1).
reach(253,L + 1) :- reach(270, L), L <= 18, L>= 13, not reach(188, L+1), not reach(269, L+1), not reach(271, L+1), not reach(287, L+1).
reach(269,L + 1) :- reach(270, L), L <= 18, L>= 13, not reach(188, L+1), not reach(253, L+1), not reach(271, L+1), not reach(287, L+1).
reach(271,L + 1) :- reach(270, L), L <= 18, L>= 13, not reach(188, L+1), not reach(253, L+1), not reach(269, L+1), not reach(287, L+1).
reach(287,L + 1) :- reach(270, L), L <= 18, L>= 13, not reach(188, L+1), not reach(253, L+1), not reach(269, L+1), not reach(271, L+1).
reach(189,L + 1) :- reach(206, L), L <= 19, L>= 15, not reach(205, L+1), not reach(207, L+1), not reach(223, L+1).
reach(205,L + 1) :- reach(206, L), L <= 19, L>= 15, not reach(189, L+1), not reach(207, L+1), not reach(223, L+1).
reach(207,L + 1) :- reach(206, L), L <= 19, L>= 15, not reach(189, L+1), not reach(205, L+1), not reach(223, L+1).
reach(223,L + 1) :- reach(206, L), L <= 19, L>= 15, not reach(189, L+1), not reach(205, L+1), not reach(207, L+1).
reach(190,L + 1) :- reach(207, L), L <= 18, L>= 14, not reach(206, L+1), not reach(208, L+1), not reach(224, L+1).
reach(206,L + 1) :- reach(207, L), L <= 18, L>= 14, not reach(190, L+1), not reach(208, L+1), not reach(224, L+1).
reach(208,L + 1) :- reach(207, L), L <= 18, L>= 14, not reach(190, L+1), not reach(206, L+1), not reach(224, L+1).
reach(224,L + 1) :- reach(207, L), L <= 18, L>= 14, not reach(190, L+1), not reach(206, L+1), not reach(208, L+1).
reach(191,L + 1) :- reach(208, L), L <= 17, L>= 13, not reach(207, L+1), not reach(209, L+1), not reach(225, L+1).
reach(207,L + 1) :- reach(208, L), L <= 17, L>= 13, not reach(191, L+1), not reach(209, L+1), not reach(225, L+1).
reach(209,L + 1) :- reach(208, L), L <= 17, L>= 13, not reach(191, L+1), not reach(207, L+1), not reach(225, L+1).
reach(225,L + 1) :- reach(208, L), L <= 17, L>= 13, not reach(191, L+1), not reach(207, L+1), not reach(209, L+1).
reach(192,L + 1) :- reach(209, L), L <= 16, L>= 14, not reach(208, L+1), not reach(210, L+1).
reach(208,L + 1) :- reach(209, L), L <= 16, L>= 14, not reach(192, L+1), not reach(210, L+1).
reach(210,L + 1) :- reach(209, L), L <= 16, L>= 14, not reach(192, L+1), not reach(208, L+1).
reach(193,L + 1) :- reach(210, L), L <= 15, L>= 14, not reach(209, L+1), not reach(211, L+1).
reach(209,L + 1) :- reach(210, L), L <= 15, L>= 14, not reach(193, L+1), not reach(211, L+1).
reach(211,L + 1) :- reach(210, L), L <= 15, L>= 14, not reach(193, L+1), not reach(209, L+1).
reach(194,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(210, L+1), not reach(212, L+1), not reach(228, L+1).
reach(210,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(194, L+1), not reach(212, L+1), not reach(228, L+1).
reach(212,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(194, L+1), not reach(210, L+1), not reach(228, L+1).
reach(228,L + 1) :- reach(211, L), L <= 15, L>= 13, not reach(194, L+1), not reach(210, L+1), not reach(212, L+1).
reach(195,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(211, L+1), not reach(213, L+1).
reach(211,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(195, L+1), not reach(213, L+1).
reach(213,L + 1) :- reach(212, L), L <= 14, L>= 12, not reach(195, L+1), not reach(211, L+1).
reach(196,L + 1) :- reach(213, L), L <= 15, L>= 11, not reach(212, L+1), not reach(214, L+1), not reach(230, L+1).
reach(212,L + 1) :- reach(213, L), L <= 15, L>= 11, not reach(196, L+1), not reach(214, L+1), not reach(230, L+1).
reach(214,L + 1) :- reach(213, L), L <= 15, L>= 11, not reach(196, L+1), not reach(212, L+1), not reach(230, L+1).
reach(230,L + 1) :- reach(213, L), L <= 15, L>= 11, not reach(196, L+1), not reach(212, L+1), not reach(214, L+1).
reach(197,L + 1) :- reach(214, L), L <= 14, L>= 10, not reach(213, L+1), not reach(215, L+1), not reach(231, L+1).
reach(213,L + 1) :- reach(214, L), L <= 14, L>= 10, not reach(197, L+1), not reach(215, L+1), not reach(231, L+1).
reach(215,L + 1) :- reach(214, L), L <= 14, L>= 10, not reach(197, L+1), not reach(213, L+1), not reach(231, L+1).
reach(231,L + 1) :- reach(214, L), L <= 14, L>= 10, not reach(197, L+1), not reach(213, L+1), not reach(215, L+1).
reach(198,L + 1) :- reach(215, L), L <= 14, L>= 11, not reach(214, L+1), not reach(216, L+1).
reach(214,L + 1) :- reach(215, L), L <= 14, L>= 11, not reach(198, L+1), not reach(216, L+1).
reach(216,L + 1) :- reach(215, L), L <= 14, L>= 11, not reach(198, L+1), not reach(214, L+1).
reach(199,L + 1) :- reach(216, L), L <= 15, L>= 12, not reach(215, L+1), not reach(217, L+1), not reach(232, L+1).
reach(215,L + 1) :- reach(216, L), L <= 15, L>= 12, not reach(199, L+1), not reach(217, L+1), not reach(232, L+1).
reach(217,L + 1) :- reach(216, L), L <= 15, L>= 12, not reach(199, L+1), not reach(215, L+1), not reach(232, L+1).
reach(232,L + 1) :- reach(216, L), L <= 15, L>= 12, not reach(199, L+1), not reach(215, L+1), not reach(217, L+1).
reach(200,L + 1) :- reach(217, L), L <= 16, L>= 13, not reach(216, L+1), not reach(218, L+1), not reach(233, L+1), not reach(219, L+1), not reach(263, L+1).
reach(216,L + 1) :- reach(217, L), L <= 16, L>= 13, not reach(200, L+1), not reach(218, L+1), not reach(233, L+1), not reach(219, L+1), not reach(263, L+1).
reach(218,L + 1) :- reach(217, L), L <= 16, L>= 13, not reach(200, L+1), not reach(216, L+1), not reach(233, L+1), not reach(219, L+1), not reach(263, L+1).
reach(233,L + 1) :- reach(217, L), L <= 16, L>= 13, not reach(200, L+1), not reach(216, L+1), not reach(218, L+1), not reach(219, L+1), not reach(263, L+1).
reach(219,L + 1) :- reach(217, L), L <= 16, L>= 13, not reach(200, L+1), not reach(216, L+1), not reach(218, L+1), not reach(233, L+1), not reach(263, L+1).
reach(263,L + 1) :- reach(217, L), L <= 16, L>= 13, not reach(200, L+1), not reach(216, L+1), not reach(218, L+1), not reach(233, L+1), not reach(219, L+1).
reach(201,L + 1) :- reach(218, L), L <= 15, L>= 13, not reach(217, L+1), not reach(219, L+1), not reach(234, L+1).
reach(217,L + 1) :- reach(218, L), L <= 15, L>= 13, not reach(201, L+1), not reach(219, L+1), not reach(234, L+1).
reach(219,L + 1) :- reach(218, L), L <= 15, L>= 13, not reach(201, L+1), not reach(217, L+1), not reach(234, L+1).
reach(234,L + 1) :- reach(218, L), L <= 15, L>= 13, not reach(201, L+1), not reach(217, L+1), not reach(219, L+1).
reach(202,L + 1) :- reach(219, L), L <= 15, L>= 12, not reach(217, L+1), not reach(218, L+1), not reach(220, L+1), not reach(235, L+1).
reach(217,L + 1) :- reach(219, L), L <= 15, L>= 12, not reach(202, L+1), not reach(218, L+1), not reach(220, L+1), not reach(235, L+1).
reach(218,L + 1) :- reach(219, L), L <= 15, L>= 12, not reach(202, L+1), not reach(217, L+1), not reach(220, L+1), not reach(235, L+1).
reach(220,L + 1) :- reach(219, L), L <= 15, L>= 12, not reach(202, L+1), not reach(217, L+1), not reach(218, L+1), not reach(235, L+1).
reach(235,L + 1) :- reach(219, L), L <= 15, L>= 12, not reach(202, L+1), not reach(217, L+1), not reach(218, L+1), not reach(220, L+1).
reach(203,L + 1) :- reach(220, L), L <= 14, L>= 11, not reach(219, L+1), not reach(221, L+1), not reach(236, L+1).
reach(219,L + 1) :- reach(220, L), L <= 14, L>= 11, not reach(203, L+1), not reach(221, L+1), not reach(236, L+1).
reach(221,L + 1) :- reach(220, L), L <= 14, L>= 11, not reach(203, L+1), not reach(219, L+1), not reach(236, L+1).
reach(236,L + 1) :- reach(220, L), L <= 14, L>= 11, not reach(203, L+1), not reach(219, L+1), not reach(221, L+1).
reach(220,L + 1) :- reach(221, L), L <= 14, L>= 12, not reach(237, L+1), not reach(203, L+1).
reach(237,L + 1) :- reach(221, L), L <= 14, L>= 12, not reach(220, L+1), not reach(203, L+1).
reach(203,L + 1) :- reach(221, L), L <= 14, L>= 12, not reach(220, L+1), not reach(237, L+1).
reach(205,L + 1) :- reach(222, L), L <= 21, L>= 16, not reach(223, L+1), not reach(238, L+1).
reach(223,L + 1) :- reach(222, L), L <= 21, L>= 16, not reach(205, L+1), not reach(238, L+1).
reach(238,L + 1) :- reach(222, L), L <= 21, L>= 16, not reach(205, L+1), not reach(223, L+1).
reach(206,L + 1) :- reach(223, L), L <= 20, L>= 16, not reach(222, L+1), not reach(224, L+1), not reach(239, L+1).
reach(222,L + 1) :- reach(223, L), L <= 20, L>= 16, not reach(206, L+1), not reach(224, L+1), not reach(239, L+1).
reach(224,L + 1) :- reach(223, L), L <= 20, L>= 16, not reach(206, L+1), not reach(222, L+1), not reach(239, L+1).
reach(239,L + 1) :- reach(223, L), L <= 20, L>= 16, not reach(206, L+1), not reach(222, L+1), not reach(224, L+1).
reach(207,L + 1) :- reach(224, L), L <= 19, L>= 15, not reach(223, L+1), not reach(225, L+1), not reach(240, L+1), not reach(245, L+1).
reach(223,L + 1) :- reach(224, L), L <= 19, L>= 15, not reach(207, L+1), not reach(225, L+1), not reach(240, L+1), not reach(245, L+1).
reach(225,L + 1) :- reach(224, L), L <= 19, L>= 15, not reach(207, L+1), not reach(223, L+1), not reach(240, L+1), not reach(245, L+1).
reach(240,L + 1) :- reach(224, L), L <= 19, L>= 15, not reach(207, L+1), not reach(223, L+1), not reach(225, L+1), not reach(245, L+1).
reach(245,L + 1) :- reach(224, L), L <= 19, L>= 15, not reach(207, L+1), not reach(223, L+1), not reach(225, L+1), not reach(240, L+1).
reach(208,L + 1) :- reach(225, L), L <= 18, L>= 14, not reach(224, L+1), not reach(226, L+1), not reach(241, L+1).
reach(224,L + 1) :- reach(225, L), L <= 18, L>= 14, not reach(208, L+1), not reach(226, L+1), not reach(241, L+1).
reach(226,L + 1) :- reach(225, L), L <= 18, L>= 14, not reach(208, L+1), not reach(224, L+1), not reach(241, L+1).
reach(241,L + 1) :- reach(225, L), L <= 18, L>= 14, not reach(208, L+1), not reach(224, L+1), not reach(226, L+1).
reach(211,L + 1) :- reach(228, L), L <= 16, L>= 14, not reach(227, L+1), not reach(229, L+1), not reach(244, L+1).
reach(227,L + 1) :- reach(228, L), L <= 16, L>= 14, not reach(211, L+1), not reach(229, L+1), not reach(244, L+1).
reach(229,L + 1) :- reach(228, L), L <= 16, L>= 14, not reach(211, L+1), not reach(227, L+1), not reach(244, L+1).
reach(244,L + 1) :- reach(228, L), L <= 16, L>= 14, not reach(211, L+1), not reach(227, L+1), not reach(229, L+1).
reach(213,L + 1) :- reach(230, L), L <= 16, L>= 12, not reach(229, L+1), not reach(231, L+1), not reach(246, L+1).
reach(229,L + 1) :- reach(230, L), L <= 16, L>= 12, not reach(213, L+1), not reach(231, L+1), not reach(246, L+1).
reach(231,L + 1) :- reach(230, L), L <= 16, L>= 12, not reach(213, L+1), not reach(229, L+1), not reach(246, L+1).
reach(246,L + 1) :- reach(230, L), L <= 16, L>= 12, not reach(213, L+1), not reach(229, L+1), not reach(231, L+1).
reach(214,L + 1) :- reach(231, L), L <= 15, L>= 11, not reach(230, L+1), not reach(247, L+1), not reach(248, L+1).
reach(230,L + 1) :- reach(231, L), L <= 15, L>= 11, not reach(214, L+1), not reach(247, L+1), not reach(248, L+1).
reach(247,L + 1) :- reach(231, L), L <= 15, L>= 11, not reach(214, L+1), not reach(230, L+1), not reach(248, L+1).
reach(248,L + 1) :- reach(231, L), L <= 15, L>= 11, not reach(214, L+1), not reach(230, L+1), not reach(247, L+1).
reach(216,L + 1) :- reach(232, L), L <= 14, L>= 13, not reach(247, L+1), not reach(233, L+1), not reach(250, L+1).
reach(247,L + 1) :- reach(232, L), L <= 14, L>= 13, not reach(216, L+1), not reach(233, L+1), not reach(250, L+1).
reach(233,L + 1) :- reach(232, L), L <= 14, L>= 13, not reach(216, L+1), not reach(247, L+1), not reach(250, L+1).
reach(250,L + 1) :- reach(232, L), L <= 14, L>= 13, not reach(216, L+1), not reach(247, L+1), not reach(233, L+1).
reach(217,L + 1) :- reach(233, L), L <= 15, L>= 12, not reach(232, L+1), not reach(234, L+1), not reach(251, L+1).
reach(232,L + 1) :- reach(233, L), L <= 15, L>= 12, not reach(217, L+1), not reach(234, L+1), not reach(251, L+1).
reach(234,L + 1) :- reach(233, L), L <= 15, L>= 12, not reach(217, L+1), not reach(232, L+1), not reach(251, L+1).
reach(251,L + 1) :- reach(233, L), L <= 15, L>= 12, not reach(217, L+1), not reach(232, L+1), not reach(234, L+1).
reach(217,L + 1) :- reach(263, L), L <= 17, L>= 14, not reach(245, L+1), not reach(262, L+1), not reach(264, L+1), not reach(280, L+1).
reach(245,L + 1) :- reach(263, L), L <= 17, L>= 14, not reach(217, L+1), not reach(262, L+1), not reach(264, L+1), not reach(280, L+1).
reach(262,L + 1) :- reach(263, L), L <= 17, L>= 14, not reach(217, L+1), not reach(245, L+1), not reach(264, L+1), not reach(280, L+1).
reach(264,L + 1) :- reach(263, L), L <= 17, L>= 14, not reach(217, L+1), not reach(245, L+1), not reach(262, L+1), not reach(280, L+1).
reach(280,L + 1) :- reach(263, L), L <= 17, L>= 14, not reach(217, L+1), not reach(245, L+1), not reach(262, L+1), not reach(264, L+1).
reach(218,L + 1) :- reach(234, L), L <= 15, L>= 12, not reach(233, L+1), not reach(235, L+1), not reach(252, L+1).
reach(233,L + 1) :- reach(234, L), L <= 15, L>= 12, not reach(218, L+1), not reach(235, L+1), not reach(252, L+1).
reach(235,L + 1) :- reach(234, L), L <= 15, L>= 12, not reach(218, L+1), not reach(233, L+1), not reach(252, L+1).
reach(252,L + 1) :- reach(234, L), L <= 15, L>= 12, not reach(218, L+1), not reach(233, L+1), not reach(235, L+1).
reach(219,L + 1) :- reach(235, L), L <= 16, L>= 11, not reach(234, L+1), not reach(236, L+1), not reach(253, L+1).
reach(234,L + 1) :- reach(235, L), L <= 16, L>= 11, not reach(219, L+1), not reach(236, L+1), not reach(253, L+1).
reach(236,L + 1) :- reach(235, L), L <= 16, L>= 11, not reach(219, L+1), not reach(234, L+1), not reach(253, L+1).
reach(253,L + 1) :- reach(235, L), L <= 16, L>= 11, not reach(219, L+1), not reach(234, L+1), not reach(236, L+1).
reach(221,L + 1) :- reach(237, L), L <= 15, L>= 11, not reach(236, L+1), not reach(255, L+1).
reach(236,L + 1) :- reach(237, L), L <= 15, L>= 11, not reach(221, L+1), not reach(255, L+1).
reach(255,L + 1) :- reach(237, L), L <= 15, L>= 11, not reach(221, L+1), not reach(236, L+1).
reach(222,L + 1) :- reach(238, L), L <= 22, L>= 17, not reach(239, L+1), not reach(256, L+1).
reach(239,L + 1) :- reach(238, L), L <= 22, L>= 17, not reach(222, L+1), not reach(256, L+1).
reach(256,L + 1) :- reach(238, L), L <= 22, L>= 17, not reach(222, L+1), not reach(239, L+1).
reach(223,L + 1) :- reach(239, L), L <= 21, L>= 17, not reach(238, L+1), not reach(240, L+1), not reach(257, L+1).
reach(238,L + 1) :- reach(239, L), L <= 21, L>= 17, not reach(223, L+1), not reach(240, L+1), not reach(257, L+1).
reach(240,L + 1) :- reach(239, L), L <= 21, L>= 17, not reach(223, L+1), not reach(238, L+1), not reach(257, L+1).
reach(257,L + 1) :- reach(239, L), L <= 21, L>= 17, not reach(223, L+1), not reach(238, L+1), not reach(240, L+1).
reach(224,L + 1) :- reach(240, L), L <= 20, L>= 16, not reach(239, L+1), not reach(241, L+1), not reach(258, L+1).
reach(239,L + 1) :- reach(240, L), L <= 20, L>= 16, not reach(224, L+1), not reach(241, L+1), not reach(258, L+1).
reach(241,L + 1) :- reach(240, L), L <= 20, L>= 16, not reach(224, L+1), not reach(239, L+1), not reach(258, L+1).
reach(258,L + 1) :- reach(240, L), L <= 20, L>= 16, not reach(224, L+1), not reach(239, L+1), not reach(241, L+1).
reach(224,L + 1) :- reach(245, L), L <= 18, L>= 14, not reach(229, L+1), not reach(244, L+1), not reach(246, L+1), not reach(263, L+1).
reach(229,L + 1) :- reach(245, L), L <= 18, L>= 14, not reach(224, L+1), not reach(244, L+1), not reach(246, L+1), not reach(263, L+1).
reach(244,L + 1) :- reach(245, L), L <= 18, L>= 14, not reach(224, L+1), not reach(229, L+1), not reach(246, L+1), not reach(263, L+1).
reach(246,L + 1) :- reach(245, L), L <= 18, L>= 14, not reach(224, L+1), not reach(229, L+1), not reach(244, L+1), not reach(263, L+1).
reach(263,L + 1) :- reach(245, L), L <= 18, L>= 14, not reach(224, L+1), not reach(229, L+1), not reach(244, L+1), not reach(246, L+1).
reach(225,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(227, L+1), not reach(242, L+1).
reach(227,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(225, L+1), not reach(242, L+1).
reach(242,L + 1) :- reach(226, L), L <= 17, L>= 15, not reach(225, L+1), not reach(227, L+1).
reach(225,L + 1) :- reach(241, L), L <= 19, L>= 15, not reach(240, L+1), not reach(242, L+1), not reach(259, L+1).
reach(240,L + 1) :- reach(241, L), L <= 19, L>= 15, not reach(225, L+1), not reach(242, L+1), not reach(259, L+1).
reach(242,L + 1) :- reach(241, L), L <= 19, L>= 15, not reach(225, L+1), not reach(240, L+1), not reach(259, L+1).
reach(259,L + 1) :- reach(241, L), L <= 19, L>= 15, not reach(225, L+1), not reach(240, L+1), not reach(242, L+1).
reach(226,L + 1) :- reach(227, L), L <= 16, L>= 15, not reach(228, L+1), not reach(243, L+1).
reach(228,L + 1) :- reach(227, L), L <= 16, L>= 15, not reach(226, L+1), not reach(243, L+1).
reach(243,L + 1) :- reach(227, L), L <= 16, L>= 15, not reach(226, L+1), not reach(228, L+1).
reach(226,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(241, L+1), not reach(243, L+1), not reach(260, L+1).
reach(241,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(226, L+1), not reach(243, L+1), not reach(260, L+1).
reach(243,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(226, L+1), not reach(241, L+1), not reach(260, L+1).
reach(260,L + 1) :- reach(242, L), L <= 18, L>= 16, not reach(226, L+1), not reach(241, L+1), not reach(243, L+1).
reach(227,L + 1) :- reach(243, L), L <= 17, L>= 16, not reach(242, L+1), not reach(244, L+1), not reach(261, L+1).
reach(242,L + 1) :- reach(243, L), L <= 17, L>= 16, not reach(227, L+1), not reach(244, L+1), not reach(261, L+1).
reach(244,L + 1) :- reach(243, L), L <= 17, L>= 16, not reach(227, L+1), not reach(242, L+1), not reach(261, L+1).
reach(261,L + 1) :- reach(243, L), L <= 17, L>= 16, not reach(227, L+1), not reach(242, L+1), not reach(244, L+1).
reach(228,L + 1) :- reach(229, L), L <= 17, L>= 13, not reach(230, L+1), not reach(245, L+1).
reach(230,L + 1) :- reach(229, L), L <= 17, L>= 13, not reach(228, L+1), not reach(245, L+1).
reach(245,L + 1) :- reach(229, L), L <= 17, L>= 13, not reach(228, L+1), not reach(230, L+1).
reach(228,L + 1) :- reach(244, L), L <= 17, L>= 15, not reach(243, L+1), not reach(245, L+1), not reach(262, L+1).
reach(243,L + 1) :- reach(244, L), L <= 17, L>= 15, not reach(228, L+1), not reach(245, L+1), not reach(262, L+1).
reach(245,L + 1) :- reach(244, L), L <= 17, L>= 15, not reach(228, L+1), not reach(243, L+1), not reach(262, L+1).
reach(262,L + 1) :- reach(244, L), L <= 17, L>= 15, not reach(228, L+1), not reach(243, L+1), not reach(245, L+1).
reach(230,L + 1) :- reach(246, L), L <= 17, L>= 13, not reach(245, L+1), not reach(248, L+1), not reach(264, L+1).
reach(245,L + 1) :- reach(246, L), L <= 17, L>= 13, not reach(230, L+1), not reach(248, L+1), not reach(264, L+1).
reach(248,L + 1) :- reach(246, L), L <= 17, L>= 13, not reach(230, L+1), not reach(245, L+1), not reach(264, L+1).
reach(264,L + 1) :- reach(246, L), L <= 17, L>= 13, not reach(230, L+1), not reach(245, L+1), not reach(248, L+1).
reach(231,L + 1) :- reach(247, L), L <= 14, L>= 12, not reach(232, L+1), not reach(249, L+1).
reach(232,L + 1) :- reach(247, L), L <= 14, L>= 12, not reach(231, L+1), not reach(249, L+1).
reach(249,L + 1) :- reach(247, L), L <= 14, L>= 12, not reach(231, L+1), not reach(232, L+1).
reach(231,L + 1) :- reach(248, L), L <= 16, L>= 12, not reach(246, L+1), not reach(249, L+1), not reach(265, L+1).
reach(246,L + 1) :- reach(248, L), L <= 16, L>= 12, not reach(231, L+1), not reach(249, L+1), not reach(265, L+1).
reach(249,L + 1) :- reach(248, L), L <= 16, L>= 12, not reach(231, L+1), not reach(246, L+1), not reach(265, L+1).
reach(265,L + 1) :- reach(248, L), L <= 16, L>= 12, not reach(231, L+1), not reach(246, L+1), not reach(249, L+1).
reach(232,L + 1) :- reach(250, L), L <= 14, L>= 14, not reach(249, L+1), not reach(267, L+1).
reach(249,L + 1) :- reach(250, L), L <= 14, L>= 14, not reach(232, L+1), not reach(267, L+1).
reach(267,L + 1) :- reach(250, L), L <= 14, L>= 14, not reach(232, L+1), not reach(249, L+1).
reach(234,L + 1) :- reach(252, L), L <= 16, L>= 12, not reach(251, L+1), not reach(253, L+1), not reach(269, L+1).
reach(251,L + 1) :- reach(252, L), L <= 16, L>= 12, not reach(234, L+1), not reach(253, L+1), not reach(269, L+1).
reach(253,L + 1) :- reach(252, L), L <= 16, L>= 12, not reach(234, L+1), not reach(251, L+1), not reach(269, L+1).
reach(269,L + 1) :- reach(252, L), L <= 16, L>= 12, not reach(234, L+1), not reach(251, L+1), not reach(253, L+1).
reach(235,L + 1) :- reach(253, L), L <= 17, L>= 12, not reach(252, L+1), not reach(254, L+1), not reach(270, L+1).
reach(252,L + 1) :- reach(253, L), L <= 17, L>= 12, not reach(235, L+1), not reach(254, L+1), not reach(270, L+1).
reach(254,L + 1) :- reach(253, L), L <= 17, L>= 12, not reach(235, L+1), not reach(252, L+1), not reach(270, L+1).
reach(270,L + 1) :- reach(253, L), L <= 17, L>= 12, not reach(235, L+1), not reach(252, L+1), not reach(254, L+1).
reach(236,L + 1) :- reach(254, L), L <= 16, L>= 11, not reach(253, L+1), not reach(255, L+1), not reach(271, L+1).
reach(253,L + 1) :- reach(254, L), L <= 16, L>= 11, not reach(236, L+1), not reach(255, L+1), not reach(271, L+1).
reach(255,L + 1) :- reach(254, L), L <= 16, L>= 11, not reach(236, L+1), not reach(253, L+1), not reach(271, L+1).
reach(271,L + 1) :- reach(254, L), L <= 16, L>= 11, not reach(236, L+1), not reach(253, L+1), not reach(255, L+1).
reach(237,L + 1) :- reach(255, L), L <= 16, L>= 12, not reach(254, L+1), not reach(271, L+1).
reach(254,L + 1) :- reach(255, L), L <= 16, L>= 12, not reach(237, L+1), not reach(271, L+1).
reach(271,L + 1) :- reach(255, L), L <= 16, L>= 12, not reach(237, L+1), not reach(254, L+1).
reach(238,L + 1) :- reach(256, L), L <= 23, L>= 18, not reach(257, L+1), not reach(273, L+1).
reach(257,L + 1) :- reach(256, L), L <= 23, L>= 18, not reach(238, L+1), not reach(273, L+1).
reach(273,L + 1) :- reach(256, L), L <= 23, L>= 18, not reach(238, L+1), not reach(257, L+1).
reach(239,L + 1) :- reach(257, L), L <= 22, L>= 18, not reach(256, L+1), not reach(258, L+1), not reach(274, L+1).
reach(256,L + 1) :- reach(257, L), L <= 22, L>= 18, not reach(239, L+1), not reach(258, L+1), not reach(274, L+1).
reach(258,L + 1) :- reach(257, L), L <= 22, L>= 18, not reach(239, L+1), not reach(256, L+1), not reach(274, L+1).
reach(274,L + 1) :- reach(257, L), L <= 22, L>= 18, not reach(239, L+1), not reach(256, L+1), not reach(258, L+1).
reach(240,L + 1) :- reach(258, L), L <= 21, L>= 17, not reach(257, L+1), not reach(259, L+1), not reach(275, L+1).
reach(257,L + 1) :- reach(258, L), L <= 21, L>= 17, not reach(240, L+1), not reach(259, L+1), not reach(275, L+1).
reach(259,L + 1) :- reach(258, L), L <= 21, L>= 17, not reach(240, L+1), not reach(257, L+1), not reach(275, L+1).
reach(275,L + 1) :- reach(258, L), L <= 21, L>= 17, not reach(240, L+1), not reach(257, L+1), not reach(259, L+1).
reach(241,L + 1) :- reach(259, L), L <= 20, L>= 16, not reach(258, L+1), not reach(260, L+1), not reach(276, L+1).
reach(258,L + 1) :- reach(259, L), L <= 20, L>= 16, not reach(241, L+1), not reach(260, L+1), not reach(276, L+1).
reach(260,L + 1) :- reach(259, L), L <= 20, L>= 16, not reach(241, L+1), not reach(258, L+1), not reach(276, L+1).
reach(276,L + 1) :- reach(259, L), L <= 20, L>= 16, not reach(241, L+1), not reach(258, L+1), not reach(260, L+1).
reach(242,L + 1) :- reach(260, L), L <= 19, L>= 17, not reach(259, L+1), not reach(261, L+1), not reach(277, L+1).
reach(259,L + 1) :- reach(260, L), L <= 19, L>= 17, not reach(242, L+1), not reach(261, L+1), not reach(277, L+1).
reach(261,L + 1) :- reach(260, L), L <= 19, L>= 17, not reach(242, L+1), not reach(259, L+1), not reach(277, L+1).
reach(277,L + 1) :- reach(260, L), L <= 19, L>= 17, not reach(242, L+1), not reach(259, L+1), not reach(261, L+1).
reach(243,L + 1) :- reach(261, L), L <= 18, L>= 16, not reach(260, L+1), not reach(262, L+1), not reach(278, L+1).
reach(260,L + 1) :- reach(261, L), L <= 18, L>= 16, not reach(243, L+1), not reach(262, L+1), not reach(278, L+1).
reach(262,L + 1) :- reach(261, L), L <= 18, L>= 16, not reach(243, L+1), not reach(260, L+1), not reach(278, L+1).
reach(278,L + 1) :- reach(261, L), L <= 18, L>= 16, not reach(243, L+1), not reach(260, L+1), not reach(262, L+1).
reach(244,L + 1) :- reach(262, L), L <= 17, L>= 15, not reach(261, L+1), not reach(263, L+1), not reach(279, L+1).
reach(261,L + 1) :- reach(262, L), L <= 17, L>= 15, not reach(244, L+1), not reach(263, L+1), not reach(279, L+1).
reach(263,L + 1) :- reach(262, L), L <= 17, L>= 15, not reach(244, L+1), not reach(261, L+1), not reach(279, L+1).
reach(279,L + 1) :- reach(262, L), L <= 17, L>= 15, not reach(244, L+1), not reach(261, L+1), not reach(263, L+1).
reach(246,L + 1) :- reach(264, L), L <= 16, L>= 14, not reach(263, L+1), not reach(265, L+1), not reach(281, L+1).
reach(263,L + 1) :- reach(264, L), L <= 16, L>= 14, not reach(246, L+1), not reach(265, L+1), not reach(281, L+1).
reach(265,L + 1) :- reach(264, L), L <= 16, L>= 14, not reach(246, L+1), not reach(263, L+1), not reach(281, L+1).
reach(281,L + 1) :- reach(264, L), L <= 16, L>= 14, not reach(246, L+1), not reach(263, L+1), not reach(265, L+1).
reach(247,L + 1) :- reach(249, L), L <= 15, L>= 13, not reach(248, L+1), not reach(250, L+1), not reach(266, L+1).
reach(248,L + 1) :- reach(249, L), L <= 15, L>= 13, not reach(247, L+1), not reach(250, L+1), not reach(266, L+1).
reach(250,L + 1) :- reach(249, L), L <= 15, L>= 13, not reach(247, L+1), not reach(248, L+1), not reach(266, L+1).
reach(266,L + 1) :- reach(249, L), L <= 15, L>= 13, not reach(247, L+1), not reach(248, L+1), not reach(250, L+1).
reach(248,L + 1) :- reach(265, L), L <= 15, L>= 13, not reach(264, L+1), not reach(266, L+1), not reach(282, L+1).
reach(264,L + 1) :- reach(265, L), L <= 15, L>= 13, not reach(248, L+1), not reach(266, L+1), not reach(282, L+1).
reach(266,L + 1) :- reach(265, L), L <= 15, L>= 13, not reach(248, L+1), not reach(264, L+1), not reach(282, L+1).
reach(282,L + 1) :- reach(265, L), L <= 15, L>= 13, not reach(248, L+1), not reach(264, L+1), not reach(266, L+1).
reach(249,L + 1) :- reach(266, L), L <= 14, L>= 14, not reach(265, L+1), not reach(267, L+1).
reach(265,L + 1) :- reach(266, L), L <= 14, L>= 14, not reach(249, L+1), not reach(267, L+1).
reach(267,L + 1) :- reach(266, L), L <= 14, L>= 14, not reach(249, L+1), not reach(265, L+1).
reach(250,L + 1) :- reach(267, L), L <= 15, L>= 13, not reach(266, L+1), not reach(268, L+1), not reach(284, L+1).
reach(266,L + 1) :- reach(267, L), L <= 15, L>= 13, not reach(250, L+1), not reach(268, L+1), not reach(284, L+1).
reach(268,L + 1) :- reach(267, L), L <= 15, L>= 13, not reach(250, L+1), not reach(266, L+1), not reach(284, L+1).
reach(284,L + 1) :- reach(267, L), L <= 15, L>= 13, not reach(250, L+1), not reach(266, L+1), not reach(268, L+1).
reach(251,L + 1) :- reach(268, L), L <= 16, L>= 12, not reach(267, L+1), not reach(269, L+1), not reach(285, L+1).
reach(267,L + 1) :- reach(268, L), L <= 16, L>= 12, not reach(251, L+1), not reach(269, L+1), not reach(285, L+1).
reach(269,L + 1) :- reach(268, L), L <= 16, L>= 12, not reach(251, L+1), not reach(267, L+1), not reach(285, L+1).
reach(285,L + 1) :- reach(268, L), L <= 16, L>= 12, not reach(251, L+1), not reach(267, L+1), not reach(269, L+1).
reach(252,L + 1) :- reach(269, L), L <= 17, L>= 13, not reach(268, L+1), not reach(270, L+1), not reach(286, L+1).
reach(268,L + 1) :- reach(269, L), L <= 17, L>= 13, not reach(252, L+1), not reach(270, L+1), not reach(286, L+1).
reach(270,L + 1) :- reach(269, L), L <= 17, L>= 13, not reach(252, L+1), not reach(268, L+1), not reach(286, L+1).
reach(286,L + 1) :- reach(269, L), L <= 17, L>= 13, not reach(252, L+1), not reach(268, L+1), not reach(270, L+1).
reach(254,L + 1) :- reach(271, L), L <= 17, L>= 12, not reach(270, L+1), not reach(255, L+1), not reach(287, L+1).
reach(270,L + 1) :- reach(271, L), L <= 17, L>= 12, not reach(254, L+1), not reach(255, L+1), not reach(287, L+1).
reach(255,L + 1) :- reach(271, L), L <= 17, L>= 12, not reach(254, L+1), not reach(270, L+1), not reach(287, L+1).
reach(287,L + 1) :- reach(271, L), L <= 17, L>= 12, not reach(254, L+1), not reach(270, L+1), not reach(255, L+1).
reach(257,L + 1) :- reach(274, L), L <= 23, L>= 19, not reach(273, L+1), not reach(275, L+1).
reach(273,L + 1) :- reach(274, L), L <= 23, L>= 19, not reach(257, L+1), not reach(275, L+1).
reach(275,L + 1) :- reach(274, L), L <= 23, L>= 19, not reach(257, L+1), not reach(273, L+1).
reach(258,L + 1) :- reach(275, L), L <= 22, L>= 18, not reach(274, L+1), not reach(276, L+1).
reach(274,L + 1) :- reach(275, L), L <= 22, L>= 18, not reach(258, L+1), not reach(276, L+1).
reach(276,L + 1) :- reach(275, L), L <= 22, L>= 18, not reach(258, L+1), not reach(274, L+1).
reach(259,L + 1) :- reach(276, L), L <= 21, L>= 17, not reach(275, L+1), not reach(277, L+1).
reach(275,L + 1) :- reach(276, L), L <= 21, L>= 17, not reach(259, L+1), not reach(277, L+1).
reach(277,L + 1) :- reach(276, L), L <= 21, L>= 17, not reach(259, L+1), not reach(275, L+1).
reach(260,L + 1) :- reach(277, L), L <= 20, L>= 18, not reach(276, L+1), not reach(278, L+1).
reach(276,L + 1) :- reach(277, L), L <= 20, L>= 18, not reach(260, L+1), not reach(278, L+1).
reach(278,L + 1) :- reach(277, L), L <= 20, L>= 18, not reach(260, L+1), not reach(276, L+1).
reach(261,L + 1) :- reach(278, L), L <= 19, L>= 17, not reach(277, L+1), not reach(279, L+1).
reach(277,L + 1) :- reach(278, L), L <= 19, L>= 17, not reach(261, L+1), not reach(279, L+1).
reach(279,L + 1) :- reach(278, L), L <= 19, L>= 17, not reach(261, L+1), not reach(277, L+1).
reach(262,L + 1) :- reach(279, L), L <= 18, L>= 16, not reach(278, L+1), not reach(280, L+1).
reach(278,L + 1) :- reach(279, L), L <= 18, L>= 16, not reach(262, L+1), not reach(280, L+1).
reach(280,L + 1) :- reach(279, L), L <= 18, L>= 16, not reach(262, L+1), not reach(278, L+1).
reach(263,L + 1) :- reach(280, L), L <= 17, L>= 15, not reach(279, L+1), not reach(281, L+1).
reach(279,L + 1) :- reach(280, L), L <= 17, L>= 15, not reach(263, L+1), not reach(281, L+1).
reach(281,L + 1) :- reach(280, L), L <= 17, L>= 15, not reach(263, L+1), not reach(279, L+1).
reach(264,L + 1) :- reach(281, L), L <= 16, L>= 15, not reach(280, L+1), not reach(282, L+1).
reach(280,L + 1) :- reach(281, L), L <= 16, L>= 15, not reach(264, L+1), not reach(282, L+1).
reach(282,L + 1) :- reach(281, L), L <= 16, L>= 15, not reach(264, L+1), not reach(280, L+1).
reach(265,L + 1) :- reach(282, L), L <= 15, L>= 14, not reach(281, L+1), not reach(284, L+1).
reach(281,L + 1) :- reach(282, L), L <= 15, L>= 14, not reach(265, L+1), not reach(284, L+1).
reach(284,L + 1) :- reach(282, L), L <= 15, L>= 14, not reach(265, L+1), not reach(281, L+1).
reach(267,L + 1) :- reach(284, L), L <= 14, L>= 14, not reach(285, L+1), not reach(282, L+1).
reach(285,L + 1) :- reach(284, L), L <= 14, L>= 14, not reach(267, L+1), not reach(282, L+1).
reach(282,L + 1) :- reach(284, L), L <= 14, L>= 14, not reach(267, L+1), not reach(285, L+1).
reach(268,L + 1) :- reach(285, L), L <= 15, L>= 13, not reach(284, L+1), not reach(286, L+1).
reach(284,L + 1) :- reach(285, L), L <= 15, L>= 13, not reach(268, L+1), not reach(286, L+1).
reach(286,L + 1) :- reach(285, L), L <= 15, L>= 13, not reach(268, L+1), not reach(284, L+1).
reach(269,L + 1) :- reach(286, L), L <= 16, L>= 14, not reach(285, L+1), not reach(287, L+1).
reach(285,L + 1) :- reach(286, L), L <= 16, L>= 14, not reach(269, L+1), not reach(287, L+1).
reach(287,L + 1) :- reach(286, L), L <= 16, L>= 14, not reach(269, L+1), not reach(285, L+1).
reach(270,L + 1) :- reach(287, L), L <= 17, L>= 13, not reach(286, L+1), not reach(271, L+1).
reach(286,L + 1) :- reach(287, L), L <= 17, L>= 13, not reach(270, L+1), not reach(271, L+1).
reach(271,L + 1) :- reach(287, L), L <= 17, L>= 13, not reach(270, L+1), not reach(286, L+1).

