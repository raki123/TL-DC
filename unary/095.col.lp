reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13).
start(18).
goal(51).
:- reach(X, L), reach(X, L'), L != L'.
reach(45,L + 1) :- reach(1, L), L <= 10, L>= 3, not reach(88, L+1), not reach(68, L+1).
reach(88,L + 1) :- reach(1, L), L <= 10, L>= 3, not reach(45, L+1), not reach(68, L+1).
reach(68,L + 1) :- reach(1, L), L <= 10, L>= 3, not reach(45, L+1), not reach(88, L+1).
reach(1,L + 1) :- reach(45, L), L <= 9, L>= 2, not reach(35, L+1), not reach(84, L+1), not reach(28, L+1).
reach(35,L + 1) :- reach(45, L), L <= 9, L>= 2, not reach(1, L+1), not reach(84, L+1), not reach(28, L+1).
reach(84,L + 1) :- reach(45, L), L <= 9, L>= 2, not reach(1, L+1), not reach(35, L+1), not reach(28, L+1).
reach(28,L + 1) :- reach(45, L), L <= 9, L>= 2, not reach(1, L+1), not reach(35, L+1), not reach(84, L+1).
reach(1,L + 1) :- reach(88, L), L <= 9, L>= 4, not reach(10, L+1), not reach(34, L+1), not reach(92, L+1).
reach(10,L + 1) :- reach(88, L), L <= 9, L>= 4, not reach(1, L+1), not reach(34, L+1), not reach(92, L+1).
reach(34,L + 1) :- reach(88, L), L <= 9, L>= 4, not reach(1, L+1), not reach(10, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(88, L), L <= 9, L>= 4, not reach(1, L+1), not reach(10, L+1), not reach(34, L+1).
reach(10,L + 1) :- reach(3, L), L <= 9, L>= 3, not reach(57, L+1), not reach(69, L+1), not reach(38, L+1).
reach(57,L + 1) :- reach(3, L), L <= 9, L>= 3, not reach(10, L+1), not reach(69, L+1), not reach(38, L+1).
reach(69,L + 1) :- reach(3, L), L <= 9, L>= 3, not reach(10, L+1), not reach(57, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(3, L), L <= 9, L>= 3, not reach(10, L+1), not reach(57, L+1), not reach(69, L+1).
reach(3,L + 1) :- reach(10, L), L <= 8, L>= 4, not reach(88, L+1), not reach(42, L+1).
reach(88,L + 1) :- reach(10, L), L <= 8, L>= 4, not reach(3, L+1), not reach(42, L+1).
reach(42,L + 1) :- reach(10, L), L <= 8, L>= 4, not reach(3, L+1), not reach(88, L+1).
reach(3,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(30, L+1), not reach(92, L+1), not reach(61, L+1).
reach(20,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(3, L+1), not reach(21, L+1), not reach(22, L+1), not reach(30, L+1), not reach(92, L+1), not reach(61, L+1).
reach(21,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(3, L+1), not reach(20, L+1), not reach(22, L+1), not reach(30, L+1), not reach(92, L+1), not reach(61, L+1).
reach(22,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(3, L+1), not reach(20, L+1), not reach(21, L+1), not reach(30, L+1), not reach(92, L+1), not reach(61, L+1).
reach(30,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(3, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(92, L+1), not reach(61, L+1).
reach(92,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(3, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(30, L+1), not reach(61, L+1).
reach(61,L + 1) :- reach(57, L), L <= 8, L>= 2, not reach(3, L+1), not reach(20, L+1), not reach(21, L+1), not reach(22, L+1), not reach(30, L+1), not reach(92, L+1).
reach(3,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(36, L+1), not reach(40, L+1), not reach(68, L+1), not reach(79, L+1), not reach(92, L+1), not reach(22, L+1), not reach(16, L+1), not reach(85, L+1).
reach(36,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(40, L+1), not reach(68, L+1), not reach(79, L+1), not reach(92, L+1), not reach(22, L+1), not reach(16, L+1), not reach(85, L+1).
reach(40,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(68, L+1), not reach(79, L+1), not reach(92, L+1), not reach(22, L+1), not reach(16, L+1), not reach(85, L+1).
reach(68,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(40, L+1), not reach(79, L+1), not reach(92, L+1), not reach(22, L+1), not reach(16, L+1), not reach(85, L+1).
reach(79,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(40, L+1), not reach(68, L+1), not reach(92, L+1), not reach(22, L+1), not reach(16, L+1), not reach(85, L+1).
reach(92,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(40, L+1), not reach(68, L+1), not reach(79, L+1), not reach(22, L+1), not reach(16, L+1), not reach(85, L+1).
reach(22,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(40, L+1), not reach(68, L+1), not reach(79, L+1), not reach(92, L+1), not reach(16, L+1), not reach(85, L+1).
reach(16,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(40, L+1), not reach(68, L+1), not reach(79, L+1), not reach(92, L+1), not reach(22, L+1), not reach(85, L+1).
reach(85,L + 1) :- reach(69, L), L <= 10, L>= 3, not reach(3, L+1), not reach(36, L+1), not reach(40, L+1), not reach(68, L+1), not reach(79, L+1), not reach(92, L+1), not reach(22, L+1), not reach(16, L+1).
reach(3,L + 1) :- reach(38, L), L <= 8, L>= 2, not reach(24, L+1), not reach(44, L+1), not reach(54, L+1).
reach(24,L + 1) :- reach(38, L), L <= 8, L>= 2, not reach(3, L+1), not reach(44, L+1), not reach(54, L+1).
reach(44,L + 1) :- reach(38, L), L <= 8, L>= 2, not reach(3, L+1), not reach(24, L+1), not reach(54, L+1).
reach(54,L + 1) :- reach(38, L), L <= 8, L>= 2, not reach(3, L+1), not reach(24, L+1), not reach(44, L+1).
reach(68,L + 1) :- reach(4, L), L <= 12, L>= 5, not reach(51, L+1), not reach(6, L+1).
reach(51,L + 1) :- reach(4, L), L <= 12, L>= 5, not reach(68, L+1), not reach(6, L+1).
reach(6,L + 1) :- reach(4, L), L <= 12, L>= 5, not reach(68, L+1), not reach(51, L+1).
reach(11,L + 1) :- reach(6, L), L <= 11, L>= 5, not reach(55, L+1), not reach(85, L+1), not reach(4, L+1).
reach(55,L + 1) :- reach(6, L), L <= 11, L>= 5, not reach(11, L+1), not reach(85, L+1), not reach(4, L+1).
reach(85,L + 1) :- reach(6, L), L <= 11, L>= 5, not reach(11, L+1), not reach(55, L+1), not reach(4, L+1).
reach(4,L + 1) :- reach(6, L), L <= 11, L>= 5, not reach(11, L+1), not reach(55, L+1), not reach(85, L+1).
reach(6,L + 1) :- reach(11, L), L <= 10, L>= 4, not reach(63, L+1), not reach(92, L+1), not reach(28, L+1).
reach(63,L + 1) :- reach(11, L), L <= 10, L>= 4, not reach(6, L+1), not reach(92, L+1), not reach(28, L+1).
reach(92,L + 1) :- reach(11, L), L <= 10, L>= 4, not reach(6, L+1), not reach(63, L+1), not reach(28, L+1).
reach(28,L + 1) :- reach(11, L), L <= 10, L>= 4, not reach(6, L+1), not reach(63, L+1), not reach(92, L+1).
reach(6,L + 1) :- reach(55, L), L <= 10, L>= 4, not reach(29, L+1), not reach(47, L+1).
reach(29,L + 1) :- reach(55, L), L <= 10, L>= 4, not reach(6, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(55, L), L <= 10, L>= 4, not reach(6, L+1), not reach(29, L+1).
reach(6,L + 1) :- reach(85, L), L <= 10, L>= 4, not reach(42, L+1), not reach(69, L+1).
reach(42,L + 1) :- reach(85, L), L <= 10, L>= 4, not reach(6, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(85, L), L <= 10, L>= 4, not reach(6, L+1), not reach(42, L+1).
reach(14,L + 1) :- reach(7, L), L <= 8, L>= 4, not reach(29, L+1), not reach(40, L+1).
reach(29,L + 1) :- reach(7, L), L <= 8, L>= 4, not reach(14, L+1), not reach(40, L+1).
reach(40,L + 1) :- reach(7, L), L <= 8, L>= 4, not reach(14, L+1), not reach(29, L+1).
reach(7,L + 1) :- reach(14, L), L <= 8, L>= 5, not reach(72, L+1), not reach(47, L+1).
reach(72,L + 1) :- reach(14, L), L <= 8, L>= 5, not reach(7, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(14, L), L <= 8, L>= 5, not reach(7, L+1), not reach(72, L+1).
reach(69,L + 1) :- reach(68, L), L <= 11, L>= 4, not reach(1, L+1), not reach(4, L+1).
reach(1,L + 1) :- reach(68, L), L <= 11, L>= 4, not reach(69, L+1), not reach(4, L+1).
reach(4,L + 1) :- reach(68, L), L <= 11, L>= 4, not reach(69, L+1), not reach(1, L+1).
reach(10,L + 1) :- reach(42, L), L <= 9, L>= 3, not reach(13, L+1), not reach(54, L+1), not reach(85, L+1), not reach(92, L+1).
reach(13,L + 1) :- reach(42, L), L <= 9, L>= 3, not reach(10, L+1), not reach(54, L+1), not reach(85, L+1), not reach(92, L+1).
reach(54,L + 1) :- reach(42, L), L <= 9, L>= 3, not reach(10, L+1), not reach(13, L+1), not reach(85, L+1), not reach(92, L+1).
reach(85,L + 1) :- reach(42, L), L <= 9, L>= 3, not reach(10, L+1), not reach(13, L+1), not reach(54, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(42, L), L <= 9, L>= 3, not reach(10, L+1), not reach(13, L+1), not reach(54, L+1), not reach(85, L+1).
reach(11,L + 1) :- reach(63, L), L <= 9, L>= 4, not reach(70, L+1), not reach(92, L+1), not reach(21, L+1).
reach(70,L + 1) :- reach(63, L), L <= 9, L>= 4, not reach(11, L+1), not reach(92, L+1), not reach(21, L+1).
reach(92,L + 1) :- reach(63, L), L <= 9, L>= 4, not reach(11, L+1), not reach(70, L+1), not reach(21, L+1).
reach(21,L + 1) :- reach(63, L), L <= 9, L>= 4, not reach(11, L+1), not reach(70, L+1), not reach(92, L+1).
reach(11,L + 1) :- reach(92, L), L <= 9, L>= 3, not reach(42, L+1), not reach(57, L+1), not reach(63, L+1), not reach(69, L+1), not reach(88, L+1).
reach(42,L + 1) :- reach(92, L), L <= 9, L>= 3, not reach(11, L+1), not reach(57, L+1), not reach(63, L+1), not reach(69, L+1), not reach(88, L+1).
reach(57,L + 1) :- reach(92, L), L <= 9, L>= 3, not reach(11, L+1), not reach(42, L+1), not reach(63, L+1), not reach(69, L+1), not reach(88, L+1).
reach(63,L + 1) :- reach(92, L), L <= 9, L>= 3, not reach(11, L+1), not reach(42, L+1), not reach(57, L+1), not reach(69, L+1), not reach(88, L+1).
reach(69,L + 1) :- reach(92, L), L <= 9, L>= 3, not reach(11, L+1), not reach(42, L+1), not reach(57, L+1), not reach(63, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(92, L), L <= 9, L>= 3, not reach(11, L+1), not reach(42, L+1), not reach(57, L+1), not reach(63, L+1), not reach(69, L+1).
reach(24,L + 1) :- reach(29, L), L <= 9, L>= 3, not reach(36, L+1), not reach(55, L+1), not reach(79, L+1), not reach(7, L+1).
reach(36,L + 1) :- reach(29, L), L <= 9, L>= 3, not reach(24, L+1), not reach(55, L+1), not reach(79, L+1), not reach(7, L+1).
reach(55,L + 1) :- reach(29, L), L <= 9, L>= 3, not reach(24, L+1), not reach(36, L+1), not reach(79, L+1), not reach(7, L+1).
reach(79,L + 1) :- reach(29, L), L <= 9, L>= 3, not reach(24, L+1), not reach(36, L+1), not reach(55, L+1), not reach(7, L+1).
reach(7,L + 1) :- reach(29, L), L <= 9, L>= 3, not reach(24, L+1), not reach(36, L+1), not reach(55, L+1), not reach(79, L+1).
reach(44,L + 1) :- reach(13, L), L <= 8, L>= 2, not reach(79, L+1), not reach(56, L+1), not reach(42, L+1), not reach(50, L+1).
reach(79,L + 1) :- reach(13, L), L <= 8, L>= 2, not reach(44, L+1), not reach(56, L+1), not reach(42, L+1), not reach(50, L+1).
reach(56,L + 1) :- reach(13, L), L <= 8, L>= 2, not reach(44, L+1), not reach(79, L+1), not reach(42, L+1), not reach(50, L+1).
reach(42,L + 1) :- reach(13, L), L <= 8, L>= 2, not reach(44, L+1), not reach(79, L+1), not reach(56, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(13, L), L <= 8, L>= 2, not reach(44, L+1), not reach(79, L+1), not reach(56, L+1), not reach(42, L+1).
reach(13,L + 1) :- reach(44, L), L <= 7, L>= 1, not reach(38, L+1), not reach(18, L+1).
reach(38,L + 1) :- reach(44, L), L <= 7, L>= 1, not reach(13, L+1), not reach(18, L+1).
reach(13,L + 1) :- reach(79, L), L <= 9, L>= 3, not reach(29, L+1), not reach(40, L+1), not reach(69, L+1).
reach(29,L + 1) :- reach(79, L), L <= 9, L>= 3, not reach(13, L+1), not reach(40, L+1), not reach(69, L+1).
reach(40,L + 1) :- reach(79, L), L <= 9, L>= 3, not reach(13, L+1), not reach(29, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(79, L), L <= 9, L>= 3, not reach(13, L+1), not reach(29, L+1), not reach(40, L+1).
reach(13,L + 1) :- reach(56, L), L <= 8, L>= 3, not reach(16, L+1), not reach(21, L+1), not reach(73, L+1).
reach(16,L + 1) :- reach(56, L), L <= 8, L>= 3, not reach(13, L+1), not reach(21, L+1), not reach(73, L+1).
reach(21,L + 1) :- reach(56, L), L <= 8, L>= 3, not reach(13, L+1), not reach(16, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(56, L), L <= 8, L>= 3, not reach(13, L+1), not reach(16, L+1), not reach(21, L+1).
reach(14,L + 1) :- reach(72, L), L <= 7, L>= 4, not reach(35, L+1), not reach(80, L+1).
reach(35,L + 1) :- reach(72, L), L <= 7, L>= 4, not reach(14, L+1), not reach(80, L+1).
reach(80,L + 1) :- reach(72, L), L <= 7, L>= 4, not reach(14, L+1), not reach(35, L+1).
reach(57,L + 1) :- reach(20, L), L <= 8, L>= 1, not reach(73, L+1), not reach(18, L+1), not reach(16, L+1).
reach(73,L + 1) :- reach(20, L), L <= 8, L>= 1, not reach(57, L+1), not reach(18, L+1), not reach(16, L+1).
reach(16,L + 1) :- reach(20, L), L <= 8, L>= 1, not reach(57, L+1), not reach(73, L+1), not reach(18, L+1).
reach(20,L + 1) :- reach(18, L), L <= 7, L>= 0, not reach(84, L+1), not reach(44, L+1).
reach(84,L + 1) :- reach(18, L), L <= 7, L>= 0, not reach(20, L+1), not reach(44, L+1).
reach(44,L + 1) :- reach(18, L), L <= 7, L>= 0, not reach(20, L+1), not reach(84, L+1).
reach(56,L + 1) :- reach(16, L), L <= 9, L>= 2, not reach(50, L+1), not reach(77, L+1), not reach(20, L+1), not reach(69, L+1).
reach(50,L + 1) :- reach(16, L), L <= 9, L>= 2, not reach(56, L+1), not reach(77, L+1), not reach(20, L+1), not reach(69, L+1).
reach(77,L + 1) :- reach(16, L), L <= 9, L>= 2, not reach(56, L+1), not reach(50, L+1), not reach(20, L+1), not reach(69, L+1).
reach(20,L + 1) :- reach(16, L), L <= 9, L>= 2, not reach(56, L+1), not reach(50, L+1), not reach(77, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(16, L), L <= 9, L>= 2, not reach(56, L+1), not reach(50, L+1), not reach(77, L+1), not reach(20, L+1).
reach(16,L + 1) :- reach(50, L), L <= 8, L>= 3, not reach(48, L+1), not reach(13, L+1), not reach(47, L+1).
reach(48,L + 1) :- reach(50, L), L <= 8, L>= 3, not reach(16, L+1), not reach(13, L+1), not reach(47, L+1).
reach(13,L + 1) :- reach(50, L), L <= 8, L>= 3, not reach(16, L+1), not reach(48, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(50, L), L <= 8, L>= 3, not reach(16, L+1), not reach(48, L+1), not reach(13, L+1).
reach(27,L + 1) :- reach(77, L), L <= 8, L>= 3, not reach(70, L+1), not reach(80, L+1), not reach(16, L+1).
reach(70,L + 1) :- reach(77, L), L <= 8, L>= 3, not reach(27, L+1), not reach(80, L+1), not reach(16, L+1).
reach(80,L + 1) :- reach(77, L), L <= 8, L>= 3, not reach(27, L+1), not reach(70, L+1), not reach(16, L+1).
reach(16,L + 1) :- reach(77, L), L <= 8, L>= 3, not reach(27, L+1), not reach(70, L+1), not reach(80, L+1).
reach(20,L + 1) :- reach(73, L), L <= 8, L>= 2, not reach(36, L+1), not reach(56, L+1), not reach(59, L+1).
reach(36,L + 1) :- reach(73, L), L <= 8, L>= 2, not reach(20, L+1), not reach(56, L+1), not reach(59, L+1).
reach(56,L + 1) :- reach(73, L), L <= 8, L>= 2, not reach(20, L+1), not reach(36, L+1), not reach(59, L+1).
reach(59,L + 1) :- reach(73, L), L <= 8, L>= 2, not reach(20, L+1), not reach(36, L+1), not reach(56, L+1).
reach(56,L + 1) :- reach(21, L), L <= 8, L>= 3, not reach(57, L+1), not reach(61, L+1), not reach(63, L+1).
reach(57,L + 1) :- reach(21, L), L <= 8, L>= 3, not reach(56, L+1), not reach(61, L+1), not reach(63, L+1).
reach(61,L + 1) :- reach(21, L), L <= 8, L>= 3, not reach(56, L+1), not reach(57, L+1), not reach(63, L+1).
reach(63,L + 1) :- reach(21, L), L <= 8, L>= 3, not reach(56, L+1), not reach(57, L+1), not reach(61, L+1).
reach(21,L + 1) :- reach(61, L), L <= 7, L>= 3, not reach(80, L+1), not reach(57, L+1).
reach(80,L + 1) :- reach(61, L), L <= 7, L>= 3, not reach(21, L+1), not reach(57, L+1).
reach(57,L + 1) :- reach(61, L), L <= 7, L>= 3, not reach(21, L+1), not reach(80, L+1).
reach(57,L + 1) :- reach(22, L), L <= 9, L>= 3, not reach(49, L+1), not reach(69, L+1).
reach(49,L + 1) :- reach(22, L), L <= 9, L>= 3, not reach(57, L+1), not reach(69, L+1).
reach(69,L + 1) :- reach(22, L), L <= 9, L>= 3, not reach(57, L+1), not reach(49, L+1).
reach(22,L + 1) :- reach(49, L), L <= 8, L>= 4, not reach(48, L+1), not reach(70, L+1).
reach(48,L + 1) :- reach(49, L), L <= 8, L>= 4, not reach(22, L+1), not reach(70, L+1).
reach(70,L + 1) :- reach(49, L), L <= 8, L>= 4, not reach(22, L+1), not reach(48, L+1).
reach(29,L + 1) :- reach(24, L), L <= 8, L>= 2, not reach(30, L+1), not reach(38, L+1), not reach(84, L+1).
reach(30,L + 1) :- reach(24, L), L <= 8, L>= 2, not reach(29, L+1), not reach(38, L+1), not reach(84, L+1).
reach(38,L + 1) :- reach(24, L), L <= 8, L>= 2, not reach(29, L+1), not reach(30, L+1), not reach(84, L+1).
reach(84,L + 1) :- reach(24, L), L <= 8, L>= 2, not reach(29, L+1), not reach(30, L+1), not reach(38, L+1).
reach(24,L + 1) :- reach(30, L), L <= 8, L>= 3, not reach(57, L+1), not reach(40, L+1), not reach(70, L+1), not reach(43, L+1).
reach(57,L + 1) :- reach(30, L), L <= 8, L>= 3, not reach(24, L+1), not reach(40, L+1), not reach(70, L+1), not reach(43, L+1).
reach(40,L + 1) :- reach(30, L), L <= 8, L>= 3, not reach(24, L+1), not reach(57, L+1), not reach(70, L+1), not reach(43, L+1).
reach(70,L + 1) :- reach(30, L), L <= 8, L>= 3, not reach(24, L+1), not reach(57, L+1), not reach(40, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(30, L), L <= 8, L>= 3, not reach(24, L+1), not reach(57, L+1), not reach(40, L+1), not reach(70, L+1).
reach(24,L + 1) :- reach(84, L), L <= 8, L>= 1, not reach(45, L+1), not reach(18, L+1).
reach(45,L + 1) :- reach(84, L), L <= 8, L>= 1, not reach(24, L+1), not reach(18, L+1).
reach(59,L + 1) :- reach(27, L), L <= 7, L>= 4, not reach(77, L+1), not reach(54, L+1).
reach(77,L + 1) :- reach(27, L), L <= 7, L>= 4, not reach(59, L+1), not reach(54, L+1).
reach(54,L + 1) :- reach(27, L), L <= 7, L>= 4, not reach(59, L+1), not reach(77, L+1).
reach(34,L + 1) :- reach(54, L), L <= 8, L>= 3, not reach(42, L+1), not reach(27, L+1), not reach(38, L+1).
reach(42,L + 1) :- reach(54, L), L <= 8, L>= 3, not reach(34, L+1), not reach(27, L+1), not reach(38, L+1).
reach(27,L + 1) :- reach(54, L), L <= 8, L>= 3, not reach(34, L+1), not reach(42, L+1), not reach(38, L+1).
reach(38,L + 1) :- reach(54, L), L <= 8, L>= 3, not reach(34, L+1), not reach(42, L+1), not reach(27, L+1).
reach(27,L + 1) :- reach(59, L), L <= 8, L>= 3, not reach(28, L+1), not reach(73, L+1).
reach(28,L + 1) :- reach(59, L), L <= 8, L>= 3, not reach(27, L+1), not reach(73, L+1).
reach(73,L + 1) :- reach(59, L), L <= 8, L>= 3, not reach(27, L+1), not reach(28, L+1).
reach(35,L + 1) :- reach(28, L), L <= 9, L>= 3, not reach(59, L+1), not reach(11, L+1), not reach(45, L+1).
reach(59,L + 1) :- reach(28, L), L <= 9, L>= 3, not reach(35, L+1), not reach(11, L+1), not reach(45, L+1).
reach(11,L + 1) :- reach(28, L), L <= 9, L>= 3, not reach(35, L+1), not reach(59, L+1), not reach(45, L+1).
reach(45,L + 1) :- reach(28, L), L <= 9, L>= 3, not reach(35, L+1), not reach(59, L+1), not reach(11, L+1).
reach(28,L + 1) :- reach(35, L), L <= 8, L>= 3, not reach(43, L+1), not reach(45, L+1), not reach(72, L+1).
reach(43,L + 1) :- reach(35, L), L <= 8, L>= 3, not reach(28, L+1), not reach(45, L+1), not reach(72, L+1).
reach(45,L + 1) :- reach(35, L), L <= 8, L>= 3, not reach(28, L+1), not reach(43, L+1), not reach(72, L+1).
reach(72,L + 1) :- reach(35, L), L <= 8, L>= 3, not reach(28, L+1), not reach(43, L+1), not reach(45, L+1).
reach(29,L + 1) :- reach(36, L), L <= 9, L>= 3, not reach(69, L+1), not reach(73, L+1), not reach(43, L+1).
reach(69,L + 1) :- reach(36, L), L <= 9, L>= 3, not reach(29, L+1), not reach(73, L+1), not reach(43, L+1).
reach(73,L + 1) :- reach(36, L), L <= 9, L>= 3, not reach(29, L+1), not reach(69, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(36, L), L <= 9, L>= 3, not reach(29, L+1), not reach(69, L+1), not reach(73, L+1).
reach(30,L + 1) :- reach(40, L), L <= 9, L>= 4, not reach(69, L+1), not reach(79, L+1), not reach(7, L+1).
reach(69,L + 1) :- reach(40, L), L <= 9, L>= 4, not reach(30, L+1), not reach(79, L+1), not reach(7, L+1).
reach(79,L + 1) :- reach(40, L), L <= 9, L>= 4, not reach(30, L+1), not reach(69, L+1), not reach(7, L+1).
reach(7,L + 1) :- reach(40, L), L <= 9, L>= 4, not reach(30, L+1), not reach(69, L+1), not reach(79, L+1).
reach(30,L + 1) :- reach(70, L), L <= 8, L>= 4, not reach(49, L+1), not reach(63, L+1), not reach(77, L+1).
reach(49,L + 1) :- reach(70, L), L <= 8, L>= 4, not reach(30, L+1), not reach(63, L+1), not reach(77, L+1).
reach(63,L + 1) :- reach(70, L), L <= 8, L>= 4, not reach(30, L+1), not reach(49, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(70, L), L <= 8, L>= 4, not reach(30, L+1), not reach(49, L+1), not reach(63, L+1).
reach(54,L + 1) :- reach(34, L), L <= 8, L>= 4, not reach(48, L+1), not reach(88, L+1).
reach(48,L + 1) :- reach(34, L), L <= 8, L>= 4, not reach(54, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(34, L), L <= 8, L>= 4, not reach(54, L+1), not reach(48, L+1).
reach(34,L + 1) :- reach(48, L), L <= 7, L>= 4, not reach(49, L+1), not reach(50, L+1).
reach(49,L + 1) :- reach(48, L), L <= 7, L>= 4, not reach(34, L+1), not reach(50, L+1).
reach(50,L + 1) :- reach(48, L), L <= 7, L>= 4, not reach(34, L+1), not reach(49, L+1).
reach(35,L + 1) :- reach(43, L), L <= 8, L>= 4, not reach(36, L+1), not reach(30, L+1).
reach(36,L + 1) :- reach(43, L), L <= 8, L>= 4, not reach(35, L+1), not reach(30, L+1).
reach(30,L + 1) :- reach(43, L), L <= 8, L>= 4, not reach(35, L+1), not reach(36, L+1).
reach(55,L + 1) :- reach(47, L), L <= 9, L>= 4, not reach(50, L+1), not reach(14, L+1).
reach(50,L + 1) :- reach(47, L), L <= 9, L>= 4, not reach(55, L+1), not reach(14, L+1).
reach(14,L + 1) :- reach(47, L), L <= 9, L>= 4, not reach(55, L+1), not reach(50, L+1).
reach(61,L + 1) :- reach(80, L), L <= 7, L>= 4, not reach(72, L+1), not reach(77, L+1).
reach(72,L + 1) :- reach(80, L), L <= 7, L>= 4, not reach(61, L+1), not reach(77, L+1).
reach(77,L + 1) :- reach(80, L), L <= 7, L>= 4, not reach(61, L+1), not reach(72, L+1).

