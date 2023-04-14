reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17), not reach(X, 18), not reach(X, 19), not reach(X, 20), not reach(X, 21), not reach(X, 22), not reach(X, 23), not reach(X, 24), not reach(X, 25), not reach(X, 26), not reach(X, 27), not reach(X, 28), not reach(X, 29), not reach(X, 30), not reach(X, 31), not reach(X, 32), not reach(X, 33), not reach(X, 34), not reach(X, 35), not reach(X, 36), not reach(X, 37), not reach(X, 38), not reach(X, 39), not reach(X, 40), not reach(X, 41), not reach(X, 42), not reach(X, 43), not reach(X, 44), not reach(X, 45), not reach(X, 46), not reach(X, 47), not reach(X, 48), not reach(X, 49), not reach(X, 50), not reach(X, 51), not reach(X, 52), not reach(X, 53), not reach(X, 54), not reach(X, 55), not reach(X, 56), not reach(X, 57), not reach(X, 58), not reach(X, 59), not reach(X, 60), not reach(X, 61).
start(25).
goal(28).
:- reach(X, L), reach(X, L'), L != L'.
reach(2,L + 1) :- reach(1, L), L <= 30, L>= 28, not reach(3, L+1).
reach(3,L + 1) :- reach(1, L), L <= 30, L>= 28, not reach(2, L+1).
reach(1,L + 1) :- reach(2, L), L <= 29, L>= 27, not reach(221, L+1), not reach(272, L+1).
reach(221,L + 1) :- reach(2, L), L <= 29, L>= 27, not reach(1, L+1), not reach(272, L+1).
reach(272,L + 1) :- reach(2, L), L <= 29, L>= 27, not reach(1, L+1), not reach(221, L+1).
reach(1,L + 1) :- reach(3, L), L <= 31, L>= 29, not reach(684, L+1).
reach(684,L + 1) :- reach(3, L), L <= 31, L>= 29, not reach(1, L+1).
reach(2,L + 1) :- reach(221, L), L <= 28, L>= 28, not reach(10, L+1).
reach(10,L + 1) :- reach(221, L), L <= 28, L>= 28, not reach(2, L+1).
reach(2,L + 1) :- reach(272, L), L <= 28, L>= 26, not reach(47, L+1), not reach(271, L+1).
reach(47,L + 1) :- reach(272, L), L <= 28, L>= 26, not reach(2, L+1), not reach(271, L+1).
reach(271,L + 1) :- reach(272, L), L <= 28, L>= 26, not reach(2, L+1), not reach(47, L+1).
reach(3,L + 1) :- reach(684, L), L <= 32, L>= 30, not reach(408, L+1).
reach(408,L + 1) :- reach(684, L), L <= 32, L>= 30, not reach(3, L+1).
reach(8,L + 1) :- reach(7, L), L <= 37, L>= 35, not reach(9, L+1).
reach(9,L + 1) :- reach(7, L), L <= 37, L>= 35, not reach(8, L+1).
reach(7,L + 1) :- reach(8, L), L <= 36, L>= 34, not reach(11, L+1).
reach(11,L + 1) :- reach(8, L), L <= 36, L>= 34, not reach(7, L+1).
reach(7,L + 1) :- reach(9, L), L <= 38, L>= 36, not reach(599, L+1).
reach(599,L + 1) :- reach(9, L), L <= 38, L>= 36, not reach(7, L+1).
reach(8,L + 1) :- reach(11, L), L <= 35, L>= 33, not reach(51, L+1).
reach(51,L + 1) :- reach(11, L), L <= 35, L>= 33, not reach(8, L+1).
reach(9,L + 1) :- reach(599, L), L <= 39, L>= 37, not reach(564, L+1).
reach(564,L + 1) :- reach(599, L), L <= 39, L>= 37, not reach(9, L+1).
reach(11,L + 1) :- reach(51, L), L <= 34, L>= 32, not reach(49, L+1), not reach(408, L+1).
reach(49,L + 1) :- reach(51, L), L <= 34, L>= 32, not reach(11, L+1), not reach(408, L+1).
reach(408,L + 1) :- reach(51, L), L <= 34, L>= 32, not reach(11, L+1), not reach(49, L+1).
reach(13,L + 1) :- reach(12, L), L <= 12, L>= 10, not reach(14, L+1).
reach(14,L + 1) :- reach(12, L), L <= 12, L>= 10, not reach(13, L+1).
reach(12,L + 1) :- reach(13, L), L <= 13, L>= 11, not reach(487, L+1).
reach(487,L + 1) :- reach(13, L), L <= 13, L>= 11, not reach(12, L+1).
reach(12,L + 1) :- reach(14, L), L <= 12, L>= 9, not reach(15, L+1), not reach(338, L+1), not reach(745, L+1).
reach(15,L + 1) :- reach(14, L), L <= 12, L>= 9, not reach(12, L+1), not reach(338, L+1), not reach(745, L+1).
reach(338,L + 1) :- reach(14, L), L <= 12, L>= 9, not reach(12, L+1), not reach(15, L+1), not reach(745, L+1).
reach(745,L + 1) :- reach(14, L), L <= 12, L>= 9, not reach(12, L+1), not reach(15, L+1), not reach(338, L+1).
reach(13,L + 1) :- reach(487, L), L <= 14, L>= 12, not reach(427, L+1), not reach(428, L+1), not reach(489, L+1).
reach(427,L + 1) :- reach(487, L), L <= 14, L>= 12, not reach(13, L+1), not reach(428, L+1), not reach(489, L+1).
reach(428,L + 1) :- reach(487, L), L <= 14, L>= 12, not reach(13, L+1), not reach(427, L+1), not reach(489, L+1).
reach(489,L + 1) :- reach(487, L), L <= 14, L>= 12, not reach(13, L+1), not reach(427, L+1), not reach(428, L+1).
reach(14,L + 1) :- reach(15, L), L <= 13, L>= 10, not reach(16, L+1).
reach(16,L + 1) :- reach(15, L), L <= 13, L>= 10, not reach(14, L+1).
reach(14,L + 1) :- reach(338, L), L <= 11, L>= 10, not reach(335, L+1).
reach(335,L + 1) :- reach(338, L), L <= 11, L>= 10, not reach(14, L+1).
reach(14,L + 1) :- reach(745, L), L <= 11, L>= 8, not reach(744, L+1).
reach(744,L + 1) :- reach(745, L), L <= 11, L>= 8, not reach(14, L+1).
reach(15,L + 1) :- reach(16, L), L <= 14, L>= 11, not reach(482, L+1).
reach(482,L + 1) :- reach(16, L), L <= 14, L>= 11, not reach(15, L+1).
reach(16,L + 1) :- reach(482, L), L <= 15, L>= 12, not reach(658, L+1).
reach(658,L + 1) :- reach(482, L), L <= 15, L>= 12, not reach(16, L+1).
reach(18,L + 1) :- reach(17, L), L <= 19, L>= 18, not reach(19, L+1).
reach(19,L + 1) :- reach(17, L), L <= 19, L>= 18, not reach(18, L+1).
reach(17,L + 1) :- reach(19, L), L <= 20, L>= 17, not reach(208, L+1), not reach(401, L+1), not reach(405, L+1), not reach(657, L+1), not reach(360, L+1).
reach(208,L + 1) :- reach(19, L), L <= 20, L>= 17, not reach(17, L+1), not reach(401, L+1), not reach(405, L+1), not reach(657, L+1), not reach(360, L+1).
reach(401,L + 1) :- reach(19, L), L <= 20, L>= 17, not reach(17, L+1), not reach(208, L+1), not reach(405, L+1), not reach(657, L+1), not reach(360, L+1).
reach(405,L + 1) :- reach(19, L), L <= 20, L>= 17, not reach(17, L+1), not reach(208, L+1), not reach(401, L+1), not reach(657, L+1), not reach(360, L+1).
reach(657,L + 1) :- reach(19, L), L <= 20, L>= 17, not reach(17, L+1), not reach(208, L+1), not reach(401, L+1), not reach(405, L+1), not reach(360, L+1).
reach(360,L + 1) :- reach(19, L), L <= 20, L>= 17, not reach(17, L+1), not reach(208, L+1), not reach(401, L+1), not reach(405, L+1), not reach(657, L+1).
reach(19,L + 1) :- reach(208, L), L <= 20, L>= 18, not reach(207, L+1).
reach(207,L + 1) :- reach(208, L), L <= 20, L>= 18, not reach(19, L+1).
reach(19,L + 1) :- reach(401, L), L <= 19, L>= 18, not reach(400, L+1).
reach(400,L + 1) :- reach(401, L), L <= 19, L>= 18, not reach(19, L+1).
reach(19,L + 1) :- reach(405, L), L <= 21, L>= 18, not reach(404, L+1).
reach(404,L + 1) :- reach(405, L), L <= 21, L>= 18, not reach(19, L+1).
reach(19,L + 1) :- reach(657, L), L <= 19, L>= 17, not reach(268, L+1), not reach(362, L+1), not reach(656, L+1).
reach(268,L + 1) :- reach(657, L), L <= 19, L>= 17, not reach(19, L+1), not reach(362, L+1), not reach(656, L+1).
reach(362,L + 1) :- reach(657, L), L <= 19, L>= 17, not reach(19, L+1), not reach(268, L+1), not reach(656, L+1).
reach(656,L + 1) :- reach(657, L), L <= 19, L>= 17, not reach(19, L+1), not reach(268, L+1), not reach(362, L+1).
reach(19,L + 1) :- reach(360, L), L <= 19, L>= 16, not reach(359, L+1).
reach(359,L + 1) :- reach(360, L), L <= 19, L>= 16, not reach(19, L+1).
reach(21,L + 1) :- reach(20, L), L <= 48, L>= 46, not reach(22, L+1).
reach(22,L + 1) :- reach(20, L), L <= 48, L>= 46, not reach(21, L+1).
reach(20,L + 1) :- reach(21, L), L <= 49, L>= 47, not reach(23, L+1).
reach(23,L + 1) :- reach(21, L), L <= 49, L>= 47, not reach(20, L+1).
reach(20,L + 1) :- reach(22, L), L <= 47, L>= 45, not reach(605, L+1).
reach(605,L + 1) :- reach(22, L), L <= 47, L>= 45, not reach(20, L+1).
reach(21,L + 1) :- reach(23, L), L <= 50, L>= 48, not reach(385, L+1).
reach(385,L + 1) :- reach(23, L), L <= 50, L>= 48, not reach(21, L+1).
reach(22,L + 1) :- reach(605, L), L <= 46, L>= 44, not reach(604, L+1).
reach(604,L + 1) :- reach(605, L), L <= 46, L>= 44, not reach(22, L+1).
reach(23,L + 1) :- reach(385, L), L <= 51, L>= 49, not reach(383, L+1).
reach(383,L + 1) :- reach(385, L), L <= 51, L>= 49, not reach(23, L+1).
reach(26,L + 1) :- reach(24, L), L <= 4, L>= 1, not reach(25, L+1).
reach(24,L + 1) :- reach(25, L), L <= 3, L>= 0, not reach(27, L+1).
reach(27,L + 1) :- reach(25, L), L <= 3, L>= 0, not reach(24, L+1).
reach(24,L + 1) :- reach(26, L), L <= 5, L>= 2, not reach(443, L+1).
reach(443,L + 1) :- reach(26, L), L <= 5, L>= 2, not reach(24, L+1).
reach(31,L + 1) :- reach(27, L), L <= 4, L>= 1, not reach(25, L+1).
reach(26,L + 1) :- reach(443, L), L <= 6, L>= 3, not reach(36, L+1).
reach(36,L + 1) :- reach(443, L), L <= 6, L>= 3, not reach(26, L+1).
reach(27,L + 1) :- reach(31, L), L <= 5, L>= 2, not reach(32, L+1).
reach(32,L + 1) :- reach(31, L), L <= 5, L>= 2, not reach(27, L+1).
reach(28,L + 1) :- reach(29, L), L <= 60, L>= 58, not reach(144, L+1).
reach(144,L + 1) :- reach(29, L), L <= 60, L>= 58, not reach(28, L+1).
reach(28,L + 1) :- reach(30, L), L <= 60, L>= 57, not reach(115, L+1).
reach(115,L + 1) :- reach(30, L), L <= 60, L>= 57, not reach(28, L+1).
reach(29,L + 1) :- reach(144, L), L <= 59, L>= 57, not reach(685, L+1).
reach(685,L + 1) :- reach(144, L), L <= 59, L>= 57, not reach(29, L+1).
reach(30,L + 1) :- reach(115, L), L <= 59, L>= 56, not reach(116, L+1).
reach(116,L + 1) :- reach(115, L), L <= 59, L>= 56, not reach(30, L+1).
reach(31,L + 1) :- reach(32, L), L <= 6, L>= 3, not reach(33, L+1).
reach(33,L + 1) :- reach(32, L), L <= 6, L>= 3, not reach(31, L+1).
reach(32,L + 1) :- reach(33, L), L <= 7, L>= 4, not reach(34, L+1).
reach(34,L + 1) :- reach(33, L), L <= 7, L>= 4, not reach(32, L+1).
reach(33,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(35, L+1), not reach(36, L+1), not reach(37, L+1).
reach(35,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(33, L+1), not reach(36, L+1), not reach(37, L+1).
reach(36,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(33, L+1), not reach(35, L+1), not reach(37, L+1).
reach(37,L + 1) :- reach(34, L), L <= 8, L>= 5, not reach(33, L+1), not reach(35, L+1), not reach(36, L+1).
reach(34,L + 1) :- reach(35, L), L <= 7, L>= 6, not reach(38, L+1).
reach(38,L + 1) :- reach(35, L), L <= 7, L>= 6, not reach(34, L+1).
reach(34,L + 1) :- reach(36, L), L <= 7, L>= 4, not reach(485, L+1), not reach(443, L+1).
reach(485,L + 1) :- reach(36, L), L <= 7, L>= 4, not reach(34, L+1), not reach(443, L+1).
reach(443,L + 1) :- reach(36, L), L <= 7, L>= 4, not reach(34, L+1), not reach(485, L+1).
reach(34,L + 1) :- reach(37, L), L <= 9, L>= 6, not reach(744, L+1).
reach(744,L + 1) :- reach(37, L), L <= 9, L>= 6, not reach(34, L+1).
reach(35,L + 1) :- reach(38, L), L <= 8, L>= 7, not reach(339, L+1).
reach(339,L + 1) :- reach(38, L), L <= 8, L>= 7, not reach(35, L+1).
reach(36,L + 1) :- reach(485, L), L <= 6, L>= 5, not reach(486, L+1).
reach(486,L + 1) :- reach(485, L), L <= 6, L>= 5, not reach(36, L+1).
reach(37,L + 1) :- reach(744, L), L <= 10, L>= 7, not reach(745, L+1).
reach(745,L + 1) :- reach(744, L), L <= 10, L>= 7, not reach(37, L+1).
reach(38,L + 1) :- reach(339, L), L <= 9, L>= 8, not reach(336, L+1).
reach(336,L + 1) :- reach(339, L), L <= 9, L>= 8, not reach(38, L+1).
reach(40,L + 1) :- reach(39, L), L <= 31, L>= 28, not reach(41, L+1), not reach(42, L+1), not reach(43, L+1).
reach(41,L + 1) :- reach(39, L), L <= 31, L>= 28, not reach(40, L+1), not reach(42, L+1), not reach(43, L+1).
reach(42,L + 1) :- reach(39, L), L <= 31, L>= 28, not reach(40, L+1), not reach(41, L+1), not reach(43, L+1).
reach(43,L + 1) :- reach(39, L), L <= 31, L>= 28, not reach(40, L+1), not reach(41, L+1), not reach(42, L+1).
reach(39,L + 1) :- reach(40, L), L <= 30, L>= 27, not reach(98, L+1).
reach(98,L + 1) :- reach(40, L), L <= 30, L>= 27, not reach(39, L+1).
reach(39,L + 1) :- reach(41, L), L <= 30, L>= 29, not reach(44, L+1).
reach(44,L + 1) :- reach(41, L), L <= 30, L>= 29, not reach(39, L+1).
reach(39,L + 1) :- reach(42, L), L <= 30, L>= 29, not reach(141, L+1), not reach(513, L+1).
reach(141,L + 1) :- reach(42, L), L <= 30, L>= 29, not reach(39, L+1), not reach(513, L+1).
reach(513,L + 1) :- reach(42, L), L <= 30, L>= 29, not reach(39, L+1), not reach(141, L+1).
reach(39,L + 1) :- reach(43, L), L <= 32, L>= 29, not reach(59, L+1), not reach(644, L+1), not reach(645, L+1).
reach(59,L + 1) :- reach(43, L), L <= 32, L>= 29, not reach(39, L+1), not reach(644, L+1), not reach(645, L+1).
reach(644,L + 1) :- reach(43, L), L <= 32, L>= 29, not reach(39, L+1), not reach(59, L+1), not reach(645, L+1).
reach(645,L + 1) :- reach(43, L), L <= 32, L>= 29, not reach(39, L+1), not reach(59, L+1), not reach(644, L+1).
reach(40,L + 1) :- reach(98, L), L <= 29, L>= 26, not reach(425, L+1), not reach(303, L+1).
reach(425,L + 1) :- reach(98, L), L <= 29, L>= 26, not reach(40, L+1), not reach(303, L+1).
reach(303,L + 1) :- reach(98, L), L <= 29, L>= 26, not reach(40, L+1), not reach(425, L+1).
reach(41,L + 1) :- reach(44, L), L <= 30, L>= 30, not reach(45, L+1).
reach(45,L + 1) :- reach(44, L), L <= 30, L>= 30, not reach(41, L+1).
reach(42,L + 1) :- reach(141, L), L <= 29, L>= 28, not reach(142, L+1), not reach(143, L+1).
reach(142,L + 1) :- reach(141, L), L <= 29, L>= 28, not reach(42, L+1), not reach(143, L+1).
reach(143,L + 1) :- reach(141, L), L <= 29, L>= 28, not reach(42, L+1), not reach(142, L+1).
reach(43,L + 1) :- reach(59, L), L <= 33, L>= 30, not reach(60, L+1), not reach(61, L+1).
reach(60,L + 1) :- reach(59, L), L <= 33, L>= 30, not reach(43, L+1), not reach(61, L+1).
reach(61,L + 1) :- reach(59, L), L <= 33, L>= 30, not reach(43, L+1), not reach(60, L+1).
reach(43,L + 1) :- reach(644, L), L <= 31, L>= 30, not reach(645, L+1).
reach(645,L + 1) :- reach(644, L), L <= 31, L>= 30, not reach(43, L+1).
reach(43,L + 1) :- reach(645, L), L <= 31, L>= 30, not reach(644, L+1), not reach(646, L+1).
reach(644,L + 1) :- reach(645, L), L <= 31, L>= 30, not reach(43, L+1), not reach(646, L+1).
reach(646,L + 1) :- reach(645, L), L <= 31, L>= 30, not reach(43, L+1), not reach(644, L+1).
reach(44,L + 1) :- reach(45, L), L <= 31, L>= 31, not reach(307, L+1).
reach(307,L + 1) :- reach(45, L), L <= 31, L>= 31, not reach(44, L+1).
reach(45,L + 1) :- reach(307, L), L <= 32, L>= 30, not reach(107, L+1), not reach(309, L+1).
reach(107,L + 1) :- reach(307, L), L <= 32, L>= 30, not reach(45, L+1), not reach(309, L+1).
reach(309,L + 1) :- reach(307, L), L <= 32, L>= 30, not reach(45, L+1), not reach(107, L+1).
reach(47,L + 1) :- reach(46, L), L <= 26, L>= 25, not reach(48, L+1).
reach(48,L + 1) :- reach(46, L), L <= 26, L>= 25, not reach(47, L+1).
reach(46,L + 1) :- reach(47, L), L <= 27, L>= 25, not reach(272, L+1), not reach(48, L+1), not reach(275, L+1).
reach(272,L + 1) :- reach(47, L), L <= 27, L>= 25, not reach(46, L+1), not reach(48, L+1), not reach(275, L+1).
reach(48,L + 1) :- reach(47, L), L <= 27, L>= 25, not reach(46, L+1), not reach(272, L+1), not reach(275, L+1).
reach(275,L + 1) :- reach(47, L), L <= 27, L>= 25, not reach(46, L+1), not reach(272, L+1), not reach(48, L+1).
reach(46,L + 1) :- reach(48, L), L <= 26, L>= 24, not reach(47, L+1), not reach(666, L+1).
reach(47,L + 1) :- reach(48, L), L <= 26, L>= 24, not reach(46, L+1), not reach(666, L+1).
reach(666,L + 1) :- reach(48, L), L <= 26, L>= 24, not reach(46, L+1), not reach(47, L+1).
reach(47,L + 1) :- reach(275, L), L <= 26, L>= 26, not reach(271, L+1).
reach(271,L + 1) :- reach(275, L), L <= 26, L>= 26, not reach(47, L+1).
reach(48,L + 1) :- reach(666, L), L <= 25, L>= 23, not reach(493, L+1).
reach(493,L + 1) :- reach(666, L), L <= 25, L>= 23, not reach(48, L+1).
reach(50,L + 1) :- reach(49, L), L <= 33, L>= 33, not reach(51, L+1).
reach(51,L + 1) :- reach(49, L), L <= 33, L>= 33, not reach(50, L+1).
reach(51,L + 1) :- reach(408, L), L <= 33, L>= 31, not reach(684, L+1).
reach(684,L + 1) :- reach(408, L), L <= 33, L>= 31, not reach(51, L+1).
reach(59,L + 1) :- reach(60, L), L <= 34, L>= 31, not reach(108, L+1).
reach(108,L + 1) :- reach(60, L), L <= 34, L>= 31, not reach(59, L+1).
reach(59,L + 1) :- reach(61, L), L <= 34, L>= 31, not reach(104, L+1), not reach(109, L+1).
reach(104,L + 1) :- reach(61, L), L <= 34, L>= 31, not reach(59, L+1), not reach(109, L+1).
reach(109,L + 1) :- reach(61, L), L <= 34, L>= 31, not reach(59, L+1), not reach(104, L+1).
reach(60,L + 1) :- reach(108, L), L <= 35, L>= 32, not reach(88, L+1), not reach(106, L+1).
reach(88,L + 1) :- reach(108, L), L <= 35, L>= 32, not reach(60, L+1), not reach(106, L+1).
reach(106,L + 1) :- reach(108, L), L <= 35, L>= 32, not reach(60, L+1), not reach(88, L+1).
reach(61,L + 1) :- reach(104, L), L <= 34, L>= 32, not reach(102, L+1), not reach(109, L+1).
reach(102,L + 1) :- reach(104, L), L <= 34, L>= 32, not reach(61, L+1), not reach(109, L+1).
reach(109,L + 1) :- reach(104, L), L <= 34, L>= 32, not reach(61, L+1), not reach(102, L+1).
reach(61,L + 1) :- reach(109, L), L <= 35, L>= 32, not reach(104, L+1), not reach(374, L+1).
reach(104,L + 1) :- reach(109, L), L <= 35, L>= 32, not reach(61, L+1), not reach(374, L+1).
reach(374,L + 1) :- reach(109, L), L <= 35, L>= 32, not reach(61, L+1), not reach(104, L+1).
reach(67,L + 1) :- reach(68, L), L <= 20, L>= 20, not reach(70, L+1).
reach(70,L + 1) :- reach(68, L), L <= 20, L>= 20, not reach(67, L+1).
reach(68,L + 1) :- reach(70, L), L <= 21, L>= 19, not reach(77, L+1), not reach(79, L+1), not reach(80, L+1), not reach(81, L+1).
reach(77,L + 1) :- reach(70, L), L <= 21, L>= 19, not reach(68, L+1), not reach(79, L+1), not reach(80, L+1), not reach(81, L+1).
reach(79,L + 1) :- reach(70, L), L <= 21, L>= 19, not reach(68, L+1), not reach(77, L+1), not reach(80, L+1), not reach(81, L+1).
reach(80,L + 1) :- reach(70, L), L <= 21, L>= 19, not reach(68, L+1), not reach(77, L+1), not reach(79, L+1), not reach(81, L+1).
reach(81,L + 1) :- reach(70, L), L <= 21, L>= 19, not reach(68, L+1), not reach(77, L+1), not reach(79, L+1), not reach(80, L+1).
reach(70,L + 1) :- reach(77, L), L <= 20, L>= 19, not reach(75, L+1).
reach(75,L + 1) :- reach(77, L), L <= 20, L>= 19, not reach(70, L+1).
reach(70,L + 1) :- reach(79, L), L <= 20, L>= 18, not reach(90, L+1).
reach(90,L + 1) :- reach(79, L), L <= 20, L>= 18, not reach(70, L+1).
reach(70,L + 1) :- reach(80, L), L <= 22, L>= 20, not reach(663, L+1).
reach(663,L + 1) :- reach(80, L), L <= 22, L>= 20, not reach(70, L+1).
reach(70,L + 1) :- reach(81, L), L <= 20, L>= 20, not reach(403, L+1).
reach(403,L + 1) :- reach(81, L), L <= 20, L>= 20, not reach(70, L+1).
reach(72,L + 1) :- reach(71, L), L <= 36, L>= 35, not reach(73, L+1).
reach(73,L + 1) :- reach(71, L), L <= 36, L>= 35, not reach(72, L+1).
reach(71,L + 1) :- reach(72, L), L <= 36, L>= 36, not reach(95, L+1).
reach(95,L + 1) :- reach(72, L), L <= 36, L>= 36, not reach(71, L+1).
reach(71,L + 1) :- reach(73, L), L <= 37, L>= 34, not reach(74, L+1), not reach(88, L+1), not reach(391, L+1), not reach(392, L+1).
reach(74,L + 1) :- reach(73, L), L <= 37, L>= 34, not reach(71, L+1), not reach(88, L+1), not reach(391, L+1), not reach(392, L+1).
reach(88,L + 1) :- reach(73, L), L <= 37, L>= 34, not reach(71, L+1), not reach(74, L+1), not reach(391, L+1), not reach(392, L+1).
reach(391,L + 1) :- reach(73, L), L <= 37, L>= 34, not reach(71, L+1), not reach(74, L+1), not reach(88, L+1), not reach(392, L+1).
reach(392,L + 1) :- reach(73, L), L <= 37, L>= 34, not reach(71, L+1), not reach(74, L+1), not reach(88, L+1), not reach(391, L+1).
reach(72,L + 1) :- reach(95, L), L <= 37, L>= 37, not reach(93, L+1).
reach(93,L + 1) :- reach(95, L), L <= 37, L>= 37, not reach(72, L+1).
reach(73,L + 1) :- reach(74, L), L <= 36, L>= 35.
reach(73,L + 1) :- reach(88, L), L <= 36, L>= 33, not reach(87, L+1), not reach(100, L+1), not reach(108, L+1), not reach(106, L+1), not reach(110, L+1).
reach(87,L + 1) :- reach(88, L), L <= 36, L>= 33, not reach(73, L+1), not reach(100, L+1), not reach(108, L+1), not reach(106, L+1), not reach(110, L+1).
reach(100,L + 1) :- reach(88, L), L <= 36, L>= 33, not reach(73, L+1), not reach(87, L+1), not reach(108, L+1), not reach(106, L+1), not reach(110, L+1).
reach(108,L + 1) :- reach(88, L), L <= 36, L>= 33, not reach(73, L+1), not reach(87, L+1), not reach(100, L+1), not reach(106, L+1), not reach(110, L+1).
reach(106,L + 1) :- reach(88, L), L <= 36, L>= 33, not reach(73, L+1), not reach(87, L+1), not reach(100, L+1), not reach(108, L+1), not reach(110, L+1).
reach(110,L + 1) :- reach(88, L), L <= 36, L>= 33, not reach(73, L+1), not reach(87, L+1), not reach(100, L+1), not reach(108, L+1), not reach(106, L+1).
reach(73,L + 1) :- reach(391, L), L <= 36, L>= 35.
reach(73,L + 1) :- reach(392, L), L <= 38, L>= 35, not reach(387, L+1), not reach(393, L+1).
reach(387,L + 1) :- reach(392, L), L <= 38, L>= 35, not reach(73, L+1), not reach(393, L+1).
reach(393,L + 1) :- reach(392, L), L <= 38, L>= 35, not reach(73, L+1), not reach(387, L+1).
reach(76,L + 1) :- reach(75, L), L <= 19, L>= 18, not reach(77, L+1).
reach(77,L + 1) :- reach(75, L), L <= 19, L>= 18, not reach(76, L+1).
reach(75,L + 1) :- reach(76, L), L <= 18, L>= 17, not reach(78, L+1).
reach(78,L + 1) :- reach(76, L), L <= 18, L>= 17, not reach(75, L+1).
reach(76,L + 1) :- reach(78, L), L <= 17, L>= 16, not reach(200, L+1).
reach(200,L + 1) :- reach(78, L), L <= 17, L>= 16, not reach(76, L+1).
reach(78,L + 1) :- reach(200, L), L <= 16, L>= 15, not reach(194, L+1).
reach(194,L + 1) :- reach(200, L), L <= 16, L>= 15, not reach(78, L+1).
reach(79,L + 1) :- reach(90, L), L <= 19, L>= 17, not reach(89, L+1), not reach(92, L+1).
reach(89,L + 1) :- reach(90, L), L <= 19, L>= 17, not reach(79, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(90, L), L <= 19, L>= 17, not reach(79, L+1), not reach(89, L+1).
reach(80,L + 1) :- reach(663, L), L <= 23, L>= 21, not reach(493, L+1), not reach(661, L+1), not reach(740, L+1).
reach(493,L + 1) :- reach(663, L), L <= 23, L>= 21, not reach(80, L+1), not reach(661, L+1), not reach(740, L+1).
reach(661,L + 1) :- reach(663, L), L <= 23, L>= 21, not reach(80, L+1), not reach(493, L+1), not reach(740, L+1).
reach(740,L + 1) :- reach(663, L), L <= 23, L>= 21, not reach(80, L+1), not reach(493, L+1), not reach(661, L+1).
reach(83,L + 1) :- reach(82, L), L <= 38, L>= 36, not reach(84, L+1).
reach(84,L + 1) :- reach(82, L), L <= 38, L>= 36, not reach(83, L+1).
reach(82,L + 1) :- reach(83, L), L <= 39, L>= 37, not reach(85, L+1), not reach(86, L+1).
reach(85,L + 1) :- reach(83, L), L <= 39, L>= 37, not reach(82, L+1), not reach(86, L+1).
reach(86,L + 1) :- reach(83, L), L <= 39, L>= 37, not reach(82, L+1), not reach(85, L+1).
reach(82,L + 1) :- reach(84, L), L <= 37, L>= 35, not reach(87, L+1).
reach(87,L + 1) :- reach(84, L), L <= 37, L>= 35, not reach(82, L+1).
reach(83,L + 1) :- reach(85, L), L <= 39, L>= 38, not reach(93, L+1), not reach(94, L+1).
reach(93,L + 1) :- reach(85, L), L <= 39, L>= 38, not reach(83, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(85, L), L <= 39, L>= 38, not reach(83, L+1), not reach(93, L+1).
reach(83,L + 1) :- reach(86, L), L <= 40, L>= 37, not reach(152, L+1), not reach(414, L+1), not reach(433, L+1), not reach(704, L+1).
reach(152,L + 1) :- reach(86, L), L <= 40, L>= 37, not reach(83, L+1), not reach(414, L+1), not reach(433, L+1), not reach(704, L+1).
reach(414,L + 1) :- reach(86, L), L <= 40, L>= 37, not reach(83, L+1), not reach(152, L+1), not reach(433, L+1), not reach(704, L+1).
reach(433,L + 1) :- reach(86, L), L <= 40, L>= 37, not reach(83, L+1), not reach(152, L+1), not reach(414, L+1), not reach(704, L+1).
reach(704,L + 1) :- reach(86, L), L <= 40, L>= 37, not reach(83, L+1), not reach(152, L+1), not reach(414, L+1), not reach(433, L+1).
reach(84,L + 1) :- reach(87, L), L <= 36, L>= 34, not reach(88, L+1).
reach(88,L + 1) :- reach(87, L), L <= 36, L>= 34, not reach(84, L+1).
reach(85,L + 1) :- reach(93, L), L <= 38, L>= 38, not reach(95, L+1), not reach(96, L+1), not reach(97, L+1).
reach(95,L + 1) :- reach(93, L), L <= 38, L>= 38, not reach(85, L+1), not reach(96, L+1), not reach(97, L+1).
reach(96,L + 1) :- reach(93, L), L <= 38, L>= 38, not reach(85, L+1), not reach(95, L+1), not reach(97, L+1).
reach(97,L + 1) :- reach(93, L), L <= 38, L>= 38, not reach(85, L+1), not reach(95, L+1), not reach(96, L+1).
reach(85,L + 1) :- reach(94, L), L <= 40, L>= 39, not reach(96, L+1), not reach(390, L+1).
reach(96,L + 1) :- reach(94, L), L <= 40, L>= 39, not reach(85, L+1), not reach(390, L+1).
reach(390,L + 1) :- reach(94, L), L <= 40, L>= 39, not reach(85, L+1), not reach(96, L+1).
reach(86,L + 1) :- reach(152, L), L <= 41, L>= 38, not reach(153, L+1).
reach(153,L + 1) :- reach(152, L), L <= 41, L>= 38, not reach(86, L+1).
reach(86,L + 1) :- reach(414, L), L <= 39, L>= 36, not reach(101, L+1).
reach(101,L + 1) :- reach(414, L), L <= 39, L>= 36, not reach(86, L+1).
reach(86,L + 1) :- reach(433, L), L <= 39, L>= 38, not reach(434, L+1).
reach(434,L + 1) :- reach(433, L), L <= 39, L>= 38, not reach(86, L+1).
reach(86,L + 1) :- reach(704, L), L <= 39, L>= 38, not reach(214, L+1), not reach(708, L+1).
reach(214,L + 1) :- reach(704, L), L <= 39, L>= 38, not reach(86, L+1), not reach(708, L+1).
reach(708,L + 1) :- reach(704, L), L <= 39, L>= 38, not reach(86, L+1), not reach(214, L+1).
reach(88,L + 1) :- reach(100, L), L <= 36, L>= 34, not reach(99, L+1), not reach(103, L+1).
reach(99,L + 1) :- reach(100, L), L <= 36, L>= 34, not reach(88, L+1), not reach(103, L+1).
reach(103,L + 1) :- reach(100, L), L <= 36, L>= 34, not reach(88, L+1), not reach(99, L+1).
reach(88,L + 1) :- reach(106, L), L <= 35, L>= 32, not reach(105, L+1), not reach(108, L+1).
reach(105,L + 1) :- reach(106, L), L <= 35, L>= 32, not reach(88, L+1), not reach(108, L+1).
reach(108,L + 1) :- reach(106, L), L <= 35, L>= 32, not reach(88, L+1), not reach(105, L+1).
reach(88,L + 1) :- reach(110, L), L <= 35, L>= 34.
reach(90,L + 1) :- reach(89, L), L <= 18, L>= 18, not reach(91, L+1).
reach(91,L + 1) :- reach(89, L), L <= 18, L>= 18, not reach(90, L+1).
reach(89,L + 1) :- reach(91, L), L <= 17, L>= 17, not reach(266, L+1), not reach(92, L+1).
reach(266,L + 1) :- reach(91, L), L <= 17, L>= 17, not reach(89, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(91, L), L <= 17, L>= 17, not reach(89, L+1), not reach(266, L+1).
reach(90,L + 1) :- reach(92, L), L <= 18, L>= 16, not reach(91, L+1), not reach(483, L+1), not reach(268, L+1), not reach(718, L+1).
reach(91,L + 1) :- reach(92, L), L <= 18, L>= 16, not reach(90, L+1), not reach(483, L+1), not reach(268, L+1), not reach(718, L+1).
reach(483,L + 1) :- reach(92, L), L <= 18, L>= 16, not reach(90, L+1), not reach(91, L+1), not reach(268, L+1), not reach(718, L+1).
reach(268,L + 1) :- reach(92, L), L <= 18, L>= 16, not reach(90, L+1), not reach(91, L+1), not reach(483, L+1), not reach(718, L+1).
reach(718,L + 1) :- reach(92, L), L <= 18, L>= 16, not reach(90, L+1), not reach(91, L+1), not reach(483, L+1), not reach(268, L+1).
reach(92,L + 1) :- reach(483, L), L <= 17, L>= 15, not reach(484, L+1).
reach(484,L + 1) :- reach(483, L), L <= 17, L>= 15, not reach(92, L+1).
reach(92,L + 1) :- reach(268, L), L <= 18, L>= 17, not reach(266, L+1), not reach(657, L+1).
reach(266,L + 1) :- reach(268, L), L <= 18, L>= 17, not reach(92, L+1), not reach(657, L+1).
reach(657,L + 1) :- reach(268, L), L <= 18, L>= 17, not reach(92, L+1), not reach(266, L+1).
reach(92,L + 1) :- reach(718, L), L <= 17, L>= 17, not reach(379, L+1).
reach(379,L + 1) :- reach(718, L), L <= 17, L>= 17, not reach(92, L+1).
reach(93,L + 1) :- reach(96, L), L <= 39, L>= 39, not reach(94, L+1).
reach(94,L + 1) :- reach(96, L), L <= 39, L>= 39, not reach(93, L+1).
reach(93,L + 1) :- reach(97, L), L <= 38, L>= 37, not reach(387, L+1).
reach(387,L + 1) :- reach(97, L), L <= 38, L>= 37, not reach(93, L+1).
reach(94,L + 1) :- reach(390, L), L <= 41, L>= 40, not reach(395, L+1), not reach(551, L+1).
reach(395,L + 1) :- reach(390, L), L <= 41, L>= 40, not reach(94, L+1), not reach(551, L+1).
reach(551,L + 1) :- reach(390, L), L <= 41, L>= 40, not reach(94, L+1), not reach(395, L+1).
reach(97,L + 1) :- reach(387, L), L <= 39, L>= 36, not reach(388, L+1), not reach(392, L+1).
reach(388,L + 1) :- reach(387, L), L <= 39, L>= 36, not reach(97, L+1), not reach(392, L+1).
reach(392,L + 1) :- reach(387, L), L <= 39, L>= 36, not reach(97, L+1), not reach(388, L+1).
reach(98,L + 1) :- reach(425, L), L <= 28, L>= 27, not reach(426, L+1).
reach(426,L + 1) :- reach(425, L), L <= 28, L>= 27, not reach(98, L+1).
reach(98,L + 1) :- reach(303, L), L <= 28, L>= 25, not reach(191, L+1), not reach(298, L+1).
reach(191,L + 1) :- reach(303, L), L <= 28, L>= 25, not reach(98, L+1), not reach(298, L+1).
reach(298,L + 1) :- reach(303, L), L <= 28, L>= 25, not reach(98, L+1), not reach(191, L+1).
reach(100,L + 1) :- reach(99, L), L <= 37, L>= 35, not reach(101, L+1).
reach(101,L + 1) :- reach(99, L), L <= 37, L>= 35, not reach(100, L+1).
reach(99,L + 1) :- reach(101, L), L <= 38, L>= 35, not reach(372, L+1), not reach(414, L+1).
reach(372,L + 1) :- reach(101, L), L <= 38, L>= 35, not reach(99, L+1), not reach(414, L+1).
reach(414,L + 1) :- reach(101, L), L <= 38, L>= 35, not reach(99, L+1), not reach(372, L+1).
reach(100,L + 1) :- reach(103, L), L <= 35, L>= 34, not reach(102, L+1).
reach(102,L + 1) :- reach(103, L), L <= 35, L>= 34, not reach(100, L+1).
reach(101,L + 1) :- reach(372, L), L <= 37, L>= 34, not reach(373, L+1), not reach(374, L+1).
reach(373,L + 1) :- reach(372, L), L <= 37, L>= 34, not reach(101, L+1), not reach(374, L+1).
reach(374,L + 1) :- reach(372, L), L <= 37, L>= 34, not reach(101, L+1), not reach(373, L+1).
reach(103,L + 1) :- reach(102, L), L <= 34, L>= 33, not reach(104, L+1).
reach(104,L + 1) :- reach(102, L), L <= 34, L>= 33, not reach(103, L+1).
reach(106,L + 1) :- reach(105, L), L <= 34, L>= 31, not reach(107, L+1).
reach(107,L + 1) :- reach(105, L), L <= 34, L>= 31, not reach(106, L+1).
reach(105,L + 1) :- reach(107, L), L <= 33, L>= 30, not reach(307, L+1), not reach(308, L+1).
reach(307,L + 1) :- reach(107, L), L <= 33, L>= 30, not reach(105, L+1), not reach(308, L+1).
reach(308,L + 1) :- reach(107, L), L <= 33, L>= 30, not reach(105, L+1), not reach(307, L+1).
reach(107,L + 1) :- reach(308, L), L <= 32, L>= 29, not reach(310, L+1).
reach(310,L + 1) :- reach(308, L), L <= 32, L>= 29, not reach(107, L+1).
reach(109,L + 1) :- reach(374, L), L <= 36, L>= 33, not reach(372, L+1), not reach(373, L+1).
reach(372,L + 1) :- reach(374, L), L <= 36, L>= 33, not reach(109, L+1), not reach(373, L+1).
reach(373,L + 1) :- reach(374, L), L <= 36, L>= 33, not reach(109, L+1), not reach(372, L+1).
reach(115,L + 1) :- reach(116, L), L <= 58, L>= 55, not reach(117, L+1).
reach(117,L + 1) :- reach(116, L), L <= 58, L>= 55, not reach(115, L+1).
reach(116,L + 1) :- reach(117, L), L <= 57, L>= 54, not reach(253, L+1), not reach(254, L+1).
reach(253,L + 1) :- reach(117, L), L <= 57, L>= 54, not reach(116, L+1), not reach(254, L+1).
reach(254,L + 1) :- reach(117, L), L <= 57, L>= 54, not reach(116, L+1), not reach(253, L+1).
reach(117,L + 1) :- reach(253, L), L <= 56, L>= 55, not reach(498, L+1).
reach(498,L + 1) :- reach(253, L), L <= 56, L>= 55, not reach(117, L+1).
reach(117,L + 1) :- reach(254, L), L <= 56, L>= 53, not reach(671, L+1), not reach(676, L+1).
reach(671,L + 1) :- reach(254, L), L <= 56, L>= 53, not reach(117, L+1), not reach(676, L+1).
reach(676,L + 1) :- reach(254, L), L <= 56, L>= 53, not reach(117, L+1), not reach(671, L+1).
reach(141,L + 1) :- reach(143, L), L <= 28, L>= 27, not reach(176, L+1), not reach(184, L+1), not reach(187, L+1), not reach(189, L+1), not reach(192, L+1).
reach(176,L + 1) :- reach(143, L), L <= 28, L>= 27, not reach(141, L+1), not reach(184, L+1), not reach(187, L+1), not reach(189, L+1), not reach(192, L+1).
reach(184,L + 1) :- reach(143, L), L <= 28, L>= 27, not reach(141, L+1), not reach(176, L+1), not reach(187, L+1), not reach(189, L+1), not reach(192, L+1).
reach(187,L + 1) :- reach(143, L), L <= 28, L>= 27, not reach(141, L+1), not reach(176, L+1), not reach(184, L+1), not reach(189, L+1), not reach(192, L+1).
reach(189,L + 1) :- reach(143, L), L <= 28, L>= 27, not reach(141, L+1), not reach(176, L+1), not reach(184, L+1), not reach(187, L+1), not reach(192, L+1).
reach(192,L + 1) :- reach(143, L), L <= 28, L>= 27, not reach(141, L+1), not reach(176, L+1), not reach(184, L+1), not reach(187, L+1), not reach(189, L+1).
reach(143,L + 1) :- reach(187, L), L <= 27, L>= 27, not reach(185, L+1), not reach(188, L+1).
reach(185,L + 1) :- reach(187, L), L <= 27, L>= 27, not reach(143, L+1), not reach(188, L+1).
reach(188,L + 1) :- reach(187, L), L <= 27, L>= 27, not reach(143, L+1), not reach(185, L+1).
reach(143,L + 1) :- reach(189, L), L <= 27, L>= 26, not reach(190, L+1).
reach(190,L + 1) :- reach(189, L), L <= 27, L>= 26, not reach(143, L+1).
reach(144,L + 1) :- reach(685, L), L <= 58, L>= 56, not reach(377, L+1).
reach(377,L + 1) :- reach(685, L), L <= 58, L>= 56, not reach(144, L+1).
reach(146,L + 1) :- reach(145, L), L <= 49, L>= 46, not reach(147, L+1), not reach(148, L+1).
reach(147,L + 1) :- reach(145, L), L <= 49, L>= 46, not reach(146, L+1), not reach(148, L+1).
reach(148,L + 1) :- reach(145, L), L <= 49, L>= 46, not reach(146, L+1), not reach(147, L+1).
reach(145,L + 1) :- reach(146, L), L <= 48, L>= 46, not reach(276, L+1), not reach(559, L+1).
reach(276,L + 1) :- reach(146, L), L <= 48, L>= 46, not reach(145, L+1), not reach(559, L+1).
reach(559,L + 1) :- reach(146, L), L <= 48, L>= 46, not reach(145, L+1), not reach(276, L+1).
reach(145,L + 1) :- reach(147, L), L <= 48, L>= 45, not reach(150, L+1).
reach(150,L + 1) :- reach(147, L), L <= 48, L>= 45, not reach(145, L+1).
reach(145,L + 1) :- reach(148, L), L <= 50, L>= 47, not reach(660, L+1).
reach(660,L + 1) :- reach(148, L), L <= 50, L>= 47, not reach(145, L+1).
reach(146,L + 1) :- reach(276, L), L <= 47, L>= 45, not reach(277, L+1).
reach(277,L + 1) :- reach(276, L), L <= 47, L>= 45, not reach(146, L+1).
reach(146,L + 1) :- reach(559, L), L <= 47, L>= 47.
reach(147,L + 1) :- reach(150, L), L <= 47, L>= 44, not reach(149, L+1), not reach(560, L+1).
reach(149,L + 1) :- reach(150, L), L <= 47, L>= 44, not reach(147, L+1), not reach(560, L+1).
reach(560,L + 1) :- reach(150, L), L <= 47, L>= 44, not reach(147, L+1), not reach(149, L+1).
reach(148,L + 1) :- reach(660, L), L <= 51, L>= 48, not reach(505, L+1).
reach(505,L + 1) :- reach(660, L), L <= 51, L>= 48, not reach(148, L+1).
reach(150,L + 1) :- reach(149, L), L <= 48, L>= 45, not reach(151, L+1).
reach(151,L + 1) :- reach(149, L), L <= 48, L>= 45, not reach(150, L+1).
reach(149,L + 1) :- reach(151, L), L <= 49, L>= 46, not reach(591, L+1), not reach(593, L+1), not reach(594, L+1).
reach(591,L + 1) :- reach(151, L), L <= 49, L>= 46, not reach(149, L+1), not reach(593, L+1), not reach(594, L+1).
reach(593,L + 1) :- reach(151, L), L <= 49, L>= 46, not reach(149, L+1), not reach(591, L+1), not reach(594, L+1).
reach(594,L + 1) :- reach(151, L), L <= 49, L>= 46, not reach(149, L+1), not reach(591, L+1), not reach(593, L+1).
reach(150,L + 1) :- reach(560, L), L <= 46, L>= 43, not reach(595, L+1), not reach(596, L+1), not reach(597, L+1).
reach(595,L + 1) :- reach(560, L), L <= 46, L>= 43, not reach(150, L+1), not reach(596, L+1), not reach(597, L+1).
reach(596,L + 1) :- reach(560, L), L <= 46, L>= 43, not reach(150, L+1), not reach(595, L+1), not reach(597, L+1).
reach(597,L + 1) :- reach(560, L), L <= 46, L>= 43, not reach(150, L+1), not reach(595, L+1), not reach(596, L+1).
reach(151,L + 1) :- reach(591, L), L <= 50, L>= 47, not reach(497, L+1).
reach(497,L + 1) :- reach(591, L), L <= 50, L>= 47, not reach(151, L+1).
reach(151,L + 1) :- reach(593, L), L <= 48, L>= 45, not reach(595, L+1).
reach(595,L + 1) :- reach(593, L), L <= 48, L>= 45, not reach(151, L+1).
reach(151,L + 1) :- reach(594, L), L <= 48, L>= 45, not reach(518, L+1).
reach(518,L + 1) :- reach(594, L), L <= 48, L>= 45, not reach(151, L+1).
reach(152,L + 1) :- reach(153, L), L <= 42, L>= 39, not reach(154, L+1), not reach(155, L+1), not reach(156, L+1).
reach(154,L + 1) :- reach(153, L), L <= 42, L>= 39, not reach(152, L+1), not reach(155, L+1), not reach(156, L+1).
reach(155,L + 1) :- reach(153, L), L <= 42, L>= 39, not reach(152, L+1), not reach(154, L+1), not reach(156, L+1).
reach(156,L + 1) :- reach(153, L), L <= 42, L>= 39, not reach(152, L+1), not reach(154, L+1), not reach(155, L+1).
reach(153,L + 1) :- reach(154, L), L <= 43, L>= 40, not reach(157, L+1).
reach(157,L + 1) :- reach(154, L), L <= 43, L>= 40, not reach(153, L+1).
reach(153,L + 1) :- reach(155, L), L <= 41, L>= 40, not reach(158, L+1), not reach(161, L+1).
reach(158,L + 1) :- reach(155, L), L <= 41, L>= 40, not reach(153, L+1), not reach(161, L+1).
reach(161,L + 1) :- reach(155, L), L <= 41, L>= 40, not reach(153, L+1), not reach(158, L+1).
reach(153,L + 1) :- reach(156, L), L <= 41, L>= 40, not reach(160, L+1), not reach(440, L+1).
reach(160,L + 1) :- reach(156, L), L <= 41, L>= 40, not reach(153, L+1), not reach(440, L+1).
reach(440,L + 1) :- reach(156, L), L <= 41, L>= 40, not reach(153, L+1), not reach(160, L+1).
reach(154,L + 1) :- reach(157, L), L <= 44, L>= 41, not reach(552, L+1), not reach(737, L+1).
reach(552,L + 1) :- reach(157, L), L <= 44, L>= 41, not reach(154, L+1), not reach(737, L+1).
reach(737,L + 1) :- reach(157, L), L <= 44, L>= 41, not reach(154, L+1), not reach(552, L+1).
reach(156,L + 1) :- reach(440, L), L <= 40, L>= 40, not reach(434, L+1).
reach(434,L + 1) :- reach(440, L), L <= 40, L>= 40, not reach(156, L+1).
reach(157,L + 1) :- reach(552, L), L <= 43, L>= 41, not reach(397, L+1), not reach(726, L+1).
reach(397,L + 1) :- reach(552, L), L <= 43, L>= 41, not reach(157, L+1), not reach(726, L+1).
reach(726,L + 1) :- reach(552, L), L <= 43, L>= 41, not reach(157, L+1), not reach(397, L+1).
reach(157,L + 1) :- reach(737, L), L <= 45, L>= 42, not reach(727, L+1), not reach(731, L+1), not reach(738, L+1).
reach(727,L + 1) :- reach(737, L), L <= 45, L>= 42, not reach(157, L+1), not reach(731, L+1), not reach(738, L+1).
reach(731,L + 1) :- reach(737, L), L <= 45, L>= 42, not reach(157, L+1), not reach(727, L+1), not reach(738, L+1).
reach(738,L + 1) :- reach(737, L), L <= 45, L>= 42, not reach(157, L+1), not reach(727, L+1), not reach(731, L+1).
reach(167,L + 1) :- reach(166, L), L <= 46, L>= 44, not reach(168, L+1).
reach(168,L + 1) :- reach(166, L), L <= 46, L>= 44, not reach(167, L+1).
reach(166,L + 1) :- reach(167, L), L <= 47, L>= 45, not reach(384, L+1).
reach(384,L + 1) :- reach(167, L), L <= 47, L>= 45, not reach(166, L+1).
reach(166,L + 1) :- reach(168, L), L <= 45, L>= 43, not reach(678, L+1).
reach(678,L + 1) :- reach(168, L), L <= 45, L>= 43, not reach(166, L+1).
reach(167,L + 1) :- reach(384, L), L <= 48, L>= 46, not reach(382, L+1).
reach(382,L + 1) :- reach(384, L), L <= 48, L>= 46, not reach(167, L+1).
reach(168,L + 1) :- reach(678, L), L <= 44, L>= 42, not reach(563, L+1), not reach(705, L+1).
reach(563,L + 1) :- reach(678, L), L <= 44, L>= 42, not reach(168, L+1), not reach(705, L+1).
reach(705,L + 1) :- reach(678, L), L <= 44, L>= 42, not reach(168, L+1), not reach(563, L+1).
reach(187,L + 1) :- reach(188, L), L <= 26, L>= 26, not reach(369, L+1).
reach(369,L + 1) :- reach(188, L), L <= 26, L>= 26, not reach(187, L+1).
reach(188,L + 1) :- reach(369, L), L <= 25, L>= 25, not reach(365, L+1), not reach(367, L+1), not reach(371, L+1).
reach(365,L + 1) :- reach(369, L), L <= 25, L>= 25, not reach(188, L+1), not reach(367, L+1), not reach(371, L+1).
reach(367,L + 1) :- reach(369, L), L <= 25, L>= 25, not reach(188, L+1), not reach(365, L+1), not reach(371, L+1).
reach(371,L + 1) :- reach(369, L), L <= 25, L>= 25, not reach(188, L+1), not reach(365, L+1), not reach(367, L+1).
reach(189,L + 1) :- reach(190, L), L <= 26, L>= 25, not reach(191, L+1).
reach(191,L + 1) :- reach(190, L), L <= 26, L>= 25, not reach(189, L+1).
reach(190,L + 1) :- reach(191, L), L <= 27, L>= 24, not reach(365, L+1), not reach(364, L+1), not reach(303, L+1).
reach(365,L + 1) :- reach(191, L), L <= 27, L>= 24, not reach(190, L+1), not reach(364, L+1), not reach(303, L+1).
reach(364,L + 1) :- reach(191, L), L <= 27, L>= 24, not reach(190, L+1), not reach(365, L+1), not reach(303, L+1).
reach(303,L + 1) :- reach(191, L), L <= 27, L>= 24, not reach(190, L+1), not reach(365, L+1), not reach(364, L+1).
reach(191,L + 1) :- reach(365, L), L <= 26, L>= 25, not reach(369, L+1), not reach(371, L+1).
reach(369,L + 1) :- reach(365, L), L <= 26, L>= 25, not reach(191, L+1), not reach(371, L+1).
reach(371,L + 1) :- reach(365, L), L <= 26, L>= 25, not reach(191, L+1), not reach(369, L+1).
reach(191,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(301, L+1), not reach(366, L+1), not reach(367, L+1), not reach(409, L+1), not reach(570, L+1), not reach(571, L+1).
reach(301,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(191, L+1), not reach(366, L+1), not reach(367, L+1), not reach(409, L+1), not reach(570, L+1), not reach(571, L+1).
reach(366,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(191, L+1), not reach(301, L+1), not reach(367, L+1), not reach(409, L+1), not reach(570, L+1), not reach(571, L+1).
reach(367,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(191, L+1), not reach(301, L+1), not reach(366, L+1), not reach(409, L+1), not reach(570, L+1), not reach(571, L+1).
reach(409,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(191, L+1), not reach(301, L+1), not reach(366, L+1), not reach(367, L+1), not reach(570, L+1), not reach(571, L+1).
reach(570,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(191, L+1), not reach(301, L+1), not reach(366, L+1), not reach(367, L+1), not reach(409, L+1), not reach(571, L+1).
reach(571,L + 1) :- reach(364, L), L <= 26, L>= 23, not reach(191, L+1), not reach(301, L+1), not reach(366, L+1), not reach(367, L+1), not reach(409, L+1), not reach(570, L+1).
reach(194,L + 1) :- reach(193, L), L <= 14, L>= 13, not reach(195, L+1).
reach(195,L + 1) :- reach(193, L), L <= 14, L>= 13, not reach(194, L+1).
reach(193,L + 1) :- reach(194, L), L <= 15, L>= 14, not reach(200, L+1).
reach(200,L + 1) :- reach(194, L), L <= 15, L>= 14, not reach(193, L+1).
reach(193,L + 1) :- reach(195, L), L <= 13, L>= 12, not reach(202, L+1), not reach(204, L+1), not reach(205, L+1).
reach(202,L + 1) :- reach(195, L), L <= 13, L>= 12, not reach(193, L+1), not reach(204, L+1), not reach(205, L+1).
reach(204,L + 1) :- reach(195, L), L <= 13, L>= 12, not reach(193, L+1), not reach(202, L+1), not reach(205, L+1).
reach(205,L + 1) :- reach(195, L), L <= 13, L>= 12, not reach(193, L+1), not reach(202, L+1), not reach(204, L+1).
reach(195,L + 1) :- reach(205, L), L <= 12, L>= 11, not reach(456, L+1).
reach(456,L + 1) :- reach(205, L), L <= 12, L>= 11, not reach(195, L+1).
reach(205,L + 1) :- reach(456, L), L <= 11, L>= 10, not reach(283, L+1), not reach(457, L+1).
reach(283,L + 1) :- reach(456, L), L <= 11, L>= 10, not reach(205, L+1), not reach(457, L+1).
reach(457,L + 1) :- reach(456, L), L <= 11, L>= 10, not reach(205, L+1), not reach(283, L+1).
reach(208,L + 1) :- reach(207, L), L <= 21, L>= 19, not reach(209, L+1).
reach(209,L + 1) :- reach(207, L), L <= 21, L>= 19, not reach(208, L+1).
reach(207,L + 1) :- reach(209, L), L <= 22, L>= 20, not reach(492, L+1), not reach(445, L+1), not reach(667, L+1), not reach(740, L+1).
reach(492,L + 1) :- reach(209, L), L <= 22, L>= 20, not reach(207, L+1), not reach(445, L+1), not reach(667, L+1), not reach(740, L+1).
reach(445,L + 1) :- reach(209, L), L <= 22, L>= 20, not reach(207, L+1), not reach(492, L+1), not reach(667, L+1), not reach(740, L+1).
reach(667,L + 1) :- reach(209, L), L <= 22, L>= 20, not reach(207, L+1), not reach(492, L+1), not reach(445, L+1), not reach(740, L+1).
reach(740,L + 1) :- reach(209, L), L <= 22, L>= 20, not reach(207, L+1), not reach(492, L+1), not reach(445, L+1), not reach(667, L+1).
reach(209,L + 1) :- reach(492, L), L <= 23, L>= 21, not reach(493, L+1).
reach(493,L + 1) :- reach(492, L), L <= 23, L>= 21, not reach(209, L+1).
reach(209,L + 1) :- reach(445, L), L <= 23, L>= 21, not reach(444, L+1), not reach(669, L+1).
reach(444,L + 1) :- reach(445, L), L <= 23, L>= 21, not reach(209, L+1), not reach(669, L+1).
reach(669,L + 1) :- reach(445, L), L <= 23, L>= 21, not reach(209, L+1), not reach(444, L+1).
reach(209,L + 1) :- reach(667, L), L <= 21, L>= 21, not reach(651, L+1).
reach(651,L + 1) :- reach(667, L), L <= 21, L>= 21, not reach(209, L+1).
reach(209,L + 1) :- reach(740, L), L <= 22, L>= 21, not reach(663, L+1).
reach(663,L + 1) :- reach(740, L), L <= 22, L>= 21, not reach(209, L+1).
reach(228,L + 1) :- reach(681, L), L <= 46, L>= 45, not reach(518, L+1), not reach(679, L+1).
reach(518,L + 1) :- reach(681, L), L <= 46, L>= 45, not reach(228, L+1), not reach(679, L+1).
reach(679,L + 1) :- reach(681, L), L <= 46, L>= 45, not reach(228, L+1), not reach(518, L+1).
reach(253,L + 1) :- reach(498, L), L <= 55, L>= 55, not reach(499, L+1).
reach(499,L + 1) :- reach(498, L), L <= 55, L>= 55, not reach(253, L+1).
reach(254,L + 1) :- reach(671, L), L <= 55, L>= 52, not reach(672, L+1).
reach(672,L + 1) :- reach(671, L), L <= 55, L>= 52, not reach(254, L+1).
reach(254,L + 1) :- reach(676, L), L <= 55, L>= 52, not reach(674, L+1).
reach(674,L + 1) :- reach(676, L), L <= 55, L>= 52, not reach(254, L+1).
reach(267,L + 1) :- reach(656, L), L <= 18, L>= 18, not reach(657, L+1).
reach(657,L + 1) :- reach(656, L), L <= 18, L>= 18, not reach(267, L+1).
reach(272,L + 1) :- reach(271, L), L <= 27, L>= 25, not reach(273, L+1), not reach(274, L+1), not reach(275, L+1).
reach(273,L + 1) :- reach(271, L), L <= 27, L>= 25, not reach(272, L+1), not reach(274, L+1), not reach(275, L+1).
reach(274,L + 1) :- reach(271, L), L <= 27, L>= 25, not reach(272, L+1), not reach(273, L+1), not reach(275, L+1).
reach(275,L + 1) :- reach(271, L), L <= 27, L>= 25, not reach(272, L+1), not reach(273, L+1), not reach(274, L+1).
reach(271,L + 1) :- reach(273, L), L <= 26, L>= 26, not reach(532, L+1), not reach(533, L+1).
reach(532,L + 1) :- reach(273, L), L <= 26, L>= 26, not reach(271, L+1), not reach(533, L+1).
reach(533,L + 1) :- reach(273, L), L <= 26, L>= 26, not reach(271, L+1), not reach(532, L+1).
reach(271,L + 1) :- reach(274, L), L <= 26, L>= 24, not reach(438, L+1), not reach(670, L+1).
reach(438,L + 1) :- reach(274, L), L <= 26, L>= 24, not reach(271, L+1), not reach(670, L+1).
reach(670,L + 1) :- reach(274, L), L <= 26, L>= 24, not reach(271, L+1), not reach(438, L+1).
reach(274,L + 1) :- reach(438, L), L <= 25, L>= 25, not reach(435, L+1), not reach(719, L+1).
reach(435,L + 1) :- reach(438, L), L <= 25, L>= 25, not reach(274, L+1), not reach(719, L+1).
reach(719,L + 1) :- reach(438, L), L <= 25, L>= 25, not reach(274, L+1), not reach(435, L+1).
reach(274,L + 1) :- reach(670, L), L <= 25, L>= 23, not reach(669, L+1), not reach(728, L+1).
reach(669,L + 1) :- reach(670, L), L <= 25, L>= 23, not reach(274, L+1), not reach(728, L+1).
reach(728,L + 1) :- reach(670, L), L <= 25, L>= 23, not reach(274, L+1), not reach(669, L+1).
reach(276,L + 1) :- reach(277, L), L <= 46, L>= 44, not reach(698, L+1).
reach(698,L + 1) :- reach(277, L), L <= 46, L>= 44, not reach(276, L+1).
reach(277,L + 1) :- reach(698, L), L <= 45, L>= 43, not reach(677, L+1), not reach(697, L+1), not reach(699, L+1), not reach(705, L+1).
reach(677,L + 1) :- reach(698, L), L <= 45, L>= 43, not reach(277, L+1), not reach(697, L+1), not reach(699, L+1), not reach(705, L+1).
reach(697,L + 1) :- reach(698, L), L <= 45, L>= 43, not reach(277, L+1), not reach(677, L+1), not reach(699, L+1), not reach(705, L+1).
reach(699,L + 1) :- reach(698, L), L <= 45, L>= 43, not reach(277, L+1), not reach(677, L+1), not reach(697, L+1), not reach(705, L+1).
reach(705,L + 1) :- reach(698, L), L <= 45, L>= 43, not reach(277, L+1), not reach(677, L+1), not reach(697, L+1), not reach(699, L+1).
reach(279,L + 1) :- reach(427, L), L <= 13, L>= 12, not reach(337, L+1), not reach(487, L+1).
reach(337,L + 1) :- reach(427, L), L <= 13, L>= 12, not reach(279, L+1), not reach(487, L+1).
reach(487,L + 1) :- reach(427, L), L <= 13, L>= 12, not reach(279, L+1), not reach(337, L+1).
reach(280,L + 1) :- reach(428, L), L <= 13, L>= 13, not reach(487, L+1).
reach(487,L + 1) :- reach(428, L), L <= 13, L>= 13, not reach(280, L+1).
reach(281,L + 1) :- reach(490, L), L <= 14, L>= 14, not reach(489, L+1), not reach(494, L+1).
reach(489,L + 1) :- reach(490, L), L <= 14, L>= 14, not reach(281, L+1), not reach(494, L+1).
reach(494,L + 1) :- reach(490, L), L <= 14, L>= 14, not reach(281, L+1), not reach(489, L+1).
reach(282,L + 1) :- reach(283, L), L <= 10, L>= 10, not reach(456, L+1), not reach(457, L+1).
reach(456,L + 1) :- reach(283, L), L <= 10, L>= 10, not reach(282, L+1), not reach(457, L+1).
reach(457,L + 1) :- reach(283, L), L <= 10, L>= 10, not reach(282, L+1), not reach(456, L+1).
reach(283,L + 1) :- reach(457, L), L <= 10, L>= 9, not reach(456, L+1), not reach(586, L+1), not reach(587, L+1), not reach(589, L+1).
reach(456,L + 1) :- reach(457, L), L <= 10, L>= 9, not reach(283, L+1), not reach(586, L+1), not reach(587, L+1), not reach(589, L+1).
reach(586,L + 1) :- reach(457, L), L <= 10, L>= 9, not reach(283, L+1), not reach(456, L+1), not reach(587, L+1), not reach(589, L+1).
reach(587,L + 1) :- reach(457, L), L <= 10, L>= 9, not reach(283, L+1), not reach(456, L+1), not reach(586, L+1), not reach(589, L+1).
reach(589,L + 1) :- reach(457, L), L <= 10, L>= 9, not reach(283, L+1), not reach(456, L+1), not reach(586, L+1), not reach(587, L+1).
reach(298,L + 1) :- reach(297, L), L <= 29, L>= 26, not reach(299, L+1), not reach(300, L+1).
reach(299,L + 1) :- reach(297, L), L <= 29, L>= 26, not reach(298, L+1), not reach(300, L+1).
reach(300,L + 1) :- reach(297, L), L <= 29, L>= 26, not reach(298, L+1), not reach(299, L+1).
reach(297,L + 1) :- reach(298, L), L <= 28, L>= 25, not reach(301, L+1), not reach(302, L+1), not reach(303, L+1).
reach(301,L + 1) :- reach(298, L), L <= 28, L>= 25, not reach(297, L+1), not reach(302, L+1), not reach(303, L+1).
reach(302,L + 1) :- reach(298, L), L <= 28, L>= 25, not reach(297, L+1), not reach(301, L+1), not reach(303, L+1).
reach(303,L + 1) :- reach(298, L), L <= 28, L>= 25, not reach(297, L+1), not reach(301, L+1), not reach(302, L+1).
reach(297,L + 1) :- reach(299, L), L <= 28, L>= 27.
reach(297,L + 1) :- reach(300, L), L <= 30, L>= 27, not reach(310, L+1), not reach(655, L+1).
reach(310,L + 1) :- reach(300, L), L <= 30, L>= 27, not reach(297, L+1), not reach(655, L+1).
reach(655,L + 1) :- reach(300, L), L <= 30, L>= 27, not reach(297, L+1), not reach(310, L+1).
reach(298,L + 1) :- reach(301, L), L <= 27, L>= 24, not reach(363, L+1), not reach(364, L+1).
reach(363,L + 1) :- reach(301, L), L <= 27, L>= 24, not reach(298, L+1), not reach(364, L+1).
reach(364,L + 1) :- reach(301, L), L <= 27, L>= 24, not reach(298, L+1), not reach(363, L+1).
reach(298,L + 1) :- reach(302, L), L <= 27, L>= 26, not reach(652, L+1), not reach(653, L+1).
reach(652,L + 1) :- reach(302, L), L <= 27, L>= 26, not reach(298, L+1), not reach(653, L+1).
reach(653,L + 1) :- reach(302, L), L <= 27, L>= 26, not reach(298, L+1), not reach(652, L+1).
reach(300,L + 1) :- reach(310, L), L <= 31, L>= 28, not reach(308, L+1), not reach(309, L+1).
reach(308,L + 1) :- reach(310, L), L <= 31, L>= 28, not reach(300, L+1), not reach(309, L+1).
reach(309,L + 1) :- reach(310, L), L <= 31, L>= 28, not reach(300, L+1), not reach(308, L+1).
reach(300,L + 1) :- reach(655, L), L <= 29, L>= 28, not reach(654, L+1).
reach(654,L + 1) :- reach(655, L), L <= 29, L>= 28, not reach(300, L+1).
reach(301,L + 1) :- reach(363, L), L <= 26, L>= 25, not reach(366, L+1).
reach(366,L + 1) :- reach(363, L), L <= 26, L>= 25, not reach(301, L+1).
reach(302,L + 1) :- reach(652, L), L <= 27, L>= 27, not reach(654, L+1).
reach(654,L + 1) :- reach(652, L), L <= 27, L>= 27, not reach(302, L+1).
reach(307,L + 1) :- reach(309, L), L <= 31, L>= 29, not reach(310, L+1).
reach(310,L + 1) :- reach(309, L), L <= 31, L>= 29, not reach(307, L+1).
reach(312,L + 1) :- reach(728, L), L <= 24, L>= 24, not reach(670, L+1).
reach(670,L + 1) :- reach(728, L), L <= 24, L>= 24, not reach(312, L+1).
reach(318,L + 1) :- reach(370, L), L <= 23, L>= 23, not reach(368, L+1), not reach(650, L+1), not reach(747, L+1), not reach(748, L+1).
reach(368,L + 1) :- reach(370, L), L <= 23, L>= 23, not reach(318, L+1), not reach(650, L+1), not reach(747, L+1), not reach(748, L+1).
reach(650,L + 1) :- reach(370, L), L <= 23, L>= 23, not reach(318, L+1), not reach(368, L+1), not reach(747, L+1), not reach(748, L+1).
reach(747,L + 1) :- reach(370, L), L <= 23, L>= 23, not reach(318, L+1), not reach(368, L+1), not reach(650, L+1), not reach(748, L+1).
reach(748,L + 1) :- reach(370, L), L <= 23, L>= 23, not reach(318, L+1), not reach(368, L+1), not reach(650, L+1), not reach(747, L+1).
reach(336,L + 1) :- reach(335, L), L <= 11, L>= 10, not reach(337, L+1), not reach(338, L+1).
reach(337,L + 1) :- reach(335, L), L <= 11, L>= 10, not reach(336, L+1), not reach(338, L+1).
reach(338,L + 1) :- reach(335, L), L <= 11, L>= 10, not reach(336, L+1), not reach(337, L+1).
reach(335,L + 1) :- reach(336, L), L <= 10, L>= 9, not reach(339, L+1).
reach(339,L + 1) :- reach(336, L), L <= 10, L>= 9, not reach(335, L+1).
reach(335,L + 1) :- reach(337, L), L <= 12, L>= 11, not reach(427, L+1).
reach(427,L + 1) :- reach(337, L), L <= 12, L>= 11, not reach(335, L+1).
reach(350,L + 1) :- reach(406, L), L <= 21, L>= 20, not reach(404, L+1).
reach(404,L + 1) :- reach(406, L), L <= 21, L>= 20, not reach(350, L+1).
reach(360,L + 1) :- reach(359, L), L <= 18, L>= 15, not reach(361, L+1), not reach(362, L+1).
reach(361,L + 1) :- reach(359, L), L <= 18, L>= 15, not reach(360, L+1), not reach(362, L+1).
reach(362,L + 1) :- reach(359, L), L <= 18, L>= 15, not reach(360, L+1), not reach(361, L+1).
reach(359,L + 1) :- reach(361, L), L <= 17, L>= 14, not reach(658, L+1).
reach(658,L + 1) :- reach(361, L), L <= 17, L>= 14, not reach(359, L+1).
reach(359,L + 1) :- reach(362, L), L <= 18, L>= 16, not reach(657, L+1).
reach(657,L + 1) :- reach(362, L), L <= 18, L>= 16, not reach(359, L+1).
reach(361,L + 1) :- reach(658, L), L <= 16, L>= 13, not reach(482, L+1).
reach(482,L + 1) :- reach(658, L), L <= 16, L>= 13, not reach(361, L+1).
reach(363,L + 1) :- reach(366, L), L <= 25, L>= 24, not reach(364, L+1).
reach(364,L + 1) :- reach(366, L), L <= 25, L>= 24, not reach(363, L+1).
reach(364,L + 1) :- reach(367, L), L <= 25, L>= 24, not reach(368, L+1), not reach(369, L+1).
reach(368,L + 1) :- reach(367, L), L <= 25, L>= 24, not reach(364, L+1), not reach(369, L+1).
reach(369,L + 1) :- reach(367, L), L <= 25, L>= 24, not reach(364, L+1), not reach(368, L+1).
reach(364,L + 1) :- reach(409, L), L <= 25, L>= 24, not reach(410, L+1).
reach(410,L + 1) :- reach(409, L), L <= 25, L>= 24, not reach(364, L+1).
reach(364,L + 1) :- reach(570, L), L <= 25, L>= 22, not reach(574, L+1).
reach(574,L + 1) :- reach(570, L), L <= 25, L>= 22, not reach(364, L+1).
reach(364,L + 1) :- reach(571, L), L <= 25, L>= 24.
reach(367,L + 1) :- reach(368, L), L <= 24, L>= 24, not reach(370, L+1).
reach(370,L + 1) :- reach(368, L), L <= 24, L>= 24, not reach(367, L+1).
reach(370,L + 1) :- reach(650, L), L <= 22, L>= 22, not reach(651, L+1).
reach(651,L + 1) :- reach(650, L), L <= 22, L>= 22, not reach(370, L+1).
reach(372,L + 1) :- reach(373, L), L <= 36, L>= 34, not reach(374, L+1).
reach(374,L + 1) :- reach(373, L), L <= 36, L>= 34, not reach(372, L+1).
reach(376,L + 1) :- reach(375, L), L <= 56, L>= 54, not reach(377, L+1).
reach(377,L + 1) :- reach(375, L), L <= 56, L>= 54, not reach(376, L+1).
reach(375,L + 1) :- reach(376, L), L <= 55, L>= 53, not reach(378, L+1).
reach(378,L + 1) :- reach(376, L), L <= 55, L>= 53, not reach(375, L+1).
reach(375,L + 1) :- reach(377, L), L <= 57, L>= 55, not reach(685, L+1).
reach(685,L + 1) :- reach(377, L), L <= 57, L>= 55, not reach(375, L+1).
reach(376,L + 1) :- reach(378, L), L <= 54, L>= 52, not reach(386, L+1).
reach(386,L + 1) :- reach(378, L), L <= 54, L>= 52, not reach(376, L+1).
reach(378,L + 1) :- reach(386, L), L <= 53, L>= 51, not reach(383, L+1).
reach(383,L + 1) :- reach(386, L), L <= 53, L>= 51, not reach(378, L+1).
reach(381,L + 1) :- reach(380, L), L <= 50, L>= 48, not reach(382, L+1).
reach(382,L + 1) :- reach(380, L), L <= 50, L>= 48, not reach(381, L+1).
reach(380,L + 1) :- reach(381, L), L <= 51, L>= 49, not reach(383, L+1).
reach(383,L + 1) :- reach(381, L), L <= 51, L>= 49, not reach(380, L+1).
reach(380,L + 1) :- reach(382, L), L <= 49, L>= 47, not reach(384, L+1).
reach(384,L + 1) :- reach(382, L), L <= 49, L>= 47, not reach(380, L+1).
reach(381,L + 1) :- reach(383, L), L <= 52, L>= 50, not reach(385, L+1), not reach(386, L+1), not reach(503, L+1).
reach(385,L + 1) :- reach(383, L), L <= 52, L>= 50, not reach(381, L+1), not reach(386, L+1), not reach(503, L+1).
reach(386,L + 1) :- reach(383, L), L <= 52, L>= 50, not reach(381, L+1), not reach(385, L+1), not reach(503, L+1).
reach(503,L + 1) :- reach(383, L), L <= 52, L>= 50, not reach(381, L+1), not reach(385, L+1), not reach(386, L+1).
reach(383,L + 1) :- reach(503, L), L <= 51, L>= 51, not reach(501, L+1).
reach(501,L + 1) :- reach(503, L), L <= 51, L>= 51, not reach(383, L+1).
reach(387,L + 1) :- reach(388, L), L <= 40, L>= 37, not reach(389, L+1).
reach(389,L + 1) :- reach(388, L), L <= 40, L>= 37, not reach(387, L+1).
reach(388,L + 1) :- reach(389, L), L <= 41, L>= 38, not reach(395, L+1).
reach(395,L + 1) :- reach(389, L), L <= 41, L>= 38, not reach(388, L+1).
reach(389,L + 1) :- reach(395, L), L <= 42, L>= 39, not reach(390, L+1), not reach(396, L+1), not reach(397, L+1), not reach(398, L+1), not reach(399, L+1).
reach(390,L + 1) :- reach(395, L), L <= 42, L>= 39, not reach(389, L+1), not reach(396, L+1), not reach(397, L+1), not reach(398, L+1), not reach(399, L+1).
reach(396,L + 1) :- reach(395, L), L <= 42, L>= 39, not reach(389, L+1), not reach(390, L+1), not reach(397, L+1), not reach(398, L+1), not reach(399, L+1).
reach(397,L + 1) :- reach(395, L), L <= 42, L>= 39, not reach(389, L+1), not reach(390, L+1), not reach(396, L+1), not reach(398, L+1), not reach(399, L+1).
reach(398,L + 1) :- reach(395, L), L <= 42, L>= 39, not reach(389, L+1), not reach(390, L+1), not reach(396, L+1), not reach(397, L+1), not reach(399, L+1).
reach(399,L + 1) :- reach(395, L), L <= 42, L>= 39, not reach(389, L+1), not reach(390, L+1), not reach(396, L+1), not reach(397, L+1), not reach(398, L+1).
reach(390,L + 1) :- reach(551, L), L <= 41, L>= 41, not reach(397, L+1).
reach(397,L + 1) :- reach(551, L), L <= 41, L>= 41, not reach(390, L+1).
reach(392,L + 1) :- reach(393, L), L <= 37, L>= 36, not reach(394, L+1).
reach(394,L + 1) :- reach(393, L), L <= 37, L>= 36, not reach(392, L+1).
reach(395,L + 1) :- reach(396, L), L <= 41, L>= 40, not reach(469, L+1).
reach(469,L + 1) :- reach(396, L), L <= 41, L>= 40, not reach(395, L+1).
reach(395,L + 1) :- reach(397, L), L <= 42, L>= 40, not reach(551, L+1), not reach(552, L+1).
reach(551,L + 1) :- reach(397, L), L <= 42, L>= 40, not reach(395, L+1), not reach(552, L+1).
reach(552,L + 1) :- reach(397, L), L <= 42, L>= 40, not reach(395, L+1), not reach(551, L+1).
reach(395,L + 1) :- reach(398, L), L <= 41, L>= 40, not reach(720, L+1).
reach(720,L + 1) :- reach(398, L), L <= 41, L>= 40, not reach(395, L+1).
reach(395,L + 1) :- reach(399, L), L <= 43, L>= 40, not reach(723, L+1).
reach(723,L + 1) :- reach(399, L), L <= 43, L>= 40, not reach(395, L+1).
reach(398,L + 1) :- reach(720, L), L <= 42, L>= 41, not reach(700, L+1).
reach(700,L + 1) :- reach(720, L), L <= 42, L>= 41, not reach(398, L+1).
reach(399,L + 1) :- reach(723, L), L <= 44, L>= 41, not reach(597, L+1).
reach(597,L + 1) :- reach(723, L), L <= 44, L>= 41, not reach(399, L+1).
reach(405,L + 1) :- reach(404, L), L <= 22, L>= 19, not reach(406, L+1), not reach(407, L+1).
reach(406,L + 1) :- reach(404, L), L <= 22, L>= 19, not reach(405, L+1), not reach(407, L+1).
reach(407,L + 1) :- reach(404, L), L <= 22, L>= 19, not reach(405, L+1), not reach(406, L+1).
reach(404,L + 1) :- reach(407, L), L <= 23, L>= 20, not reach(574, L+1), not reach(651, L+1).
reach(574,L + 1) :- reach(407, L), L <= 23, L>= 20, not reach(404, L+1), not reach(651, L+1).
reach(651,L + 1) :- reach(407, L), L <= 23, L>= 20, not reach(404, L+1), not reach(574, L+1).
reach(407,L + 1) :- reach(574, L), L <= 24, L>= 21, not reach(570, L+1), not reach(648, L+1).
reach(570,L + 1) :- reach(574, L), L <= 24, L>= 21, not reach(407, L+1), not reach(648, L+1).
reach(648,L + 1) :- reach(574, L), L <= 24, L>= 21, not reach(407, L+1), not reach(570, L+1).
reach(407,L + 1) :- reach(651, L), L <= 22, L>= 21, not reach(576, L+1), not reach(650, L+1), not reach(667, L+1).
reach(576,L + 1) :- reach(651, L), L <= 22, L>= 21, not reach(407, L+1), not reach(650, L+1), not reach(667, L+1).
reach(650,L + 1) :- reach(651, L), L <= 22, L>= 21, not reach(407, L+1), not reach(576, L+1), not reach(667, L+1).
reach(667,L + 1) :- reach(651, L), L <= 22, L>= 21, not reach(407, L+1), not reach(576, L+1), not reach(650, L+1).
reach(433,L + 1) :- reach(434, L), L <= 39, L>= 39, not reach(439, L+1), not reach(440, L+1).
reach(439,L + 1) :- reach(434, L), L <= 39, L>= 39, not reach(433, L+1), not reach(440, L+1).
reach(440,L + 1) :- reach(434, L), L <= 39, L>= 39, not reach(433, L+1), not reach(439, L+1).
reach(445,L + 1) :- reach(444, L), L <= 22, L>= 22, not reach(446, L+1).
reach(446,L + 1) :- reach(444, L), L <= 22, L>= 22, not reach(445, L+1).
reach(445,L + 1) :- reach(669, L), L <= 24, L>= 22, not reach(670, L+1).
reach(670,L + 1) :- reach(669, L), L <= 24, L>= 22, not reach(445, L+1).
reach(457,L + 1) :- reach(586, L), L <= 9, L>= 8, not reach(488, L+1).
reach(488,L + 1) :- reach(586, L), L <= 9, L>= 8, not reach(457, L+1).
reach(457,L + 1) :- reach(587, L), L <= 9, L>= 8, not reach(488, L+1), not reach(588, L+1).
reach(488,L + 1) :- reach(587, L), L <= 9, L>= 8, not reach(457, L+1), not reach(588, L+1).
reach(588,L + 1) :- reach(587, L), L <= 9, L>= 8, not reach(457, L+1), not reach(488, L+1).
reach(483,L + 1) :- reach(484, L), L <= 16, L>= 14, not reach(489, L+1).
reach(489,L + 1) :- reach(484, L), L <= 16, L>= 14, not reach(483, L+1).
reach(484,L + 1) :- reach(489, L), L <= 15, L>= 13, not reach(487, L+1), not reach(490, L+1).
reach(487,L + 1) :- reach(489, L), L <= 15, L>= 13, not reach(484, L+1), not reach(490, L+1).
reach(490,L + 1) :- reach(489, L), L <= 15, L>= 13, not reach(484, L+1), not reach(487, L+1).
reach(485,L + 1) :- reach(486, L), L <= 7, L>= 6, not reach(488, L+1).
reach(488,L + 1) :- reach(486, L), L <= 7, L>= 6, not reach(485, L+1).
reach(486,L + 1) :- reach(488, L), L <= 8, L>= 7, not reach(586, L+1), not reach(587, L+1), not reach(583, L+1).
reach(586,L + 1) :- reach(488, L), L <= 8, L>= 7, not reach(486, L+1), not reach(587, L+1), not reach(583, L+1).
reach(587,L + 1) :- reach(488, L), L <= 8, L>= 7, not reach(486, L+1), not reach(586, L+1), not reach(583, L+1).
reach(583,L + 1) :- reach(488, L), L <= 8, L>= 7, not reach(486, L+1), not reach(586, L+1), not reach(587, L+1).
reach(492,L + 1) :- reach(493, L), L <= 24, L>= 22, not reach(666, L+1), not reach(668, L+1), not reach(663, L+1).
reach(666,L + 1) :- reach(493, L), L <= 24, L>= 22, not reach(492, L+1), not reach(668, L+1), not reach(663, L+1).
reach(668,L + 1) :- reach(493, L), L <= 24, L>= 22, not reach(492, L+1), not reach(666, L+1), not reach(663, L+1).
reach(663,L + 1) :- reach(493, L), L <= 24, L>= 22, not reach(492, L+1), not reach(666, L+1), not reach(668, L+1).
reach(493,L + 1) :- reach(668, L), L <= 23, L>= 23, not reach(665, L+1).
reach(665,L + 1) :- reach(668, L), L <= 23, L>= 23, not reach(493, L+1).
reach(496,L + 1) :- reach(495, L), L <= 50, L>= 49, not reach(497, L+1).
reach(497,L + 1) :- reach(495, L), L <= 50, L>= 49, not reach(496, L+1).
reach(495,L + 1) :- reach(497, L), L <= 51, L>= 48, not reach(590, L+1), not reach(591, L+1).
reach(590,L + 1) :- reach(497, L), L <= 51, L>= 48, not reach(495, L+1), not reach(591, L+1).
reach(591,L + 1) :- reach(497, L), L <= 51, L>= 48, not reach(495, L+1), not reach(590, L+1).
reach(497,L + 1) :- reach(590, L), L <= 52, L>= 49, not reach(592, L+1).
reach(592,L + 1) :- reach(590, L), L <= 52, L>= 49, not reach(497, L+1).
reach(498,L + 1) :- reach(499, L), L <= 54, L>= 54, not reach(500, L+1).
reach(500,L + 1) :- reach(499, L), L <= 54, L>= 54, not reach(498, L+1).
reach(499,L + 1) :- reach(500, L), L <= 53, L>= 53, not reach(501, L+1).
reach(501,L + 1) :- reach(500, L), L <= 53, L>= 53, not reach(499, L+1).
reach(500,L + 1) :- reach(501, L), L <= 52, L>= 52, not reach(502, L+1), not reach(503, L+1).
reach(502,L + 1) :- reach(501, L), L <= 52, L>= 52, not reach(500, L+1), not reach(503, L+1).
reach(503,L + 1) :- reach(501, L), L <= 52, L>= 52, not reach(500, L+1), not reach(502, L+1).
reach(501,L + 1) :- reach(502, L), L <= 51, L>= 51, not reach(504, L+1).
reach(504,L + 1) :- reach(502, L), L <= 51, L>= 51, not reach(501, L+1).
reach(502,L + 1) :- reach(504, L), L <= 51, L>= 50, not reach(505, L+1).
reach(505,L + 1) :- reach(504, L), L <= 51, L>= 50, not reach(502, L+1).
reach(504,L + 1) :- reach(505, L), L <= 52, L>= 49, not reach(660, L+1), not reach(673, L+1).
reach(660,L + 1) :- reach(505, L), L <= 52, L>= 49, not reach(504, L+1), not reach(673, L+1).
reach(673,L + 1) :- reach(505, L), L <= 52, L>= 49, not reach(504, L+1), not reach(660, L+1).
reach(505,L + 1) :- reach(673, L), L <= 53, L>= 50, not reach(672, L+1).
reach(672,L + 1) :- reach(673, L), L <= 53, L>= 50, not reach(505, L+1).
reach(518,L + 1) :- reach(517, L), L <= 46, L>= 45, not reach(519, L+1).
reach(519,L + 1) :- reach(517, L), L <= 46, L>= 45, not reach(518, L+1).
reach(517,L + 1) :- reach(518, L), L <= 47, L>= 44, not reach(598, L+1), not reach(594, L+1), not reach(596, L+1), not reach(681, L+1), not reach(731, L+1).
reach(598,L + 1) :- reach(518, L), L <= 47, L>= 44, not reach(517, L+1), not reach(594, L+1), not reach(596, L+1), not reach(681, L+1), not reach(731, L+1).
reach(594,L + 1) :- reach(518, L), L <= 47, L>= 44, not reach(517, L+1), not reach(598, L+1), not reach(596, L+1), not reach(681, L+1), not reach(731, L+1).
reach(596,L + 1) :- reach(518, L), L <= 47, L>= 44, not reach(517, L+1), not reach(598, L+1), not reach(594, L+1), not reach(681, L+1), not reach(731, L+1).
reach(681,L + 1) :- reach(518, L), L <= 47, L>= 44, not reach(517, L+1), not reach(598, L+1), not reach(594, L+1), not reach(596, L+1), not reach(731, L+1).
reach(731,L + 1) :- reach(518, L), L <= 47, L>= 44, not reach(517, L+1), not reach(598, L+1), not reach(594, L+1), not reach(596, L+1), not reach(681, L+1).
reach(517,L + 1) :- reach(519, L), L <= 45, L>= 44, not reach(738, L+1).
reach(738,L + 1) :- reach(519, L), L <= 45, L>= 44, not reach(517, L+1).
reach(518,L + 1) :- reach(598, L), L <= 46, L>= 45.
reach(518,L + 1) :- reach(596, L), L <= 46, L>= 44, not reach(560, L+1).
reach(560,L + 1) :- reach(596, L), L <= 46, L>= 44, not reach(518, L+1).
reach(518,L + 1) :- reach(731, L), L <= 46, L>= 43, not reach(737, L+1).
reach(737,L + 1) :- reach(731, L), L <= 46, L>= 43, not reach(518, L+1).
reach(519,L + 1) :- reach(738, L), L <= 44, L>= 43, not reach(737, L+1).
reach(737,L + 1) :- reach(738, L), L <= 44, L>= 43, not reach(519, L+1).
reach(552,L + 1) :- reach(726, L), L <= 43, L>= 42, not reach(727, L+1).
reach(727,L + 1) :- reach(726, L), L <= 43, L>= 42, not reach(552, L+1).
reach(560,L + 1) :- reach(595, L), L <= 47, L>= 44, not reach(593, L+1).
reach(593,L + 1) :- reach(595, L), L <= 47, L>= 44, not reach(560, L+1).
reach(560,L + 1) :- reach(597, L), L <= 45, L>= 42, not reach(723, L+1).
reach(723,L + 1) :- reach(597, L), L <= 45, L>= 42, not reach(560, L+1).
reach(562,L + 1) :- reach(561, L), L <= 42, L>= 40, not reach(563, L+1).
reach(563,L + 1) :- reach(561, L), L <= 42, L>= 40, not reach(562, L+1).
reach(561,L + 1) :- reach(562, L), L <= 41, L>= 39, not reach(564, L+1).
reach(564,L + 1) :- reach(562, L), L <= 41, L>= 39, not reach(561, L+1).
reach(561,L + 1) :- reach(563, L), L <= 43, L>= 41, not reach(677, L+1), not reach(678, L+1).
reach(677,L + 1) :- reach(563, L), L <= 43, L>= 41, not reach(561, L+1), not reach(678, L+1).
reach(678,L + 1) :- reach(563, L), L <= 43, L>= 41, not reach(561, L+1), not reach(677, L+1).
reach(562,L + 1) :- reach(564, L), L <= 40, L>= 38, not reach(599, L+1), not reach(600, L+1).
reach(599,L + 1) :- reach(564, L), L <= 40, L>= 38, not reach(562, L+1), not reach(600, L+1).
reach(600,L + 1) :- reach(564, L), L <= 40, L>= 38, not reach(562, L+1), not reach(599, L+1).
reach(563,L + 1) :- reach(677, L), L <= 44, L>= 42, not reach(697, L+1), not reach(698, L+1).
reach(697,L + 1) :- reach(677, L), L <= 44, L>= 42, not reach(563, L+1), not reach(698, L+1).
reach(698,L + 1) :- reach(677, L), L <= 44, L>= 42, not reach(563, L+1), not reach(697, L+1).
reach(564,L + 1) :- reach(600, L), L <= 41, L>= 39, not reach(601, L+1).
reach(601,L + 1) :- reach(600, L), L <= 41, L>= 39, not reach(564, L+1).
reach(574,L + 1) :- reach(648, L), L <= 23, L>= 22, not reach(649, L+1).
reach(649,L + 1) :- reach(648, L), L <= 23, L>= 22, not reach(574, L+1).
reach(590,L + 1) :- reach(592, L), L <= 53, L>= 50, not reach(674, L+1).
reach(674,L + 1) :- reach(592, L), L <= 53, L>= 50, not reach(590, L+1).
reach(592,L + 1) :- reach(674, L), L <= 54, L>= 51, not reach(675, L+1), not reach(676, L+1).
reach(675,L + 1) :- reach(674, L), L <= 54, L>= 51, not reach(592, L+1), not reach(676, L+1).
reach(676,L + 1) :- reach(674, L), L <= 54, L>= 51, not reach(592, L+1), not reach(675, L+1).
reach(600,L + 1) :- reach(601, L), L <= 42, L>= 40, not reach(602, L+1).
reach(602,L + 1) :- reach(601, L), L <= 42, L>= 40, not reach(600, L+1).
reach(601,L + 1) :- reach(602, L), L <= 43, L>= 41, not reach(603, L+1).
reach(603,L + 1) :- reach(602, L), L <= 43, L>= 41, not reach(601, L+1).
reach(602,L + 1) :- reach(603, L), L <= 44, L>= 42, not reach(604, L+1).
reach(604,L + 1) :- reach(603, L), L <= 44, L>= 42, not reach(602, L+1).
reach(603,L + 1) :- reach(604, L), L <= 45, L>= 43, not reach(605, L+1).
reach(605,L + 1) :- reach(604, L), L <= 45, L>= 43, not reach(603, L+1).
reach(652,L + 1) :- reach(654, L), L <= 28, L>= 28, not reach(655, L+1).
reach(655,L + 1) :- reach(654, L), L <= 28, L>= 28, not reach(652, L+1).
reach(662,L + 1) :- reach(661, L), L <= 22, L>= 22, not reach(663, L+1).
reach(663,L + 1) :- reach(661, L), L <= 22, L>= 22, not reach(662, L+1).
reach(671,L + 1) :- reach(672, L), L <= 54, L>= 51, not reach(673, L+1), not reach(675, L+1).
reach(673,L + 1) :- reach(672, L), L <= 54, L>= 51, not reach(671, L+1), not reach(675, L+1).
reach(675,L + 1) :- reach(672, L), L <= 54, L>= 51, not reach(671, L+1), not reach(673, L+1).
reach(672,L + 1) :- reach(675, L), L <= 53, L>= 52, not reach(674, L+1).
reach(674,L + 1) :- reach(675, L), L <= 53, L>= 52, not reach(672, L+1).
reach(677,L + 1) :- reach(697, L), L <= 44, L>= 43, not reach(698, L+1).
reach(698,L + 1) :- reach(697, L), L <= 44, L>= 43, not reach(677, L+1).
reach(678,L + 1) :- reach(705, L), L <= 44, L>= 43, not reach(698, L+1).
reach(698,L + 1) :- reach(705, L), L <= 44, L>= 43, not reach(678, L+1).
reach(698,L + 1) :- reach(699, L), L <= 44, L>= 43, not reach(700, L+1).
reach(700,L + 1) :- reach(699, L), L <= 44, L>= 43, not reach(698, L+1).
reach(699,L + 1) :- reach(700, L), L <= 43, L>= 42, not reach(720, L+1), not reach(721, L+1).
reach(720,L + 1) :- reach(700, L), L <= 43, L>= 42, not reach(699, L+1), not reach(721, L+1).
reach(721,L + 1) :- reach(700, L), L <= 43, L>= 42, not reach(699, L+1), not reach(720, L+1).
reach(726,L + 1) :- reach(727, L), L <= 44, L>= 43, not reach(737, L+1).
reach(737,L + 1) :- reach(727, L), L <= 44, L>= 43, not reach(726, L+1).

