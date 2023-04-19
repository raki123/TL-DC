reach(X, 0) :- start(X).
:- goal(X), not reach(X, 1), not reach(X, 2), not reach(X, 3), not reach(X, 4), not reach(X, 5), not reach(X, 6), not reach(X, 7), not reach(X, 8), not reach(X, 9), not reach(X, 10), not reach(X, 11), not reach(X, 12), not reach(X, 13), not reach(X, 14), not reach(X, 15), not reach(X, 16), not reach(X, 17), not reach(X, 18), not reach(X, 19), not reach(X, 20), not reach(X, 21), not reach(X, 22), not reach(X, 23), not reach(X, 24), not reach(X, 25), not reach(X, 26), not reach(X, 27), not reach(X, 28), not reach(X, 29), not reach(X, 30), not reach(X, 31), not reach(X, 32), not reach(X, 33), not reach(X, 34), not reach(X, 35), not reach(X, 36), not reach(X, 37), not reach(X, 38), not reach(X, 39), not reach(X, 40), not reach(X, 41), not reach(X, 42), not reach(X, 43), not reach(X, 44), not reach(X, 45), not reach(X, 46), not reach(X, 47), not reach(X, 48), not reach(X, 49), not reach(X, 50), not reach(X, 51), not reach(X, 52), not reach(X, 53), not reach(X, 54), not reach(X, 55), not reach(X, 56), not reach(X, 57), not reach(X, 58), not reach(X, 59), not reach(X, 60), not reach(X, 61).
start(25).
goal(28).
:- reach(X, L), reach(X, L'), L != L'.
reach(47,L + 1) :- reach(272, L), L <= 58, L>= 8, not reach(271, L+1), not reach(564, L+1).
reach(271,L + 1) :- reach(272, L), L <= 58, L>= 8, not reach(47, L+1), not reach(564, L+1).
reach(564,L + 1) :- reach(272, L), L <= 58, L>= 8, not reach(47, L+1), not reach(271, L+1).
reach(487,L + 1) :- reach(14, L), L <= 52, L>= 2, not reach(335, L+1), not reach(34, L+1), not reach(359, L+1).
reach(335,L + 1) :- reach(14, L), L <= 52, L>= 2, not reach(487, L+1), not reach(34, L+1), not reach(359, L+1).
reach(34,L + 1) :- reach(14, L), L <= 52, L>= 2, not reach(487, L+1), not reach(335, L+1), not reach(359, L+1).
reach(359,L + 1) :- reach(14, L), L <= 52, L>= 2, not reach(487, L+1), not reach(335, L+1), not reach(34, L+1).
reach(427,L + 1) :- reach(487, L), L <= 51, L>= 3, not reach(489, L+1), not reach(14, L+1), not reach(280, L+1).
reach(489,L + 1) :- reach(487, L), L <= 51, L>= 3, not reach(427, L+1), not reach(14, L+1), not reach(280, L+1).
reach(14,L + 1) :- reach(487, L), L <= 51, L>= 3, not reach(427, L+1), not reach(489, L+1), not reach(280, L+1).
reach(280,L + 1) :- reach(487, L), L <= 51, L>= 3, not reach(427, L+1), not reach(489, L+1), not reach(14, L+1).
reach(657,L + 1) :- reach(19, L), L <= 54, L>= 4, not reach(348, L+1), not reach(404, L+1), not reach(359, L+1), not reach(209, L+1), not reach(70, L+1).
reach(348,L + 1) :- reach(19, L), L <= 54, L>= 4, not reach(657, L+1), not reach(404, L+1), not reach(359, L+1), not reach(209, L+1), not reach(70, L+1).
reach(404,L + 1) :- reach(19, L), L <= 54, L>= 4, not reach(657, L+1), not reach(348, L+1), not reach(359, L+1), not reach(209, L+1), not reach(70, L+1).
reach(359,L + 1) :- reach(19, L), L <= 54, L>= 4, not reach(657, L+1), not reach(348, L+1), not reach(404, L+1), not reach(209, L+1), not reach(70, L+1).
reach(209,L + 1) :- reach(19, L), L <= 54, L>= 4, not reach(657, L+1), not reach(348, L+1), not reach(404, L+1), not reach(359, L+1), not reach(70, L+1).
reach(70,L + 1) :- reach(19, L), L <= 54, L>= 4, not reach(657, L+1), not reach(348, L+1), not reach(404, L+1), not reach(359, L+1), not reach(209, L+1).
reach(19,L + 1) :- reach(348, L), L <= 53, L>= 5, not reach(352, L+1), not reach(404, L+1), not reach(746, L+1), not reach(55, L+1).
reach(352,L + 1) :- reach(348, L), L <= 53, L>= 5, not reach(19, L+1), not reach(404, L+1), not reach(746, L+1), not reach(55, L+1).
reach(404,L + 1) :- reach(348, L), L <= 53, L>= 5, not reach(19, L+1), not reach(352, L+1), not reach(746, L+1), not reach(55, L+1).
reach(746,L + 1) :- reach(348, L), L <= 53, L>= 5, not reach(19, L+1), not reach(352, L+1), not reach(404, L+1), not reach(55, L+1).
reach(55,L + 1) :- reach(348, L), L <= 53, L>= 5, not reach(19, L+1), not reach(352, L+1), not reach(404, L+1), not reach(746, L+1).
reach(19,L + 1) :- reach(657, L), L <= 53, L>= 4, not reach(268, L+1), not reach(266, L+1), not reach(359, L+1).
reach(268,L + 1) :- reach(657, L), L <= 53, L>= 4, not reach(19, L+1), not reach(266, L+1), not reach(359, L+1).
reach(266,L + 1) :- reach(657, L), L <= 53, L>= 4, not reach(19, L+1), not reach(268, L+1), not reach(359, L+1).
reach(359,L + 1) :- reach(657, L), L <= 53, L>= 4, not reach(19, L+1), not reach(268, L+1), not reach(266, L+1).
reach(36,L + 1) :- reach(25, L), L <= 50, L>= 0, not reach(34, L+1).
reach(34,L + 1) :- reach(25, L), L <= 50, L>= 0, not reach(36, L+1).
reach(36,L + 1) :- reach(34, L), L <= 51, L>= 1, not reach(25, L+1), not reach(14, L+1), not reach(335, L+1).
reach(14,L + 1) :- reach(34, L), L <= 51, L>= 1, not reach(36, L+1), not reach(25, L+1), not reach(335, L+1).
reach(335,L + 1) :- reach(34, L), L <= 51, L>= 1, not reach(36, L+1), not reach(25, L+1), not reach(14, L+1).
reach(34,L + 1) :- reach(36, L), L <= 50, L>= 1, not reach(25, L+1), not reach(488, L+1).
reach(488,L + 1) :- reach(36, L), L <= 50, L>= 1, not reach(34, L+1), not reach(25, L+1).
reach(43,L + 1) :- reach(39, L), L <= 50, L>= 11, not reach(98, L+1), not reach(307, L+1), not reach(141, L+1).
reach(98,L + 1) :- reach(39, L), L <= 50, L>= 11, not reach(43, L+1), not reach(307, L+1), not reach(141, L+1).
reach(307,L + 1) :- reach(39, L), L <= 50, L>= 11, not reach(43, L+1), not reach(98, L+1), not reach(141, L+1).
reach(141,L + 1) :- reach(39, L), L <= 50, L>= 11, not reach(43, L+1), not reach(98, L+1), not reach(307, L+1).
reach(39,L + 1) :- reach(43, L), L <= 51, L>= 12, not reach(59, L+1), not reach(466, L+1).
reach(59,L + 1) :- reach(43, L), L <= 51, L>= 12, not reach(39, L+1), not reach(466, L+1).
reach(466,L + 1) :- reach(43, L), L <= 51, L>= 12, not reach(39, L+1), not reach(59, L+1).
reach(303,L + 1) :- reach(98, L), L <= 51, L>= 10, not reach(39, L+1), not reach(302, L+1).
reach(39,L + 1) :- reach(98, L), L <= 51, L>= 10, not reach(303, L+1), not reach(302, L+1).
reach(302,L + 1) :- reach(98, L), L <= 51, L>= 10, not reach(303, L+1), not reach(39, L+1).
reach(143,L + 1) :- reach(141, L), L <= 51, L>= 10, not reach(176, L+1), not reach(39, L+1).
reach(176,L + 1) :- reach(141, L), L <= 51, L>= 10, not reach(143, L+1), not reach(39, L+1).
reach(39,L + 1) :- reach(141, L), L <= 51, L>= 10, not reach(143, L+1), not reach(176, L+1).
reach(43,L + 1) :- reach(59, L), L <= 52, L>= 13, not reach(61, L+1), not reach(108, L+1).
reach(61,L + 1) :- reach(59, L), L <= 52, L>= 13, not reach(43, L+1), not reach(108, L+1).
reach(108,L + 1) :- reach(59, L), L <= 52, L>= 13, not reach(43, L+1), not reach(61, L+1).
reach(107,L + 1) :- reach(307, L), L <= 51, L>= 11, not reach(39, L+1), not reach(310, L+1).
reach(39,L + 1) :- reach(307, L), L <= 51, L>= 11, not reach(107, L+1), not reach(310, L+1).
reach(310,L + 1) :- reach(307, L), L <= 51, L>= 11, not reach(107, L+1), not reach(39, L+1).
reach(272,L + 1) :- reach(47, L), L <= 57, L>= 7, not reach(271, L+1), not reach(493, L+1).
reach(271,L + 1) :- reach(47, L), L <= 57, L>= 7, not reach(272, L+1), not reach(493, L+1).
reach(493,L + 1) :- reach(47, L), L <= 57, L>= 7, not reach(272, L+1), not reach(271, L+1).
reach(352,L + 1) :- reach(55, L), L <= 52, L>= 6, not reach(57, L+1), not reach(348, L+1).
reach(57,L + 1) :- reach(55, L), L <= 52, L>= 6, not reach(352, L+1), not reach(348, L+1).
reach(348,L + 1) :- reach(55, L), L <= 52, L>= 6, not reach(352, L+1), not reach(57, L+1).
reach(55,L + 1) :- reach(352, L), L <= 52, L>= 6, not reach(315, L+1), not reach(58, L+1), not reach(348, L+1).
reach(315,L + 1) :- reach(352, L), L <= 52, L>= 6, not reach(55, L+1), not reach(58, L+1), not reach(348, L+1).
reach(58,L + 1) :- reach(352, L), L <= 52, L>= 6, not reach(55, L+1), not reach(315, L+1), not reach(348, L+1).
reach(348,L + 1) :- reach(352, L), L <= 52, L>= 6, not reach(55, L+1), not reach(315, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(57, L), L <= 51, L>= 7, not reach(314, L+1), not reach(542, L+1), not reach(55, L+1).
reach(314,L + 1) :- reach(57, L), L <= 51, L>= 7, not reach(58, L+1), not reach(542, L+1), not reach(55, L+1).
reach(542,L + 1) :- reach(57, L), L <= 51, L>= 7, not reach(58, L+1), not reach(314, L+1), not reach(55, L+1).
reach(55,L + 1) :- reach(57, L), L <= 51, L>= 7, not reach(58, L+1), not reach(314, L+1), not reach(542, L+1).
reach(57,L + 1) :- reach(58, L), L <= 51, L>= 7, not reach(352, L+1), not reach(746, L+1), not reach(318, L+1).
reach(352,L + 1) :- reach(58, L), L <= 51, L>= 7, not reach(57, L+1), not reach(746, L+1), not reach(318, L+1).
reach(746,L + 1) :- reach(58, L), L <= 51, L>= 7, not reach(57, L+1), not reach(352, L+1), not reach(318, L+1).
reach(318,L + 1) :- reach(58, L), L <= 51, L>= 7, not reach(57, L+1), not reach(352, L+1), not reach(746, L+1).
reach(59,L + 1) :- reach(61, L), L <= 51, L>= 14, not reach(104, L+1), not reach(109, L+1).
reach(104,L + 1) :- reach(61, L), L <= 51, L>= 14, not reach(59, L+1), not reach(109, L+1).
reach(109,L + 1) :- reach(61, L), L <= 51, L>= 14, not reach(59, L+1), not reach(104, L+1).
reach(88,L + 1) :- reach(108, L), L <= 53, L>= 13, not reach(106, L+1), not reach(59, L+1).
reach(106,L + 1) :- reach(108, L), L <= 53, L>= 13, not reach(88, L+1), not reach(59, L+1).
reach(59,L + 1) :- reach(108, L), L <= 53, L>= 13, not reach(88, L+1), not reach(106, L+1).
reach(61,L + 1) :- reach(104, L), L <= 52, L>= 15, not reach(109, L+1), not reach(100, L+1).
reach(109,L + 1) :- reach(104, L), L <= 52, L>= 15, not reach(61, L+1), not reach(100, L+1).
reach(100,L + 1) :- reach(104, L), L <= 52, L>= 15, not reach(61, L+1), not reach(109, L+1).
reach(61,L + 1) :- reach(109, L), L <= 51, L>= 15, not reach(104, L+1), not reach(101, L+1).
reach(104,L + 1) :- reach(109, L), L <= 51, L>= 15, not reach(61, L+1), not reach(101, L+1).
reach(101,L + 1) :- reach(109, L), L <= 51, L>= 15, not reach(61, L+1), not reach(104, L+1).
reach(90,L + 1) :- reach(70, L), L <= 54, L>= 5, not reach(663, L+1), not reach(92, L+1), not reach(195, L+1), not reach(19, L+1).
reach(663,L + 1) :- reach(70, L), L <= 54, L>= 5, not reach(90, L+1), not reach(92, L+1), not reach(195, L+1), not reach(19, L+1).
reach(92,L + 1) :- reach(70, L), L <= 54, L>= 5, not reach(90, L+1), not reach(663, L+1), not reach(195, L+1), not reach(19, L+1).
reach(195,L + 1) :- reach(70, L), L <= 54, L>= 5, not reach(90, L+1), not reach(663, L+1), not reach(92, L+1), not reach(19, L+1).
reach(19,L + 1) :- reach(70, L), L <= 54, L>= 5, not reach(90, L+1), not reach(663, L+1), not reach(92, L+1), not reach(195, L+1).
reach(88,L + 1) :- reach(73, L), L <= 55, L>= 14, not reach(392, L+1), not reach(93, L+1).
reach(392,L + 1) :- reach(73, L), L <= 55, L>= 14, not reach(88, L+1), not reach(93, L+1).
reach(93,L + 1) :- reach(73, L), L <= 55, L>= 14, not reach(88, L+1), not reach(392, L+1).
reach(73,L + 1) :- reach(88, L), L <= 54, L>= 13, not reach(100, L+1), not reach(108, L+1), not reach(106, L+1), not reach(83, L+1).
reach(100,L + 1) :- reach(88, L), L <= 54, L>= 13, not reach(73, L+1), not reach(108, L+1), not reach(106, L+1), not reach(83, L+1).
reach(108,L + 1) :- reach(88, L), L <= 54, L>= 13, not reach(73, L+1), not reach(100, L+1), not reach(106, L+1), not reach(83, L+1).
reach(106,L + 1) :- reach(88, L), L <= 54, L>= 13, not reach(73, L+1), not reach(100, L+1), not reach(108, L+1), not reach(83, L+1).
reach(83,L + 1) :- reach(88, L), L <= 54, L>= 13, not reach(73, L+1), not reach(100, L+1), not reach(108, L+1), not reach(106, L+1).
reach(73,L + 1) :- reach(392, L), L <= 56, L>= 13, not reach(387, L+1), not reach(395, L+1).
reach(387,L + 1) :- reach(392, L), L <= 56, L>= 13, not reach(73, L+1), not reach(395, L+1).
reach(395,L + 1) :- reach(392, L), L <= 56, L>= 13, not reach(73, L+1), not reach(387, L+1).
reach(92,L + 1) :- reach(90, L), L <= 53, L>= 6, not reach(70, L+1), not reach(91, L+1).
reach(70,L + 1) :- reach(90, L), L <= 53, L>= 6, not reach(92, L+1), not reach(91, L+1).
reach(91,L + 1) :- reach(90, L), L <= 53, L>= 6, not reach(92, L+1), not reach(70, L+1).
reach(493,L + 1) :- reach(663, L), L <= 55, L>= 6, not reach(70, L+1), not reach(209, L+1).
reach(70,L + 1) :- reach(663, L), L <= 55, L>= 6, not reach(493, L+1), not reach(209, L+1).
reach(209,L + 1) :- reach(663, L), L <= 55, L>= 6, not reach(493, L+1), not reach(70, L+1).
reach(85,L + 1) :- reach(83, L), L <= 53, L>= 14, not reach(86, L+1), not reach(88, L+1).
reach(86,L + 1) :- reach(83, L), L <= 53, L>= 14, not reach(85, L+1), not reach(88, L+1).
reach(88,L + 1) :- reach(83, L), L <= 53, L>= 14, not reach(85, L+1), not reach(86, L+1).
reach(83,L + 1) :- reach(85, L), L <= 54, L>= 15, not reach(93, L+1), not reach(94, L+1).
reach(93,L + 1) :- reach(85, L), L <= 54, L>= 15, not reach(83, L+1), not reach(94, L+1).
reach(94,L + 1) :- reach(85, L), L <= 54, L>= 15, not reach(83, L+1), not reach(93, L+1).
reach(83,L + 1) :- reach(86, L), L <= 52, L>= 15, not reach(704, L+1), not reach(153, L+1), not reach(101, L+1), not reach(434, L+1).
reach(704,L + 1) :- reach(86, L), L <= 52, L>= 15, not reach(83, L+1), not reach(153, L+1), not reach(101, L+1), not reach(434, L+1).
reach(153,L + 1) :- reach(86, L), L <= 52, L>= 15, not reach(83, L+1), not reach(704, L+1), not reach(101, L+1), not reach(434, L+1).
reach(101,L + 1) :- reach(86, L), L <= 52, L>= 15, not reach(83, L+1), not reach(704, L+1), not reach(153, L+1), not reach(434, L+1).
reach(434,L + 1) :- reach(86, L), L <= 52, L>= 15, not reach(83, L+1), not reach(704, L+1), not reach(153, L+1), not reach(101, L+1).
reach(85,L + 1) :- reach(93, L), L <= 55, L>= 14, not reach(73, L+1), not reach(94, L+1), not reach(387, L+1).
reach(73,L + 1) :- reach(93, L), L <= 55, L>= 14, not reach(85, L+1), not reach(94, L+1), not reach(387, L+1).
reach(94,L + 1) :- reach(93, L), L <= 55, L>= 14, not reach(85, L+1), not reach(73, L+1), not reach(387, L+1).
reach(387,L + 1) :- reach(93, L), L <= 55, L>= 14, not reach(85, L+1), not reach(73, L+1), not reach(94, L+1).
reach(85,L + 1) :- reach(94, L), L <= 55, L>= 14, not reach(390, L+1), not reach(93, L+1).
reach(390,L + 1) :- reach(94, L), L <= 55, L>= 14, not reach(85, L+1), not reach(93, L+1).
reach(93,L + 1) :- reach(94, L), L <= 55, L>= 14, not reach(85, L+1), not reach(390, L+1).
reach(86,L + 1) :- reach(704, L), L <= 51, L>= 16, not reach(214, L+1), not reach(703, L+1).
reach(214,L + 1) :- reach(704, L), L <= 51, L>= 16, not reach(86, L+1), not reach(703, L+1).
reach(703,L + 1) :- reach(704, L), L <= 51, L>= 16, not reach(86, L+1), not reach(214, L+1).
reach(88,L + 1) :- reach(100, L), L <= 53, L>= 14, not reach(101, L+1), not reach(104, L+1).
reach(101,L + 1) :- reach(100, L), L <= 53, L>= 14, not reach(88, L+1), not reach(104, L+1).
reach(104,L + 1) :- reach(100, L), L <= 53, L>= 14, not reach(88, L+1), not reach(101, L+1).
reach(88,L + 1) :- reach(106, L), L <= 53, L>= 12, not reach(108, L+1), not reach(107, L+1).
reach(108,L + 1) :- reach(106, L), L <= 53, L>= 12, not reach(88, L+1), not reach(107, L+1).
reach(107,L + 1) :- reach(106, L), L <= 53, L>= 12, not reach(88, L+1), not reach(108, L+1).
reach(266,L + 1) :- reach(91, L), L <= 52, L>= 6, not reach(92, L+1), not reach(90, L+1).
reach(92,L + 1) :- reach(91, L), L <= 52, L>= 6, not reach(266, L+1), not reach(90, L+1).
reach(90,L + 1) :- reach(91, L), L <= 52, L>= 6, not reach(266, L+1), not reach(92, L+1).
reach(90,L + 1) :- reach(92, L), L <= 53, L>= 5, not reach(91, L+1), not reach(268, L+1), not reach(70, L+1), not reach(489, L+1).
reach(91,L + 1) :- reach(92, L), L <= 53, L>= 5, not reach(90, L+1), not reach(268, L+1), not reach(70, L+1), not reach(489, L+1).
reach(268,L + 1) :- reach(92, L), L <= 53, L>= 5, not reach(90, L+1), not reach(91, L+1), not reach(70, L+1), not reach(489, L+1).
reach(70,L + 1) :- reach(92, L), L <= 53, L>= 5, not reach(90, L+1), not reach(91, L+1), not reach(268, L+1), not reach(489, L+1).
reach(489,L + 1) :- reach(92, L), L <= 53, L>= 5, not reach(90, L+1), not reach(91, L+1), not reach(268, L+1), not reach(70, L+1).
reach(91,L + 1) :- reach(266, L), L <= 52, L>= 5, not reach(268, L+1), not reach(657, L+1).
reach(268,L + 1) :- reach(266, L), L <= 52, L>= 5, not reach(91, L+1), not reach(657, L+1).
reach(657,L + 1) :- reach(266, L), L <= 52, L>= 5, not reach(91, L+1), not reach(268, L+1).
reach(92,L + 1) :- reach(268, L), L <= 52, L>= 5, not reach(266, L+1), not reach(657, L+1).
reach(266,L + 1) :- reach(268, L), L <= 52, L>= 5, not reach(92, L+1), not reach(657, L+1).
reach(657,L + 1) :- reach(268, L), L <= 52, L>= 5, not reach(92, L+1), not reach(266, L+1).
reach(94,L + 1) :- reach(390, L), L <= 56, L>= 13, not reach(395, L+1), not reach(397, L+1).
reach(395,L + 1) :- reach(390, L), L <= 56, L>= 13, not reach(94, L+1), not reach(397, L+1).
reach(397,L + 1) :- reach(390, L), L <= 56, L>= 13, not reach(94, L+1), not reach(395, L+1).
reach(392,L + 1) :- reach(387, L), L <= 56, L>= 13, not reach(93, L+1), not reach(395, L+1).
reach(93,L + 1) :- reach(387, L), L <= 56, L>= 13, not reach(392, L+1), not reach(395, L+1).
reach(395,L + 1) :- reach(387, L), L <= 56, L>= 13, not reach(392, L+1), not reach(93, L+1).
reach(98,L + 1) :- reach(303, L), L <= 52, L>= 9, not reach(191, L+1), not reach(298, L+1).
reach(191,L + 1) :- reach(303, L), L <= 52, L>= 9, not reach(98, L+1), not reach(298, L+1).
reach(298,L + 1) :- reach(303, L), L <= 52, L>= 9, not reach(98, L+1), not reach(191, L+1).
reach(86,L + 1) :- reach(101, L), L <= 52, L>= 15, not reach(100, L+1), not reach(109, L+1).
reach(100,L + 1) :- reach(101, L), L <= 52, L>= 15, not reach(86, L+1), not reach(109, L+1).
reach(109,L + 1) :- reach(101, L), L <= 52, L>= 15, not reach(86, L+1), not reach(100, L+1).
reach(307,L + 1) :- reach(107, L), L <= 52, L>= 11, not reach(106, L+1), not reach(310, L+1).
reach(106,L + 1) :- reach(107, L), L <= 52, L>= 11, not reach(307, L+1), not reach(310, L+1).
reach(310,L + 1) :- reach(107, L), L <= 52, L>= 11, not reach(307, L+1), not reach(106, L+1).
reach(291,L + 1) :- reach(327, L), L <= 50, L>= 11, not reach(325, L+1), not reach(329, L+1), not reach(250, L+1).
reach(325,L + 1) :- reach(327, L), L <= 50, L>= 11, not reach(291, L+1), not reach(329, L+1), not reach(250, L+1).
reach(329,L + 1) :- reach(327, L), L <= 50, L>= 11, not reach(291, L+1), not reach(325, L+1), not reach(250, L+1).
reach(250,L + 1) :- reach(327, L), L <= 50, L>= 11, not reach(291, L+1), not reach(325, L+1), not reach(329, L+1).
reach(254,L + 1) :- reach(117, L), L <= 60, L>= 12, not reach(28, L+1), not reach(501, L+1).
reach(28,L + 1) :- reach(117, L), L <= 60, L>= 12, not reach(254, L+1), not reach(501, L+1).
reach(501,L + 1) :- reach(117, L), L <= 60, L>= 12, not reach(254, L+1), not reach(28, L+1).
reach(117,L + 1) :- reach(254, L), L <= 59, L>= 13, not reach(672, L+1), not reach(674, L+1).
reach(672,L + 1) :- reach(254, L), L <= 59, L>= 13, not reach(117, L+1), not reach(674, L+1).
reach(674,L + 1) :- reach(254, L), L <= 59, L>= 13, not reach(117, L+1), not reach(672, L+1).
reach(688,L + 1) :- reach(547, L), L <= 50, L>= 11, not reach(548, L+1), not reach(127, L+1), not reach(241, L+1).
reach(548,L + 1) :- reach(547, L), L <= 50, L>= 11, not reach(688, L+1), not reach(127, L+1), not reach(241, L+1).
reach(127,L + 1) :- reach(547, L), L <= 50, L>= 11, not reach(688, L+1), not reach(548, L+1), not reach(241, L+1).
reach(241,L + 1) :- reach(547, L), L <= 50, L>= 11, not reach(688, L+1), not reach(548, L+1), not reach(127, L+1).
reach(213,L + 1) :- reach(125, L), L <= 51, L>= 16, not reach(136, L+1), not reach(516, L+1).
reach(136,L + 1) :- reach(125, L), L <= 51, L>= 16, not reach(213, L+1), not reach(516, L+1).
reach(516,L + 1) :- reach(125, L), L <= 51, L>= 16, not reach(213, L+1), not reach(136, L+1).
reach(524,L + 1) :- reach(516, L), L <= 52, L>= 15, not reach(630, L+1), not reach(125, L+1).
reach(630,L + 1) :- reach(516, L), L <= 52, L>= 15, not reach(524, L+1), not reach(125, L+1).
reach(125,L + 1) :- reach(516, L), L <= 52, L>= 15, not reach(524, L+1), not reach(630, L+1).
reach(125,L + 1) :- reach(213, L), L <= 50, L>= 17, not reach(134, L+1), not reach(136, L+1).
reach(134,L + 1) :- reach(213, L), L <= 50, L>= 17, not reach(125, L+1), not reach(136, L+1).
reach(136,L + 1) :- reach(213, L), L <= 50, L>= 17, not reach(125, L+1), not reach(134, L+1).
reach(125,L + 1) :- reach(136, L), L <= 50, L>= 17, not reach(133, L+1), not reach(213, L+1).
reach(133,L + 1) :- reach(136, L), L <= 50, L>= 17, not reach(125, L+1), not reach(213, L+1).
reach(213,L + 1) :- reach(136, L), L <= 50, L>= 17, not reach(125, L+1), not reach(133, L+1).
reach(128,L + 1) :- reach(127, L), L <= 51, L>= 12, not reach(292, L+1), not reach(217, L+1), not reach(547, L+1).
reach(292,L + 1) :- reach(127, L), L <= 51, L>= 12, not reach(128, L+1), not reach(217, L+1), not reach(547, L+1).
reach(217,L + 1) :- reach(127, L), L <= 51, L>= 12, not reach(128, L+1), not reach(292, L+1), not reach(547, L+1).
reach(547,L + 1) :- reach(127, L), L <= 51, L>= 12, not reach(128, L+1), not reach(292, L+1), not reach(217, L+1).
reach(127,L + 1) :- reach(128, L), L <= 51, L>= 13, not reach(217, L+1), not reach(285, L+1).
reach(217,L + 1) :- reach(128, L), L <= 51, L>= 13, not reach(127, L+1), not reach(285, L+1).
reach(285,L + 1) :- reach(128, L), L <= 51, L>= 13, not reach(127, L+1), not reach(217, L+1).
reach(134,L + 1) :- reach(133, L), L <= 49, L>= 18, not reach(136, L+1), not reach(452, L+1).
reach(136,L + 1) :- reach(133, L), L <= 49, L>= 18, not reach(134, L+1), not reach(452, L+1).
reach(452,L + 1) :- reach(133, L), L <= 49, L>= 18, not reach(134, L+1), not reach(136, L+1).
reach(133,L + 1) :- reach(134, L), L <= 49, L>= 18, not reach(213, L+1), not reach(214, L+1).
reach(213,L + 1) :- reach(134, L), L <= 49, L>= 18, not reach(133, L+1), not reach(214, L+1).
reach(214,L + 1) :- reach(134, L), L <= 49, L>= 18, not reach(133, L+1), not reach(213, L+1).
reach(134,L + 1) :- reach(214, L), L <= 50, L>= 17, not reach(703, L+1), not reach(704, L+1).
reach(703,L + 1) :- reach(214, L), L <= 50, L>= 17, not reach(134, L+1), not reach(704, L+1).
reach(704,L + 1) :- reach(214, L), L <= 50, L>= 17, not reach(134, L+1), not reach(703, L+1).
reach(515,L + 1) :- reach(452, L), L <= 50, L>= 19, not reach(133, L+1), not reach(453, L+1).
reach(133,L + 1) :- reach(452, L), L <= 50, L>= 19, not reach(515, L+1), not reach(453, L+1).
reach(453,L + 1) :- reach(452, L), L <= 50, L>= 19, not reach(515, L+1), not reach(133, L+1).
reach(128,L + 1) :- reach(217, L), L <= 52, L>= 13, not reach(127, L+1), not reach(630, L+1).
reach(127,L + 1) :- reach(217, L), L <= 52, L>= 13, not reach(128, L+1), not reach(630, L+1).
reach(630,L + 1) :- reach(217, L), L <= 52, L>= 13, not reach(128, L+1), not reach(127, L+1).
reach(170,L + 1) :- reach(250, L), L <= 51, L>= 10, not reach(355, L+1), not reach(249, L+1), not reach(329, L+1), not reach(143, L+1), not reach(327, L+1).
reach(355,L + 1) :- reach(250, L), L <= 51, L>= 10, not reach(170, L+1), not reach(249, L+1), not reach(329, L+1), not reach(143, L+1), not reach(327, L+1).
reach(249,L + 1) :- reach(250, L), L <= 51, L>= 10, not reach(170, L+1), not reach(355, L+1), not reach(329, L+1), not reach(143, L+1), not reach(327, L+1).
reach(329,L + 1) :- reach(250, L), L <= 51, L>= 10, not reach(170, L+1), not reach(355, L+1), not reach(249, L+1), not reach(143, L+1), not reach(327, L+1).
reach(143,L + 1) :- reach(250, L), L <= 51, L>= 10, not reach(170, L+1), not reach(355, L+1), not reach(249, L+1), not reach(329, L+1), not reach(327, L+1).
reach(327,L + 1) :- reach(250, L), L <= 51, L>= 10, not reach(170, L+1), not reach(355, L+1), not reach(249, L+1), not reach(329, L+1), not reach(143, L+1).
reach(141,L + 1) :- reach(143, L), L <= 52, L>= 9, not reach(176, L+1), not reach(191, L+1), not reach(250, L+1), not reach(369, L+1).
reach(176,L + 1) :- reach(143, L), L <= 52, L>= 9, not reach(141, L+1), not reach(191, L+1), not reach(250, L+1), not reach(369, L+1).
reach(191,L + 1) :- reach(143, L), L <= 52, L>= 9, not reach(141, L+1), not reach(176, L+1), not reach(250, L+1), not reach(369, L+1).
reach(250,L + 1) :- reach(143, L), L <= 52, L>= 9, not reach(141, L+1), not reach(176, L+1), not reach(191, L+1), not reach(369, L+1).
reach(369,L + 1) :- reach(143, L), L <= 52, L>= 9, not reach(141, L+1), not reach(176, L+1), not reach(191, L+1), not reach(250, L+1).
reach(143,L + 1) :- reach(176, L), L <= 51, L>= 10, not reach(177, L+1), not reach(141, L+1).
reach(177,L + 1) :- reach(176, L), L <= 51, L>= 10, not reach(143, L+1), not reach(141, L+1).
reach(141,L + 1) :- reach(176, L), L <= 51, L>= 10, not reach(143, L+1), not reach(177, L+1).
reach(150,L + 1) :- reach(145, L), L <= 57, L>= 12, not reach(505, L+1), not reach(698, L+1).
reach(505,L + 1) :- reach(145, L), L <= 57, L>= 12, not reach(150, L+1), not reach(698, L+1).
reach(698,L + 1) :- reach(145, L), L <= 57, L>= 12, not reach(150, L+1), not reach(505, L+1).
reach(560,L + 1) :- reach(150, L), L <= 56, L>= 13, not reach(145, L+1), not reach(151, L+1).
reach(145,L + 1) :- reach(150, L), L <= 56, L>= 13, not reach(560, L+1), not reach(151, L+1).
reach(151,L + 1) :- reach(150, L), L <= 56, L>= 13, not reach(560, L+1), not reach(145, L+1).
reach(150,L + 1) :- reach(151, L), L <= 56, L>= 14, not reach(497, L+1), not reach(518, L+1), not reach(560, L+1).
reach(497,L + 1) :- reach(151, L), L <= 56, L>= 14, not reach(150, L+1), not reach(518, L+1), not reach(560, L+1).
reach(518,L + 1) :- reach(151, L), L <= 56, L>= 14, not reach(150, L+1), not reach(497, L+1), not reach(560, L+1).
reach(560,L + 1) :- reach(151, L), L <= 56, L>= 14, not reach(150, L+1), not reach(497, L+1), not reach(518, L+1).
reach(150,L + 1) :- reach(560, L), L <= 56, L>= 13, not reach(518, L+1), not reach(151, L+1), not reach(395, L+1).
reach(518,L + 1) :- reach(560, L), L <= 56, L>= 13, not reach(150, L+1), not reach(151, L+1), not reach(395, L+1).
reach(151,L + 1) :- reach(560, L), L <= 56, L>= 13, not reach(150, L+1), not reach(518, L+1), not reach(395, L+1).
reach(395,L + 1) :- reach(560, L), L <= 56, L>= 13, not reach(150, L+1), not reach(518, L+1), not reach(151, L+1).
reach(86,L + 1) :- reach(153, L), L <= 53, L>= 16, not reach(157, L+1), not reach(159, L+1), not reach(434, L+1).
reach(157,L + 1) :- reach(153, L), L <= 53, L>= 16, not reach(86, L+1), not reach(159, L+1), not reach(434, L+1).
reach(159,L + 1) :- reach(153, L), L <= 53, L>= 16, not reach(86, L+1), not reach(157, L+1), not reach(434, L+1).
reach(434,L + 1) :- reach(153, L), L <= 53, L>= 16, not reach(86, L+1), not reach(157, L+1), not reach(159, L+1).
reach(552,L + 1) :- reach(157, L), L <= 54, L>= 15, not reach(737, L+1), not reach(153, L+1).
reach(737,L + 1) :- reach(157, L), L <= 54, L>= 15, not reach(552, L+1), not reach(153, L+1).
reach(153,L + 1) :- reach(157, L), L <= 54, L>= 15, not reach(552, L+1), not reach(737, L+1).
reach(157,L + 1) :- reach(552, L), L <= 55, L>= 14, not reach(397, L+1), not reach(737, L+1).
reach(397,L + 1) :- reach(552, L), L <= 55, L>= 14, not reach(157, L+1), not reach(737, L+1).
reach(737,L + 1) :- reach(552, L), L <= 55, L>= 14, not reach(157, L+1), not reach(397, L+1).
reach(157,L + 1) :- reach(737, L), L <= 54, L>= 15, not reach(518, L+1), not reach(552, L+1).
reach(518,L + 1) :- reach(737, L), L <= 54, L>= 15, not reach(157, L+1), not reach(552, L+1).
reach(552,L + 1) :- reach(737, L), L <= 54, L>= 15, not reach(157, L+1), not reach(518, L+1).
reach(262,L + 1) :- reach(159, L), L <= 52, L>= 17, not reach(259, L+1), not reach(153, L+1).
reach(259,L + 1) :- reach(159, L), L <= 52, L>= 17, not reach(262, L+1), not reach(153, L+1).
reach(153,L + 1) :- reach(159, L), L <= 52, L>= 17, not reach(262, L+1), not reach(259, L+1).
reach(259,L + 1) :- reach(262, L), L <= 53, L>= 16, not reach(159, L+1), not reach(233, L+1), not reach(681, L+1).
reach(159,L + 1) :- reach(262, L), L <= 53, L>= 16, not reach(259, L+1), not reach(233, L+1), not reach(681, L+1).
reach(233,L + 1) :- reach(262, L), L <= 53, L>= 16, not reach(259, L+1), not reach(159, L+1), not reach(681, L+1).
reach(681,L + 1) :- reach(262, L), L <= 53, L>= 16, not reach(259, L+1), not reach(159, L+1), not reach(233, L+1).
reach(563,L + 1) :- reach(678, L), L <= 59, L>= 11, not reach(383, L+1), not reach(698, L+1).
reach(383,L + 1) :- reach(678, L), L <= 59, L>= 11, not reach(563, L+1), not reach(698, L+1).
reach(698,L + 1) :- reach(678, L), L <= 59, L>= 11, not reach(563, L+1), not reach(383, L+1).
reach(249,L + 1) :- reach(170, L), L <= 50, L>= 9, not reach(250, L+1), not reach(314, L+1).
reach(250,L + 1) :- reach(170, L), L <= 50, L>= 9, not reach(249, L+1), not reach(314, L+1).
reach(314,L + 1) :- reach(170, L), L <= 50, L>= 9, not reach(249, L+1), not reach(250, L+1).
reach(170,L + 1) :- reach(249, L), L <= 50, L>= 9, not reach(250, L+1), not reach(542, L+1).
reach(250,L + 1) :- reach(249, L), L <= 50, L>= 9, not reach(170, L+1), not reach(542, L+1).
reach(542,L + 1) :- reach(249, L), L <= 50, L>= 9, not reach(170, L+1), not reach(250, L+1).
reach(315,L + 1) :- reach(314, L), L <= 50, L>= 8, not reach(57, L+1), not reach(170, L+1).
reach(57,L + 1) :- reach(314, L), L <= 50, L>= 8, not reach(315, L+1), not reach(170, L+1).
reach(170,L + 1) :- reach(314, L), L <= 50, L>= 8, not reach(315, L+1), not reach(57, L+1).
reach(176,L + 1) :- reach(177, L), L <= 50, L>= 11, not reach(463, L+1), not reach(466, L+1), not reach(329, L+1).
reach(463,L + 1) :- reach(177, L), L <= 50, L>= 11, not reach(176, L+1), not reach(466, L+1), not reach(329, L+1).
reach(466,L + 1) :- reach(177, L), L <= 50, L>= 11, not reach(176, L+1), not reach(463, L+1), not reach(329, L+1).
reach(329,L + 1) :- reach(177, L), L <= 50, L>= 11, not reach(176, L+1), not reach(463, L+1), not reach(466, L+1).
reach(177,L + 1) :- reach(463, L), L <= 49, L>= 12, not reach(467, L+1), not reach(420, L+1), not reach(466, L+1).
reach(467,L + 1) :- reach(463, L), L <= 49, L>= 12, not reach(177, L+1), not reach(420, L+1), not reach(466, L+1).
reach(420,L + 1) :- reach(463, L), L <= 49, L>= 12, not reach(177, L+1), not reach(467, L+1), not reach(466, L+1).
reach(466,L + 1) :- reach(463, L), L <= 49, L>= 12, not reach(177, L+1), not reach(467, L+1), not reach(420, L+1).
reach(177,L + 1) :- reach(466, L), L <= 50, L>= 12, not reach(43, L+1), not reach(463, L+1).
reach(43,L + 1) :- reach(466, L), L <= 50, L>= 12, not reach(177, L+1), not reach(463, L+1).
reach(463,L + 1) :- reach(466, L), L <= 50, L>= 12, not reach(177, L+1), not reach(43, L+1).
reach(367,L + 1) :- reach(369, L), L <= 52, L>= 9, not reach(143, L+1), not reach(191, L+1).
reach(143,L + 1) :- reach(369, L), L <= 52, L>= 9, not reach(367, L+1), not reach(191, L+1).
reach(191,L + 1) :- reach(369, L), L <= 52, L>= 9, not reach(367, L+1), not reach(143, L+1).
reach(364,L + 1) :- reach(191, L), L <= 53, L>= 8, not reach(303, L+1), not reach(143, L+1), not reach(369, L+1).
reach(303,L + 1) :- reach(191, L), L <= 53, L>= 8, not reach(364, L+1), not reach(143, L+1), not reach(369, L+1).
reach(143,L + 1) :- reach(191, L), L <= 53, L>= 8, not reach(364, L+1), not reach(303, L+1), not reach(369, L+1).
reach(369,L + 1) :- reach(191, L), L <= 53, L>= 8, not reach(364, L+1), not reach(303, L+1), not reach(143, L+1).
reach(191,L + 1) :- reach(364, L), L <= 54, L>= 7, not reach(367, L+1), not reach(574, L+1), not reach(445, L+1), not reach(298, L+1).
reach(367,L + 1) :- reach(364, L), L <= 54, L>= 7, not reach(191, L+1), not reach(574, L+1), not reach(445, L+1), not reach(298, L+1).
reach(574,L + 1) :- reach(364, L), L <= 54, L>= 7, not reach(191, L+1), not reach(367, L+1), not reach(445, L+1), not reach(298, L+1).
reach(445,L + 1) :- reach(364, L), L <= 54, L>= 7, not reach(191, L+1), not reach(367, L+1), not reach(574, L+1), not reach(298, L+1).
reach(298,L + 1) :- reach(364, L), L <= 54, L>= 7, not reach(191, L+1), not reach(367, L+1), not reach(574, L+1), not reach(445, L+1).
reach(70,L + 1) :- reach(195, L), L <= 53, L>= 4, not reach(456, L+1), not reach(618, L+1), not reach(589, L+1).
reach(456,L + 1) :- reach(195, L), L <= 53, L>= 4, not reach(70, L+1), not reach(618, L+1), not reach(589, L+1).
reach(618,L + 1) :- reach(195, L), L <= 53, L>= 4, not reach(70, L+1), not reach(456, L+1), not reach(589, L+1).
reach(589,L + 1) :- reach(195, L), L <= 53, L>= 4, not reach(70, L+1), not reach(456, L+1), not reach(618, L+1).
reach(457,L + 1) :- reach(589, L), L <= 52, L>= 3, not reach(488, L+1), not reach(195, L+1).
reach(488,L + 1) :- reach(589, L), L <= 52, L>= 3, not reach(457, L+1), not reach(195, L+1).
reach(195,L + 1) :- reach(589, L), L <= 52, L>= 3, not reach(457, L+1), not reach(488, L+1).
reach(283,L + 1) :- reach(456, L), L <= 52, L>= 4, not reach(457, L+1), not reach(195, L+1).
reach(457,L + 1) :- reach(456, L), L <= 52, L>= 4, not reach(283, L+1), not reach(195, L+1).
reach(195,L + 1) :- reach(456, L), L <= 52, L>= 4, not reach(283, L+1), not reach(457, L+1).
reach(494,L + 1) :- reach(618, L), L <= 52, L>= 5, not reach(195, L+1), not reach(283, L+1).
reach(195,L + 1) :- reach(618, L), L <= 52, L>= 5, not reach(494, L+1), not reach(283, L+1).
reach(283,L + 1) :- reach(618, L), L <= 52, L>= 5, not reach(494, L+1), not reach(195, L+1).
reach(445,L + 1) :- reach(209, L), L <= 55, L>= 5, not reach(19, L+1), not reach(493, L+1), not reach(651, L+1), not reach(663, L+1).
reach(19,L + 1) :- reach(209, L), L <= 55, L>= 5, not reach(445, L+1), not reach(493, L+1), not reach(651, L+1), not reach(663, L+1).
reach(493,L + 1) :- reach(209, L), L <= 55, L>= 5, not reach(445, L+1), not reach(19, L+1), not reach(651, L+1), not reach(663, L+1).
reach(651,L + 1) :- reach(209, L), L <= 55, L>= 5, not reach(445, L+1), not reach(19, L+1), not reach(493, L+1), not reach(663, L+1).
reach(663,L + 1) :- reach(209, L), L <= 55, L>= 5, not reach(445, L+1), not reach(19, L+1), not reach(493, L+1), not reach(651, L+1).
reach(209,L + 1) :- reach(445, L), L <= 55, L>= 6, not reach(364, L+1), not reach(274, L+1).
reach(364,L + 1) :- reach(445, L), L <= 55, L>= 6, not reach(209, L+1), not reach(274, L+1).
reach(274,L + 1) :- reach(445, L), L <= 55, L>= 6, not reach(209, L+1), not reach(364, L+1).
reach(214,L + 1) :- reach(703, L), L <= 50, L>= 17, not reach(439, L+1), not reach(704, L+1).
reach(439,L + 1) :- reach(703, L), L <= 50, L>= 17, not reach(214, L+1), not reach(704, L+1).
reach(704,L + 1) :- reach(703, L), L <= 50, L>= 17, not reach(214, L+1), not reach(439, L+1).
reach(228,L + 1) :- reach(224, L), L <= 55, L>= 17, not reach(257, L+1), not reach(227, L+1).
reach(257,L + 1) :- reach(224, L), L <= 55, L>= 17, not reach(228, L+1), not reach(227, L+1).
reach(227,L + 1) :- reach(224, L), L <= 55, L>= 17, not reach(228, L+1), not reach(257, L+1).
reach(224,L + 1) :- reach(228, L), L <= 54, L>= 16, not reach(681, L+1), not reach(233, L+1).
reach(681,L + 1) :- reach(228, L), L <= 54, L>= 16, not reach(224, L+1), not reach(233, L+1).
reach(233,L + 1) :- reach(228, L), L <= 54, L>= 16, not reach(224, L+1), not reach(681, L+1).
reach(224,L + 1) :- reach(257, L), L <= 56, L>= 16, not reach(227, L+1), not reach(497, L+1), not reach(231, L+1).
reach(227,L + 1) :- reach(257, L), L <= 56, L>= 16, not reach(224, L+1), not reach(497, L+1), not reach(231, L+1).
reach(497,L + 1) :- reach(257, L), L <= 56, L>= 16, not reach(224, L+1), not reach(227, L+1), not reach(231, L+1).
reach(231,L + 1) :- reach(257, L), L <= 56, L>= 16, not reach(224, L+1), not reach(227, L+1), not reach(497, L+1).
reach(224,L + 1) :- reach(227, L), L <= 55, L>= 17, not reach(233, L+1), not reach(257, L+1).
reach(233,L + 1) :- reach(227, L), L <= 55, L>= 17, not reach(224, L+1), not reach(257, L+1).
reach(257,L + 1) :- reach(227, L), L <= 55, L>= 17, not reach(224, L+1), not reach(233, L+1).
reach(228,L + 1) :- reach(681, L), L <= 54, L>= 15, not reach(518, L+1), not reach(262, L+1).
reach(518,L + 1) :- reach(681, L), L <= 54, L>= 15, not reach(228, L+1), not reach(262, L+1).
reach(262,L + 1) :- reach(681, L), L <= 54, L>= 15, not reach(228, L+1), not reach(518, L+1).
reach(234,L + 1) :- reach(231, L), L <= 55, L>= 16, not reach(512, L+1), not reach(257, L+1).
reach(512,L + 1) :- reach(231, L), L <= 55, L>= 16, not reach(234, L+1), not reach(257, L+1).
reach(257,L + 1) :- reach(231, L), L <= 55, L>= 16, not reach(234, L+1), not reach(512, L+1).
reach(233,L + 1) :- reach(234, L), L <= 54, L>= 17, not reach(231, L+1), not reach(620, L+1).
reach(231,L + 1) :- reach(234, L), L <= 54, L>= 17, not reach(233, L+1), not reach(620, L+1).
reach(620,L + 1) :- reach(234, L), L <= 54, L>= 17, not reach(233, L+1), not reach(231, L+1).
reach(234,L + 1) :- reach(233, L), L <= 54, L>= 17, not reach(238, L+1), not reach(227, L+1), not reach(262, L+1), not reach(228, L+1).
reach(238,L + 1) :- reach(233, L), L <= 54, L>= 17, not reach(234, L+1), not reach(227, L+1), not reach(262, L+1), not reach(228, L+1).
reach(227,L + 1) :- reach(233, L), L <= 54, L>= 17, not reach(234, L+1), not reach(238, L+1), not reach(262, L+1), not reach(228, L+1).
reach(262,L + 1) :- reach(233, L), L <= 54, L>= 17, not reach(234, L+1), not reach(238, L+1), not reach(227, L+1), not reach(228, L+1).
reach(228,L + 1) :- reach(233, L), L <= 54, L>= 17, not reach(234, L+1), not reach(238, L+1), not reach(227, L+1), not reach(262, L+1).
reach(233,L + 1) :- reach(238, L), L <= 53, L>= 18, not reach(620, L+1), not reach(259, L+1).
reach(620,L + 1) :- reach(238, L), L <= 53, L>= 18, not reach(233, L+1), not reach(259, L+1).
reach(259,L + 1) :- reach(238, L), L <= 53, L>= 18, not reach(233, L+1), not reach(620, L+1).
reach(238,L + 1) :- reach(620, L), L <= 53, L>= 18, not reach(624, L+1), not reach(234, L+1).
reach(624,L + 1) :- reach(620, L), L <= 53, L>= 18, not reach(238, L+1), not reach(234, L+1).
reach(234,L + 1) :- reach(620, L), L <= 53, L>= 18, not reach(238, L+1), not reach(624, L+1).
reach(548,L + 1) :- reach(241, L), L <= 49, L>= 10, not reach(346, L+1), not reach(547, L+1).
reach(346,L + 1) :- reach(241, L), L <= 49, L>= 10, not reach(548, L+1), not reach(547, L+1).
reach(547,L + 1) :- reach(241, L), L <= 49, L>= 10, not reach(548, L+1), not reach(346, L+1).
reach(262,L + 1) :- reach(259, L), L <= 52, L>= 17, not reach(238, L+1), not reach(159, L+1).
reach(238,L + 1) :- reach(259, L), L <= 52, L>= 17, not reach(262, L+1), not reach(159, L+1).
reach(159,L + 1) :- reach(259, L), L <= 52, L>= 17, not reach(262, L+1), not reach(238, L+1).
reach(272,L + 1) :- reach(271, L), L <= 57, L>= 8, not reach(274, L+1), not reach(47, L+1), not reach(530, L+1).
reach(274,L + 1) :- reach(271, L), L <= 57, L>= 8, not reach(272, L+1), not reach(47, L+1), not reach(530, L+1).
reach(47,L + 1) :- reach(271, L), L <= 57, L>= 8, not reach(272, L+1), not reach(274, L+1), not reach(530, L+1).
reach(530,L + 1) :- reach(271, L), L <= 57, L>= 8, not reach(272, L+1), not reach(274, L+1), not reach(47, L+1).
reach(271,L + 1) :- reach(274, L), L <= 56, L>= 7, not reach(438, L+1), not reach(445, L+1).
reach(438,L + 1) :- reach(274, L), L <= 56, L>= 7, not reach(271, L+1), not reach(445, L+1).
reach(445,L + 1) :- reach(274, L), L <= 56, L>= 7, not reach(271, L+1), not reach(438, L+1).
reach(274,L + 1) :- reach(438, L), L <= 55, L>= 8, not reach(435, L+1), not reach(306, L+1).
reach(435,L + 1) :- reach(438, L), L <= 55, L>= 8, not reach(274, L+1), not reach(306, L+1).
reach(306,L + 1) :- reach(438, L), L <= 55, L>= 8, not reach(274, L+1), not reach(435, L+1).
reach(678,L + 1) :- reach(698, L), L <= 58, L>= 11, not reach(145, L+1), not reach(563, L+1), not reach(395, L+1).
reach(145,L + 1) :- reach(698, L), L <= 58, L>= 11, not reach(678, L+1), not reach(563, L+1), not reach(395, L+1).
reach(563,L + 1) :- reach(698, L), L <= 58, L>= 11, not reach(678, L+1), not reach(145, L+1), not reach(395, L+1).
reach(395,L + 1) :- reach(698, L), L <= 58, L>= 11, not reach(678, L+1), not reach(145, L+1), not reach(563, L+1).
reach(280,L + 1) :- reach(278, L), L <= 50, L>= 4, not reach(427, L+1), not reach(490, L+1).
reach(427,L + 1) :- reach(278, L), L <= 50, L>= 4, not reach(280, L+1), not reach(490, L+1).
reach(490,L + 1) :- reach(278, L), L <= 50, L>= 4, not reach(280, L+1), not reach(427, L+1).
reach(278,L + 1) :- reach(280, L), L <= 50, L>= 4, not reach(487, L+1), not reach(494, L+1).
reach(487,L + 1) :- reach(280, L), L <= 50, L>= 4, not reach(278, L+1), not reach(494, L+1).
reach(494,L + 1) :- reach(280, L), L <= 50, L>= 4, not reach(278, L+1), not reach(487, L+1).
reach(487,L + 1) :- reach(427, L), L <= 50, L>= 3, not reach(278, L+1), not reach(335, L+1).
reach(278,L + 1) :- reach(427, L), L <= 50, L>= 3, not reach(487, L+1), not reach(335, L+1).
reach(335,L + 1) :- reach(427, L), L <= 50, L>= 3, not reach(487, L+1), not reach(278, L+1).
reach(489,L + 1) :- reach(490, L), L <= 51, L>= 5, not reach(494, L+1), not reach(278, L+1).
reach(494,L + 1) :- reach(490, L), L <= 51, L>= 5, not reach(489, L+1), not reach(278, L+1).
reach(278,L + 1) :- reach(490, L), L <= 51, L>= 5, not reach(489, L+1), not reach(494, L+1).
reach(456,L + 1) :- reach(283, L), L <= 51, L>= 4, not reach(457, L+1), not reach(618, L+1).
reach(457,L + 1) :- reach(283, L), L <= 51, L>= 4, not reach(456, L+1), not reach(618, L+1).
reach(618,L + 1) :- reach(283, L), L <= 51, L>= 4, not reach(456, L+1), not reach(457, L+1).
reach(283,L + 1) :- reach(457, L), L <= 51, L>= 3, not reach(456, L+1), not reach(589, L+1), not reach(488, L+1).
reach(456,L + 1) :- reach(457, L), L <= 51, L>= 3, not reach(283, L+1), not reach(589, L+1), not reach(488, L+1).
reach(589,L + 1) :- reach(457, L), L <= 51, L>= 3, not reach(283, L+1), not reach(456, L+1), not reach(488, L+1).
reach(488,L + 1) :- reach(457, L), L <= 51, L>= 3, not reach(283, L+1), not reach(456, L+1), not reach(589, L+1).
reach(288,L + 1) :- reach(285, L), L <= 50, L>= 14, not reach(289, L+1), not reach(290, L+1), not reach(128, L+1).
reach(289,L + 1) :- reach(285, L), L <= 50, L>= 14, not reach(288, L+1), not reach(290, L+1), not reach(128, L+1).
reach(290,L + 1) :- reach(285, L), L <= 50, L>= 14, not reach(288, L+1), not reach(289, L+1), not reach(128, L+1).
reach(128,L + 1) :- reach(285, L), L <= 50, L>= 14, not reach(288, L+1), not reach(289, L+1), not reach(290, L+1).
reach(285,L + 1) :- reach(288, L), L <= 49, L>= 13, not reach(291, L+1), not reach(289, L+1).
reach(291,L + 1) :- reach(288, L), L <= 49, L>= 13, not reach(285, L+1), not reach(289, L+1).
reach(289,L + 1) :- reach(288, L), L <= 49, L>= 13, not reach(285, L+1), not reach(291, L+1).
reach(285,L + 1) :- reach(289, L), L <= 50, L>= 14, not reach(524, L+1), not reach(288, L+1), not reach(528, L+1).
reach(524,L + 1) :- reach(289, L), L <= 50, L>= 14, not reach(285, L+1), not reach(288, L+1), not reach(528, L+1).
reach(288,L + 1) :- reach(289, L), L <= 50, L>= 14, not reach(285, L+1), not reach(524, L+1), not reach(528, L+1).
reach(528,L + 1) :- reach(289, L), L <= 50, L>= 14, not reach(285, L+1), not reach(524, L+1), not reach(288, L+1).
reach(291,L + 1) :- reach(290, L), L <= 49, L>= 13, not reach(292, L+1), not reach(285, L+1).
reach(292,L + 1) :- reach(290, L), L <= 49, L>= 13, not reach(291, L+1), not reach(285, L+1).
reach(285,L + 1) :- reach(290, L), L <= 49, L>= 13, not reach(291, L+1), not reach(292, L+1).
reach(288,L + 1) :- reach(291, L), L <= 49, L>= 12, not reach(290, L+1), not reach(327, L+1).
reach(290,L + 1) :- reach(291, L), L <= 49, L>= 12, not reach(288, L+1), not reach(327, L+1).
reach(327,L + 1) :- reach(291, L), L <= 49, L>= 12, not reach(288, L+1), not reach(290, L+1).
reach(289,L + 1) :- reach(524, L), L <= 51, L>= 15, not reach(528, L+1), not reach(516, L+1).
reach(528,L + 1) :- reach(524, L), L <= 51, L>= 15, not reach(289, L+1), not reach(516, L+1).
reach(516,L + 1) :- reach(524, L), L <= 51, L>= 15, not reach(289, L+1), not reach(528, L+1).
reach(290,L + 1) :- reach(292, L), L <= 50, L>= 13, not reach(127, L+1), not reach(688, L+1).
reach(127,L + 1) :- reach(292, L), L <= 50, L>= 13, not reach(290, L+1), not reach(688, L+1).
reach(688,L + 1) :- reach(292, L), L <= 50, L>= 13, not reach(290, L+1), not reach(127, L+1).
reach(302,L + 1) :- reach(298, L), L <= 53, L>= 8, not reach(303, L+1), not reach(300, L+1), not reach(364, L+1).
reach(303,L + 1) :- reach(298, L), L <= 53, L>= 8, not reach(302, L+1), not reach(300, L+1), not reach(364, L+1).
reach(300,L + 1) :- reach(298, L), L <= 53, L>= 8, not reach(302, L+1), not reach(303, L+1), not reach(364, L+1).
reach(364,L + 1) :- reach(298, L), L <= 53, L>= 8, not reach(302, L+1), not reach(303, L+1), not reach(300, L+1).
reach(310,L + 1) :- reach(300, L), L <= 52, L>= 9, not reach(302, L+1), not reach(298, L+1).
reach(302,L + 1) :- reach(300, L), L <= 52, L>= 9, not reach(310, L+1), not reach(298, L+1).
reach(298,L + 1) :- reach(300, L), L <= 52, L>= 9, not reach(310, L+1), not reach(302, L+1).
reach(298,L + 1) :- reach(302, L), L <= 52, L>= 9, not reach(98, L+1), not reach(300, L+1).
reach(98,L + 1) :- reach(302, L), L <= 52, L>= 9, not reach(298, L+1), not reach(300, L+1).
reach(300,L + 1) :- reach(302, L), L <= 52, L>= 9, not reach(298, L+1), not reach(98, L+1).
reach(300,L + 1) :- reach(310, L), L <= 51, L>= 10, not reach(107, L+1), not reach(307, L+1).
reach(107,L + 1) :- reach(310, L), L <= 51, L>= 10, not reach(300, L+1), not reach(307, L+1).
reach(307,L + 1) :- reach(310, L), L <= 51, L>= 10, not reach(300, L+1), not reach(107, L+1).
reach(438,L + 1) :- reach(306, L), L <= 55, L>= 9, not reach(530, L+1), not reach(435, L+1).
reach(530,L + 1) :- reach(306, L), L <= 55, L>= 9, not reach(438, L+1), not reach(435, L+1).
reach(435,L + 1) :- reach(306, L), L <= 55, L>= 9, not reach(438, L+1), not reach(530, L+1).
reach(314,L + 1) :- reach(315, L), L <= 51, L>= 7, not reach(352, L+1), not reach(318, L+1).
reach(352,L + 1) :- reach(315, L), L <= 51, L>= 7, not reach(314, L+1), not reach(318, L+1).
reach(318,L + 1) :- reach(315, L), L <= 51, L>= 7, not reach(314, L+1), not reach(352, L+1).
reach(315,L + 1) :- reach(318, L), L <= 52, L>= 8, not reach(370, L+1), not reach(58, L+1).
reach(370,L + 1) :- reach(318, L), L <= 52, L>= 8, not reach(315, L+1), not reach(58, L+1).
reach(58,L + 1) :- reach(318, L), L <= 52, L>= 8, not reach(315, L+1), not reach(370, L+1).
reach(318,L + 1) :- reach(370, L), L <= 53, L>= 7, not reach(747, L+1), not reach(367, L+1), not reach(651, L+1), not reach(746, L+1).
reach(747,L + 1) :- reach(370, L), L <= 53, L>= 7, not reach(318, L+1), not reach(367, L+1), not reach(651, L+1), not reach(746, L+1).
reach(367,L + 1) :- reach(370, L), L <= 53, L>= 7, not reach(318, L+1), not reach(747, L+1), not reach(651, L+1), not reach(746, L+1).
reach(651,L + 1) :- reach(370, L), L <= 53, L>= 7, not reach(318, L+1), not reach(747, L+1), not reach(367, L+1), not reach(746, L+1).
reach(746,L + 1) :- reach(370, L), L <= 53, L>= 7, not reach(318, L+1), not reach(747, L+1), not reach(367, L+1), not reach(651, L+1).
reach(467,L + 1) :- reach(325, L), L <= 49, L>= 12, not reach(327, L+1), not reach(420, L+1).
reach(327,L + 1) :- reach(325, L), L <= 49, L>= 12, not reach(467, L+1), not reach(420, L+1).
reach(420,L + 1) :- reach(325, L), L <= 49, L>= 12, not reach(467, L+1), not reach(327, L+1).
reach(325,L + 1) :- reach(467, L), L <= 48, L>= 13, not reach(463, L+1), not reach(420, L+1).
reach(463,L + 1) :- reach(467, L), L <= 48, L>= 13, not reach(325, L+1), not reach(420, L+1).
reach(420,L + 1) :- reach(467, L), L <= 48, L>= 13, not reach(325, L+1), not reach(463, L+1).
reach(327,L + 1) :- reach(329, L), L <= 50, L>= 11, not reach(250, L+1), not reach(177, L+1).
reach(250,L + 1) :- reach(329, L), L <= 50, L>= 11, not reach(327, L+1), not reach(177, L+1).
reach(177,L + 1) :- reach(329, L), L <= 50, L>= 11, not reach(327, L+1), not reach(250, L+1).
reach(14,L + 1) :- reach(335, L), L <= 51, L>= 2, not reach(34, L+1), not reach(427, L+1).
reach(34,L + 1) :- reach(335, L), L <= 51, L>= 2, not reach(14, L+1), not reach(427, L+1).
reach(427,L + 1) :- reach(335, L), L <= 51, L>= 2, not reach(14, L+1), not reach(34, L+1).
reach(241,L + 1) :- reach(548, L), L <= 49, L>= 11, not reach(688, L+1), not reach(547, L+1).
reach(688,L + 1) :- reach(548, L), L <= 49, L>= 11, not reach(241, L+1), not reach(547, L+1).
reach(547,L + 1) :- reach(548, L), L <= 49, L>= 11, not reach(241, L+1), not reach(688, L+1).
reach(241,L + 1) :- reach(346, L), L <= 49, L>= 9, not reach(355, L+1), not reach(542, L+1).
reach(355,L + 1) :- reach(346, L), L <= 49, L>= 9, not reach(241, L+1), not reach(542, L+1).
reach(542,L + 1) :- reach(346, L), L <= 49, L>= 9, not reach(241, L+1), not reach(355, L+1).
reach(250,L + 1) :- reach(355, L), L <= 50, L>= 9, not reach(346, L+1), not reach(542, L+1).
reach(346,L + 1) :- reach(355, L), L <= 50, L>= 9, not reach(250, L+1), not reach(542, L+1).
reach(542,L + 1) :- reach(355, L), L <= 50, L>= 9, not reach(250, L+1), not reach(346, L+1).
reach(19,L + 1) :- reach(359, L), L <= 53, L>= 3, not reach(657, L+1), not reach(14, L+1).
reach(657,L + 1) :- reach(359, L), L <= 53, L>= 3, not reach(19, L+1), not reach(14, L+1).
reach(14,L + 1) :- reach(359, L), L <= 53, L>= 3, not reach(19, L+1), not reach(657, L+1).
reach(364,L + 1) :- reach(367, L), L <= 53, L>= 8, not reach(369, L+1), not reach(370, L+1).
reach(369,L + 1) :- reach(367, L), L <= 53, L>= 8, not reach(364, L+1), not reach(370, L+1).
reach(370,L + 1) :- reach(367, L), L <= 53, L>= 8, not reach(364, L+1), not reach(369, L+1).
reach(370,L + 1) :- reach(747, L), L <= 52, L>= 7, not reach(746, L+1), not reach(574, L+1).
reach(746,L + 1) :- reach(747, L), L <= 52, L>= 7, not reach(370, L+1), not reach(574, L+1).
reach(574,L + 1) :- reach(747, L), L <= 52, L>= 7, not reach(370, L+1), not reach(746, L+1).
reach(28,L + 1) :- reach(383, L), L <= 60, L>= 10, not reach(678, L+1), not reach(501, L+1), not reach(564, L+1).
reach(678,L + 1) :- reach(383, L), L <= 60, L>= 10, not reach(28, L+1), not reach(501, L+1), not reach(564, L+1).
reach(501,L + 1) :- reach(383, L), L <= 60, L>= 10, not reach(28, L+1), not reach(678, L+1), not reach(564, L+1).
reach(564,L + 1) :- reach(383, L), L <= 60, L>= 10, not reach(28, L+1), not reach(678, L+1), not reach(501, L+1).
reach(390,L + 1) :- reach(395, L), L <= 57, L>= 12, not reach(397, L+1), not reach(387, L+1), not reach(560, L+1), not reach(392, L+1), not reach(698, L+1).
reach(397,L + 1) :- reach(395, L), L <= 57, L>= 12, not reach(390, L+1), not reach(387, L+1), not reach(560, L+1), not reach(392, L+1), not reach(698, L+1).
reach(387,L + 1) :- reach(395, L), L <= 57, L>= 12, not reach(390, L+1), not reach(397, L+1), not reach(560, L+1), not reach(392, L+1), not reach(698, L+1).
reach(560,L + 1) :- reach(395, L), L <= 57, L>= 12, not reach(390, L+1), not reach(397, L+1), not reach(387, L+1), not reach(392, L+1), not reach(698, L+1).
reach(392,L + 1) :- reach(395, L), L <= 57, L>= 12, not reach(390, L+1), not reach(397, L+1), not reach(387, L+1), not reach(560, L+1), not reach(698, L+1).
reach(698,L + 1) :- reach(395, L), L <= 57, L>= 12, not reach(390, L+1), not reach(397, L+1), not reach(387, L+1), not reach(560, L+1), not reach(392, L+1).
reach(395,L + 1) :- reach(397, L), L <= 56, L>= 13, not reach(552, L+1), not reach(390, L+1).
reach(552,L + 1) :- reach(397, L), L <= 56, L>= 13, not reach(395, L+1), not reach(390, L+1).
reach(390,L + 1) :- reach(397, L), L <= 56, L>= 13, not reach(395, L+1), not reach(552, L+1).
reach(407,L + 1) :- reach(404, L), L <= 53, L>= 5, not reach(19, L+1), not reach(348, L+1).
reach(19,L + 1) :- reach(404, L), L <= 53, L>= 5, not reach(407, L+1), not reach(348, L+1).
reach(348,L + 1) :- reach(404, L), L <= 53, L>= 5, not reach(407, L+1), not reach(19, L+1).
reach(404,L + 1) :- reach(407, L), L <= 53, L>= 6, not reach(574, L+1), not reach(651, L+1).
reach(574,L + 1) :- reach(407, L), L <= 53, L>= 6, not reach(404, L+1), not reach(651, L+1).
reach(651,L + 1) :- reach(407, L), L <= 53, L>= 6, not reach(404, L+1), not reach(574, L+1).
reach(407,L + 1) :- reach(574, L), L <= 53, L>= 7, not reach(364, L+1), not reach(747, L+1).
reach(364,L + 1) :- reach(574, L), L <= 53, L>= 7, not reach(407, L+1), not reach(747, L+1).
reach(747,L + 1) :- reach(574, L), L <= 53, L>= 7, not reach(407, L+1), not reach(364, L+1).
reach(407,L + 1) :- reach(651, L), L <= 54, L>= 6, not reach(209, L+1), not reach(370, L+1).
reach(209,L + 1) :- reach(651, L), L <= 54, L>= 6, not reach(407, L+1), not reach(370, L+1).
reach(370,L + 1) :- reach(651, L), L <= 54, L>= 6, not reach(407, L+1), not reach(209, L+1).
reach(528,L + 1) :- reach(420, L), L <= 49, L>= 13, not reach(463, L+1), not reach(467, L+1), not reach(325, L+1).
reach(463,L + 1) :- reach(420, L), L <= 49, L>= 13, not reach(528, L+1), not reach(467, L+1), not reach(325, L+1).
reach(467,L + 1) :- reach(420, L), L <= 49, L>= 13, not reach(528, L+1), not reach(463, L+1), not reach(325, L+1).
reach(325,L + 1) :- reach(420, L), L <= 49, L>= 13, not reach(528, L+1), not reach(463, L+1), not reach(467, L+1).
reach(420,L + 1) :- reach(528, L), L <= 50, L>= 14, not reach(524, L+1), not reach(289, L+1).
reach(524,L + 1) :- reach(528, L), L <= 50, L>= 14, not reach(420, L+1), not reach(289, L+1).
reach(289,L + 1) :- reach(528, L), L <= 50, L>= 14, not reach(420, L+1), not reach(524, L+1).
reach(490,L + 1) :- reach(494, L), L <= 51, L>= 5, not reach(618, L+1), not reach(280, L+1).
reach(618,L + 1) :- reach(494, L), L <= 51, L>= 5, not reach(490, L+1), not reach(280, L+1).
reach(280,L + 1) :- reach(494, L), L <= 51, L>= 5, not reach(490, L+1), not reach(618, L+1).
reach(439,L + 1) :- reach(434, L), L <= 52, L>= 16, not reach(86, L+1), not reach(153, L+1).
reach(86,L + 1) :- reach(434, L), L <= 52, L>= 16, not reach(439, L+1), not reach(153, L+1).
reach(153,L + 1) :- reach(434, L), L <= 52, L>= 16, not reach(439, L+1), not reach(86, L+1).
reach(434,L + 1) :- reach(439, L), L <= 51, L>= 17, not reach(515, L+1), not reach(703, L+1).
reach(515,L + 1) :- reach(439, L), L <= 51, L>= 17, not reach(434, L+1), not reach(703, L+1).
reach(703,L + 1) :- reach(439, L), L <= 51, L>= 17, not reach(434, L+1), not reach(515, L+1).
reach(438,L + 1) :- reach(435, L), L <= 55, L>= 9, not reach(306, L+1), not reach(530, L+1).
reach(306,L + 1) :- reach(435, L), L <= 55, L>= 9, not reach(438, L+1), not reach(530, L+1).
reach(530,L + 1) :- reach(435, L), L <= 55, L>= 9, not reach(438, L+1), not reach(306, L+1).
reach(306,L + 1) :- reach(530, L), L <= 56, L>= 9, not reach(435, L+1), not reach(271, L+1).
reach(435,L + 1) :- reach(530, L), L <= 56, L>= 9, not reach(306, L+1), not reach(271, L+1).
reach(271,L + 1) :- reach(530, L), L <= 56, L>= 9, not reach(306, L+1), not reach(435, L+1).
reach(439,L + 1) :- reach(515, L), L <= 50, L>= 18, not reach(452, L+1), not reach(453, L+1), not reach(619, L+1).
reach(452,L + 1) :- reach(515, L), L <= 50, L>= 18, not reach(439, L+1), not reach(453, L+1), not reach(619, L+1).
reach(453,L + 1) :- reach(515, L), L <= 50, L>= 18, not reach(439, L+1), not reach(452, L+1), not reach(619, L+1).
reach(619,L + 1) :- reach(515, L), L <= 50, L>= 18, not reach(439, L+1), not reach(452, L+1), not reach(453, L+1).
reach(619,L + 1) :- reach(453, L), L <= 51, L>= 18, not reach(515, L+1), not reach(452, L+1), not reach(616, L+1).
reach(515,L + 1) :- reach(453, L), L <= 51, L>= 18, not reach(619, L+1), not reach(452, L+1), not reach(616, L+1).
reach(452,L + 1) :- reach(453, L), L <= 51, L>= 18, not reach(619, L+1), not reach(515, L+1), not reach(616, L+1).
reach(616,L + 1) :- reach(453, L), L <= 51, L>= 18, not reach(619, L+1), not reach(515, L+1), not reach(452, L+1).
reach(453,L + 1) :- reach(619, L), L <= 50, L>= 19, not reach(515, L+1), not reach(613, L+1).
reach(515,L + 1) :- reach(619, L), L <= 50, L>= 19, not reach(453, L+1), not reach(613, L+1).
reach(613,L + 1) :- reach(619, L), L <= 50, L>= 19, not reach(453, L+1), not reach(515, L+1).
reach(487,L + 1) :- reach(489, L), L <= 52, L>= 4, not reach(490, L+1), not reach(92, L+1).
reach(490,L + 1) :- reach(489, L), L <= 52, L>= 4, not reach(487, L+1), not reach(92, L+1).
reach(92,L + 1) :- reach(489, L), L <= 52, L>= 4, not reach(487, L+1), not reach(490, L+1).
reach(457,L + 1) :- reach(488, L), L <= 51, L>= 2, not reach(36, L+1), not reach(589, L+1).
reach(36,L + 1) :- reach(488, L), L <= 51, L>= 2, not reach(457, L+1), not reach(589, L+1).
reach(589,L + 1) :- reach(488, L), L <= 51, L>= 2, not reach(457, L+1), not reach(36, L+1).
reach(663,L + 1) :- reach(493, L), L <= 56, L>= 6, not reach(47, L+1), not reach(209, L+1).
reach(47,L + 1) :- reach(493, L), L <= 56, L>= 6, not reach(663, L+1), not reach(209, L+1).
reach(209,L + 1) :- reach(493, L), L <= 56, L>= 6, not reach(663, L+1), not reach(47, L+1).
reach(151,L + 1) :- reach(497, L), L <= 57, L>= 15, not reach(257, L+1), not reach(674, L+1).
reach(257,L + 1) :- reach(497, L), L <= 57, L>= 15, not reach(151, L+1), not reach(674, L+1).
reach(674,L + 1) :- reach(497, L), L <= 57, L>= 15, not reach(151, L+1), not reach(257, L+1).
reach(383,L + 1) :- reach(501, L), L <= 59, L>= 11, not reach(117, L+1), not reach(505, L+1).
reach(117,L + 1) :- reach(501, L), L <= 59, L>= 11, not reach(383, L+1), not reach(505, L+1).
reach(505,L + 1) :- reach(501, L), L <= 59, L>= 11, not reach(383, L+1), not reach(117, L+1).
reach(145,L + 1) :- reach(505, L), L <= 58, L>= 12, not reach(501, L+1), not reach(672, L+1).
reach(501,L + 1) :- reach(505, L), L <= 58, L>= 12, not reach(145, L+1), not reach(672, L+1).
reach(672,L + 1) :- reach(505, L), L <= 58, L>= 12, not reach(145, L+1), not reach(501, L+1).
reach(231,L + 1) :- reach(512, L), L <= 54, L>= 15, not reach(632, L+1), not reach(630, L+1).
reach(632,L + 1) :- reach(512, L), L <= 54, L>= 15, not reach(231, L+1), not reach(630, L+1).
reach(630,L + 1) :- reach(512, L), L <= 54, L>= 15, not reach(231, L+1), not reach(632, L+1).
reach(681,L + 1) :- reach(518, L), L <= 55, L>= 14, not reach(151, L+1), not reach(560, L+1), not reach(737, L+1).
reach(151,L + 1) :- reach(518, L), L <= 55, L>= 14, not reach(681, L+1), not reach(560, L+1), not reach(737, L+1).
reach(560,L + 1) :- reach(518, L), L <= 55, L>= 14, not reach(681, L+1), not reach(151, L+1), not reach(737, L+1).
reach(737,L + 1) :- reach(518, L), L <= 55, L>= 14, not reach(681, L+1), not reach(151, L+1), not reach(560, L+1).
reach(217,L + 1) :- reach(630, L), L <= 53, L>= 14, not reach(516, L+1), not reach(512, L+1).
reach(516,L + 1) :- reach(630, L), L <= 53, L>= 14, not reach(217, L+1), not reach(512, L+1).
reach(512,L + 1) :- reach(630, L), L <= 53, L>= 14, not reach(217, L+1), not reach(516, L+1).
reach(57,L + 1) :- reach(542, L), L <= 50, L>= 8, not reach(249, L+1), not reach(346, L+1), not reach(355, L+1).
reach(249,L + 1) :- reach(542, L), L <= 50, L>= 8, not reach(57, L+1), not reach(346, L+1), not reach(355, L+1).
reach(346,L + 1) :- reach(542, L), L <= 50, L>= 8, not reach(57, L+1), not reach(249, L+1), not reach(355, L+1).
reach(355,L + 1) :- reach(542, L), L <= 50, L>= 8, not reach(57, L+1), not reach(249, L+1), not reach(346, L+1).
reach(547,L + 1) :- reach(688, L), L <= 49, L>= 12, not reach(548, L+1), not reach(292, L+1).
reach(548,L + 1) :- reach(688, L), L <= 49, L>= 12, not reach(547, L+1), not reach(292, L+1).
reach(292,L + 1) :- reach(688, L), L <= 49, L>= 12, not reach(547, L+1), not reach(548, L+1).
reach(747,L + 1) :- reach(746, L), L <= 52, L>= 6, not reach(58, L+1), not reach(370, L+1), not reach(348, L+1).
reach(58,L + 1) :- reach(746, L), L <= 52, L>= 6, not reach(747, L+1), not reach(370, L+1), not reach(348, L+1).
reach(370,L + 1) :- reach(746, L), L <= 52, L>= 6, not reach(747, L+1), not reach(58, L+1), not reach(348, L+1).
reach(348,L + 1) :- reach(746, L), L <= 52, L>= 6, not reach(747, L+1), not reach(58, L+1), not reach(370, L+1).
reach(678,L + 1) :- reach(563, L), L <= 58, L>= 10, not reach(564, L+1), not reach(698, L+1).
reach(564,L + 1) :- reach(563, L), L <= 58, L>= 10, not reach(678, L+1), not reach(698, L+1).
reach(698,L + 1) :- reach(563, L), L <= 58, L>= 10, not reach(678, L+1), not reach(564, L+1).
reach(563,L + 1) :- reach(564, L), L <= 59, L>= 9, not reach(383, L+1), not reach(272, L+1).
reach(383,L + 1) :- reach(564, L), L <= 59, L>= 9, not reach(563, L+1), not reach(272, L+1).
reach(272,L + 1) :- reach(564, L), L <= 59, L>= 9, not reach(563, L+1), not reach(383, L+1).
reach(254,L + 1) :- reach(674, L), L <= 58, L>= 14, not reach(497, L+1), not reach(672, L+1).
reach(497,L + 1) :- reach(674, L), L <= 58, L>= 14, not reach(254, L+1), not reach(672, L+1).
reach(672,L + 1) :- reach(674, L), L <= 58, L>= 14, not reach(254, L+1), not reach(497, L+1).
reach(616,L + 1) :- reach(613, L), L <= 51, L>= 18, not reach(619, L+1), not reach(642, L+1).
reach(619,L + 1) :- reach(613, L), L <= 51, L>= 18, not reach(616, L+1), not reach(642, L+1).
reach(642,L + 1) :- reach(613, L), L <= 51, L>= 18, not reach(616, L+1), not reach(619, L+1).
reach(613,L + 1) :- reach(616, L), L <= 52, L>= 17, not reach(632, L+1), not reach(453, L+1).
reach(632,L + 1) :- reach(616, L), L <= 52, L>= 17, not reach(613, L+1), not reach(453, L+1).
reach(453,L + 1) :- reach(616, L), L <= 52, L>= 17, not reach(613, L+1), not reach(632, L+1).
reach(624,L + 1) :- reach(642, L), L <= 51, L>= 18, not reach(635, L+1), not reach(613, L+1).
reach(635,L + 1) :- reach(642, L), L <= 51, L>= 18, not reach(624, L+1), not reach(613, L+1).
reach(613,L + 1) :- reach(642, L), L <= 51, L>= 18, not reach(624, L+1), not reach(635, L+1).
reach(620,L + 1) :- reach(624, L), L <= 52, L>= 18, not reach(642, L+1), not reach(635, L+1).
reach(642,L + 1) :- reach(624, L), L <= 52, L>= 18, not reach(620, L+1), not reach(635, L+1).
reach(635,L + 1) :- reach(624, L), L <= 52, L>= 18, not reach(620, L+1), not reach(642, L+1).
reach(616,L + 1) :- reach(632, L), L <= 53, L>= 16, not reach(512, L+1), not reach(635, L+1).
reach(512,L + 1) :- reach(632, L), L <= 53, L>= 16, not reach(616, L+1), not reach(635, L+1).
reach(635,L + 1) :- reach(632, L), L <= 53, L>= 16, not reach(616, L+1), not reach(512, L+1).
reach(642,L + 1) :- reach(635, L), L <= 52, L>= 17, not reach(624, L+1), not reach(632, L+1).
reach(624,L + 1) :- reach(635, L), L <= 52, L>= 17, not reach(642, L+1), not reach(632, L+1).
reach(632,L + 1) :- reach(635, L), L <= 52, L>= 17, not reach(642, L+1), not reach(624, L+1).
reach(254,L + 1) :- reach(672, L), L <= 58, L>= 13, not reach(505, L+1), not reach(674, L+1).
reach(505,L + 1) :- reach(672, L), L <= 58, L>= 13, not reach(254, L+1), not reach(674, L+1).
reach(674,L + 1) :- reach(672, L), L <= 58, L>= 13, not reach(254, L+1), not reach(505, L+1).

