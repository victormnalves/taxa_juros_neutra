function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 516);

T = dynare_ss_1970_alt_cobb.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(234) = y(114)/params(105)/T(172);
T(235) = y(115)/params(106)/T(173);
T(236) = y(116)/params(107)/T(174);
T(237) = y(117)/params(108)/T(175);
T(238) = y(118)/params(109)/T(176);
T(239) = y(119)/params(110)/T(177);
T(240) = y(120)/params(111)/T(178);
T(241) = y(121)/params(112)/T(179);
T(242) = y(122)/params(113)/T(180);
T(243) = y(123)/params(114)/T(181);
T(244) = y(124)/params(115)/T(182);
T(245) = y(125)/params(116)/T(183);
T(246) = y(126)/params(117)/T(184);
T(247) = y(127)/params(118)/T(185);
T(248) = y(128)/params(119)/T(186);
T(249) = y(129)/params(120)/T(187);
T(250) = y(130)/params(121)/T(188);
T(251) = y(131)/params(122)/T(189);
T(252) = y(132)/params(123)/T(190);
T(253) = y(133)/params(124)/T(191);
T(254) = y(134)/params(125)/T(192);
T(255) = y(135)/params(126)/T(193);
T(256) = y(136)/params(127)/T(194);
T(257) = y(137)/params(128)/T(195);
T(258) = y(138)/params(129)/T(196);
T(259) = y(139)/params(130)/T(197);
T(260) = y(140)/params(131)/T(198);
T(261) = y(141)/params(132)/T(199);
T(262) = y(142)/params(133)/T(200);
T(263) = y(143)/params(134)/T(201);
T(264) = y(144)/params(135)/T(202);
T(265) = y(145)/params(136)/T(204);
T(266) = y(146)/params(137)/T(205);
T(267) = y(147)/params(138)/T(206);
T(268) = y(148)/params(139)/T(207);
T(269) = y(149)/params(140)/T(208);
T(270) = y(150)/params(141)/T(209);
T(271) = y(151)/params(142)/T(210);
T(272) = y(152)/params(143)/T(211);
T(273) = y(153)/params(144)/T(216);
T(274) = y(154)/params(145)/T(217);
T(275) = y(155)/params(146)/T(218);
T(276) = y(156)/params(147)/T(219);
T(277) = y(157)/params(148)/T(220);
T(278) = y(158)/params(149)/T(221);
T(279) = y(159)/params(150)/T(222);
T(280) = y(160)/params(151)/T(223);
T(281) = y(161)/params(152)/T(224);
T(282) = y(162)/params(153)/T(225);
T(283) = y(163)/params(154)/T(226);
T(284) = y(164)/params(155)/T(227);
T(285) = y(165)/params(156)/T(228);
T(286) = y(166)/params(157)/T(229);
T(287) = y(167)/params(158)/T(230);
T(288) = y(168)/params(159)/T(231);
T(289) = (-y(58))/(y(57)*y(57));
T(290) = getPowerDeriv(y(58)/y(57),T(3),1);
T(291) = getPowerDeriv(y(57),T(6),1);
T(292) = (-y(59))/(y(58)*y(58));
T(293) = getPowerDeriv(y(59)/y(58),T(3),1);
T(294) = getPowerDeriv(y(58),T(6),1);
T(295) = (-y(60))/(y(59)*y(59));
T(296) = getPowerDeriv(y(60)/y(59),T(3),1);
T(297) = getPowerDeriv(y(59),T(6),1);
T(298) = (-y(61))/(y(60)*y(60));
T(299) = getPowerDeriv(y(61)/y(60),T(3),1);
T(300) = getPowerDeriv(y(60),T(6),1);
T(301) = (-y(62))/(y(61)*y(61));
T(302) = getPowerDeriv(y(62)/y(61),T(3),1);
T(303) = getPowerDeriv(y(61),T(6),1);
T(304) = (-y(63))/(y(62)*y(62));
T(305) = getPowerDeriv(y(63)/y(62),T(3),1);
T(306) = getPowerDeriv(y(62),T(6),1);
T(307) = (-y(64))/(y(63)*y(63));
T(308) = getPowerDeriv(y(64)/y(63),T(3),1);
T(309) = getPowerDeriv(y(63),T(6),1);
T(310) = (-y(65))/(y(64)*y(64));
T(311) = getPowerDeriv(y(65)/y(64),T(3),1);
T(312) = getPowerDeriv(y(64),T(6),1);
T(313) = (-y(66))/(y(65)*y(65));
T(314) = getPowerDeriv(y(66)/y(65),T(3),1);
T(315) = getPowerDeriv(y(65),T(6),1);
T(316) = (-y(67))/(y(66)*y(66));
T(317) = getPowerDeriv(y(67)/y(66),T(3),1);
T(318) = getPowerDeriv(y(66),T(6),1);
T(319) = (-y(68))/(y(67)*y(67));
T(320) = getPowerDeriv(y(68)/y(67),T(3),1);
T(321) = getPowerDeriv(y(67),T(6),1);
T(322) = (-y(69))/(y(68)*y(68));
T(323) = getPowerDeriv(y(69)/y(68),T(3),1);
T(324) = getPowerDeriv(y(68),T(6),1);
T(325) = (-y(70))/(y(69)*y(69));
T(326) = getPowerDeriv(y(70)/y(69),T(3),1);
T(327) = getPowerDeriv(y(69),T(6),1);
T(328) = (-y(71))/(y(70)*y(70));
T(329) = getPowerDeriv(y(71)/y(70),T(3),1);
T(330) = getPowerDeriv(y(70),T(6),1);
T(331) = (-y(72))/(y(71)*y(71));
T(332) = getPowerDeriv(y(72)/y(71),T(3),1);
T(333) = getPowerDeriv(y(71),T(6),1);
T(334) = (-y(73))/(y(72)*y(72));
T(335) = getPowerDeriv(y(73)/y(72),T(3),1);
T(336) = getPowerDeriv(y(72),T(6),1);
T(337) = (-y(74))/(y(73)*y(73));
T(338) = getPowerDeriv(y(74)/y(73),T(3),1);
T(339) = getPowerDeriv(y(73),T(6),1);
T(340) = (-y(75))/(y(74)*y(74));
T(341) = getPowerDeriv(y(75)/y(74),T(3),1);
T(342) = getPowerDeriv(y(74),T(6),1);
T(343) = (-y(76))/(y(75)*y(75));
T(344) = getPowerDeriv(y(76)/y(75),T(3),1);
T(345) = getPowerDeriv(y(75),T(6),1);
T(346) = (-y(77))/(y(76)*y(76));
T(347) = getPowerDeriv(y(77)/y(76),T(3),1);
T(348) = getPowerDeriv(y(76),T(6),1);
T(349) = (-y(78))/(y(77)*y(77));
T(350) = getPowerDeriv(y(78)/y(77),T(3),1);
T(351) = getPowerDeriv(y(77),T(6),1);
T(352) = (-y(79))/(y(78)*y(78));
T(353) = getPowerDeriv(y(79)/y(78),T(3),1);
T(354) = getPowerDeriv(y(78),T(6),1);
T(355) = (-y(80))/(y(79)*y(79));
T(356) = getPowerDeriv(y(80)/y(79),T(3),1);
T(357) = getPowerDeriv(y(79),T(6),1);
T(358) = (-y(81))/(y(80)*y(80));
T(359) = getPowerDeriv(y(81)/y(80),T(3),1);
T(360) = getPowerDeriv(y(80),T(6),1);
T(361) = (-y(82))/(y(81)*y(81));
T(362) = getPowerDeriv(y(82)/y(81),T(3),1);
T(363) = getPowerDeriv(y(81),T(6),1);
T(364) = (-y(83))/(y(82)*y(82));
T(365) = getPowerDeriv(y(83)/y(82),T(3),1);
T(366) = getPowerDeriv(y(82),T(6),1);
T(367) = (-y(84))/(y(83)*y(83));
T(368) = getPowerDeriv(y(84)/y(83),T(3),1);
T(369) = getPowerDeriv(y(83),T(6),1);
T(370) = (-y(85))/(y(84)*y(84));
T(371) = getPowerDeriv(y(85)/y(84),T(3),1);
T(372) = getPowerDeriv(y(84),T(6),1);
T(373) = (-y(86))/(y(85)*y(85));
T(374) = getPowerDeriv(y(86)/y(85),T(3),1);
T(375) = getPowerDeriv(y(85),T(6),1);
T(376) = (-y(87))/(y(86)*y(86));
T(377) = getPowerDeriv(y(87)/y(86),T(3),1);
T(378) = getPowerDeriv(y(86),T(6),1);
T(379) = (-y(88))/(y(87)*y(87));
T(380) = getPowerDeriv(y(88)/y(87),T(3),1);
T(381) = getPowerDeriv(y(87),T(6),1);
T(382) = (-y(89))/(y(88)*y(88));
T(383) = getPowerDeriv(y(89)/y(88),T(3),1);
T(384) = getPowerDeriv(y(88),T(6),1);
T(385) = (-y(90))/(y(89)*y(89));
T(386) = getPowerDeriv(y(90)/y(89),T(3),1);
T(387) = getPowerDeriv(y(89),T(6),1);
T(388) = (-y(91))/(y(90)*y(90));
T(389) = getPowerDeriv(y(91)/y(90),T(3),1);
T(390) = getPowerDeriv(y(90),T(6),1);
T(391) = (-y(92))/(y(91)*y(91));
T(392) = getPowerDeriv(y(92)/y(91),T(3),1);
T(393) = getPowerDeriv(y(91),T(6),1);
T(394) = (-y(93))/(y(92)*y(92));
T(395) = getPowerDeriv(y(93)/y(92),T(3),1);
T(396) = getPowerDeriv(y(92),T(6),1);
T(397) = (-y(94))/(y(93)*y(93));
T(398) = getPowerDeriv(y(94)/y(93),T(3),1);
T(399) = getPowerDeriv(y(93),T(6),1);
T(400) = (-y(95))/(y(94)*y(94));
T(401) = getPowerDeriv(y(95)/y(94),T(3),1);
T(402) = getPowerDeriv(y(94),T(6),1);
T(403) = (-y(96))/(y(95)*y(95));
T(404) = getPowerDeriv(y(96)/y(95),T(3),1);
T(405) = getPowerDeriv(y(95),T(6),1);
T(406) = (-y(97))/(y(96)*y(96));
T(407) = getPowerDeriv(y(97)/y(96),T(3),1);
T(408) = getPowerDeriv(y(96),T(6),1);
T(409) = (-y(98))/(y(97)*y(97));
T(410) = getPowerDeriv(y(98)/y(97),T(3),1);
T(411) = getPowerDeriv(y(97),T(6),1);
T(412) = (-y(99))/(y(98)*y(98));
T(413) = getPowerDeriv(y(99)/y(98),T(3),1);
T(414) = getPowerDeriv(y(98),T(6),1);
T(415) = (-y(100))/(y(99)*y(99));
T(416) = getPowerDeriv(y(100)/y(99),T(3),1);
T(417) = getPowerDeriv(y(99),T(6),1);
T(418) = (-y(101))/(y(100)*y(100));
T(419) = getPowerDeriv(y(101)/y(100),T(3),1);
T(420) = getPowerDeriv(y(100),T(6),1);
T(421) = (-y(102))/(y(101)*y(101));
T(422) = getPowerDeriv(y(102)/y(101),T(3),1);
T(423) = getPowerDeriv(y(101),T(6),1);
T(424) = (-y(103))/(y(102)*y(102));
T(425) = getPowerDeriv(y(103)/y(102),T(3),1);
T(426) = getPowerDeriv(y(102),T(6),1);
T(427) = (-y(104))/(y(103)*y(103));
T(428) = getPowerDeriv(y(104)/y(103),T(3),1);
T(429) = getPowerDeriv(y(103),T(6),1);
T(430) = (-y(105))/(y(104)*y(104));
T(431) = getPowerDeriv(y(105)/y(104),T(3),1);
T(432) = getPowerDeriv(y(104),T(6),1);
T(433) = (-y(106))/(y(105)*y(105));
T(434) = getPowerDeriv(y(106)/y(105),T(3),1);
T(435) = getPowerDeriv(y(105),T(6),1);
T(436) = (-y(107))/(y(106)*y(106));
T(437) = getPowerDeriv(y(107)/y(106),T(3),1);
T(438) = getPowerDeriv(y(106),T(6),1);
T(439) = (-y(108))/(y(107)*y(107));
T(440) = getPowerDeriv(y(108)/y(107),T(3),1);
T(441) = getPowerDeriv(y(107),T(6),1);
T(442) = (-y(109))/(y(108)*y(108));
T(443) = getPowerDeriv(y(109)/y(108),T(3),1);
T(444) = getPowerDeriv(y(108),T(6),1);
T(445) = (-y(110))/(y(109)*y(109));
T(446) = getPowerDeriv(y(110)/y(109),T(3),1);
T(447) = getPowerDeriv(y(109),T(6),1);
T(448) = (-y(111))/(y(110)*y(110));
T(449) = getPowerDeriv(y(111)/y(110),T(3),1);
T(450) = getPowerDeriv(y(110),T(6),1);
T(451) = (-y(112))/(y(111)*y(111));
T(452) = getPowerDeriv(y(112)/y(111),T(3),1);
T(453) = getPowerDeriv(y(111),T(6),1);
T(454) = y(2)/params(105)/T(172);
T(455) = y(3)/params(106)/T(173);
T(456) = y(4)/params(107)/T(174);
T(457) = y(5)/params(108)/T(175);
T(458) = y(6)/params(109)/T(176);
T(459) = y(7)/params(110)/T(177);
T(460) = y(8)/params(111)/T(178);
T(461) = y(9)/params(112)/T(179);
T(462) = y(10)/params(113)/T(180);
T(463) = y(11)/params(114)/T(181);
T(464) = y(12)/params(115)/T(182);
T(465) = y(13)/params(116)/T(183);
T(466) = y(14)/params(117)/T(184);
T(467) = y(15)/params(118)/T(185);
T(468) = y(16)/params(119)/T(186);
T(469) = y(17)/params(120)/T(187);
T(470) = y(18)/params(121)/T(188);
T(471) = y(19)/params(122)/T(189);
T(472) = y(20)/params(123)/T(190);
T(473) = y(21)/params(124)/T(191);
T(474) = y(22)/params(125)/T(192);
T(475) = y(23)/params(126)/T(193);
T(476) = y(24)/params(127)/T(194);
T(477) = y(25)/params(128)/T(195);
T(478) = y(26)/params(129)/T(196);
T(479) = y(27)/params(130)/T(197);
T(480) = y(28)/params(131)/T(198);
T(481) = y(29)/params(132)/T(199);
T(482) = y(30)/params(133)/T(200);
T(483) = y(31)/params(134)/T(201);
T(484) = y(32)/params(135)/T(202);
T(485) = y(33)/params(136)/T(204);
T(486) = y(34)/params(137)/T(205);
T(487) = y(35)/params(138)/T(206);
T(488) = y(36)/params(139)/T(207);
T(489) = y(37)/params(140)/T(208);
T(490) = y(38)/params(141)/T(209);
T(491) = y(39)/params(142)/T(210);
T(492) = y(40)/params(143)/T(211);
T(493) = y(41)/params(144)/T(216);
T(494) = y(42)/params(145)/T(217);
T(495) = y(43)/params(146)/T(218);
T(496) = y(44)/params(147)/T(219);
T(497) = y(45)/params(148)/T(220);
T(498) = y(46)/params(149)/T(221);
T(499) = y(47)/params(150)/T(222);
T(500) = y(48)/params(151)/T(223);
T(501) = y(49)/params(152)/T(224);
T(502) = y(50)/params(153)/T(225);
T(503) = y(51)/params(154)/T(226);
T(504) = y(52)/params(155)/T(227);
T(505) = y(53)/params(156)/T(228);
T(506) = y(54)/params(157)/T(229);
T(507) = y(55)/params(158)/T(230);
T(508) = y(56)/params(159)/T(231);
T(509) = (y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(273)+y(276)+y(270)*y(278)+y(278)*params(3));
T(510) = y(276)*y(276);
T(511) = params(223)*getPowerDeriv(y(276)*params(223),(-params(1)),1);
T(512) = (-(params(223)*params(224)*y(278)))/(y(276)*params(223)*y(276)*params(223));
T(513) = getPowerDeriv(params(224)*y(278)/(y(276)*params(223)),params(1)-1,1);
T(514) = params(223)*getPowerDeriv(y(276)*params(223),1-params(1),1);
T(515) = params(224)*getPowerDeriv(params(224)*y(278),params(1),1);
T(516) = T(232)*T(232);

end