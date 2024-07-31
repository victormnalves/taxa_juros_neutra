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

assert(length(T) >= 412);

T = dynare_ss_1970.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(238) = (-y(58))/(y(57)*y(57));
T(239) = getPowerDeriv(y(58)/y(57),T(3),1);
T(240) = getPowerDeriv(y(57),T(6),1);
T(241) = (-y(59))/(y(58)*y(58));
T(242) = getPowerDeriv(y(59)/y(58),T(3),1);
T(243) = getPowerDeriv(y(58),T(6),1);
T(244) = (-y(60))/(y(59)*y(59));
T(245) = getPowerDeriv(y(60)/y(59),T(3),1);
T(246) = getPowerDeriv(y(59),T(6),1);
T(247) = (-y(61))/(y(60)*y(60));
T(248) = getPowerDeriv(y(61)/y(60),T(3),1);
T(249) = getPowerDeriv(y(60),T(6),1);
T(250) = (-y(62))/(y(61)*y(61));
T(251) = getPowerDeriv(y(62)/y(61),T(3),1);
T(252) = getPowerDeriv(y(61),T(6),1);
T(253) = (-y(63))/(y(62)*y(62));
T(254) = getPowerDeriv(y(63)/y(62),T(3),1);
T(255) = getPowerDeriv(y(62),T(6),1);
T(256) = (-y(64))/(y(63)*y(63));
T(257) = getPowerDeriv(y(64)/y(63),T(3),1);
T(258) = getPowerDeriv(y(63),T(6),1);
T(259) = (-y(65))/(y(64)*y(64));
T(260) = getPowerDeriv(y(65)/y(64),T(3),1);
T(261) = getPowerDeriv(y(64),T(6),1);
T(262) = (-y(66))/(y(65)*y(65));
T(263) = getPowerDeriv(y(66)/y(65),T(3),1);
T(264) = getPowerDeriv(y(65),T(6),1);
T(265) = (-y(67))/(y(66)*y(66));
T(266) = getPowerDeriv(y(67)/y(66),T(3),1);
T(267) = getPowerDeriv(y(66),T(6),1);
T(268) = (-y(68))/(y(67)*y(67));
T(269) = getPowerDeriv(y(68)/y(67),T(3),1);
T(270) = getPowerDeriv(y(67),T(6),1);
T(271) = (-y(69))/(y(68)*y(68));
T(272) = getPowerDeriv(y(69)/y(68),T(3),1);
T(273) = getPowerDeriv(y(68),T(6),1);
T(274) = (-y(70))/(y(69)*y(69));
T(275) = getPowerDeriv(y(70)/y(69),T(3),1);
T(276) = getPowerDeriv(y(69),T(6),1);
T(277) = (-y(71))/(y(70)*y(70));
T(278) = getPowerDeriv(y(71)/y(70),T(3),1);
T(279) = getPowerDeriv(y(70),T(6),1);
T(280) = (-y(72))/(y(71)*y(71));
T(281) = getPowerDeriv(y(72)/y(71),T(3),1);
T(282) = getPowerDeriv(y(71),T(6),1);
T(283) = (-y(73))/(y(72)*y(72));
T(284) = getPowerDeriv(y(73)/y(72),T(3),1);
T(285) = getPowerDeriv(y(72),T(6),1);
T(286) = (-y(74))/(y(73)*y(73));
T(287) = getPowerDeriv(y(74)/y(73),T(3),1);
T(288) = getPowerDeriv(y(73),T(6),1);
T(289) = (-y(75))/(y(74)*y(74));
T(290) = getPowerDeriv(y(75)/y(74),T(3),1);
T(291) = getPowerDeriv(y(74),T(6),1);
T(292) = (-y(76))/(y(75)*y(75));
T(293) = getPowerDeriv(y(76)/y(75),T(3),1);
T(294) = getPowerDeriv(y(75),T(6),1);
T(295) = (-y(77))/(y(76)*y(76));
T(296) = getPowerDeriv(y(77)/y(76),T(3),1);
T(297) = getPowerDeriv(y(76),T(6),1);
T(298) = (-y(78))/(y(77)*y(77));
T(299) = getPowerDeriv(y(78)/y(77),T(3),1);
T(300) = getPowerDeriv(y(77),T(6),1);
T(301) = (-y(79))/(y(78)*y(78));
T(302) = getPowerDeriv(y(79)/y(78),T(3),1);
T(303) = getPowerDeriv(y(78),T(6),1);
T(304) = (-y(80))/(y(79)*y(79));
T(305) = getPowerDeriv(y(80)/y(79),T(3),1);
T(306) = getPowerDeriv(y(79),T(6),1);
T(307) = (-y(81))/(y(80)*y(80));
T(308) = getPowerDeriv(y(81)/y(80),T(3),1);
T(309) = getPowerDeriv(y(80),T(6),1);
T(310) = (-y(82))/(y(81)*y(81));
T(311) = getPowerDeriv(y(82)/y(81),T(3),1);
T(312) = getPowerDeriv(y(81),T(6),1);
T(313) = (-y(83))/(y(82)*y(82));
T(314) = getPowerDeriv(y(83)/y(82),T(3),1);
T(315) = getPowerDeriv(y(82),T(6),1);
T(316) = (-y(84))/(y(83)*y(83));
T(317) = getPowerDeriv(y(84)/y(83),T(3),1);
T(318) = getPowerDeriv(y(83),T(6),1);
T(319) = (-y(85))/(y(84)*y(84));
T(320) = getPowerDeriv(y(85)/y(84),T(3),1);
T(321) = getPowerDeriv(y(84),T(6),1);
T(322) = (-y(86))/(y(85)*y(85));
T(323) = getPowerDeriv(y(86)/y(85),T(3),1);
T(324) = getPowerDeriv(y(85),T(6),1);
T(325) = (-y(87))/(y(86)*y(86));
T(326) = getPowerDeriv(y(87)/y(86),T(3),1);
T(327) = getPowerDeriv(y(86),T(6),1);
T(328) = (-y(88))/(y(87)*y(87));
T(329) = getPowerDeriv(y(88)/y(87),T(3),1);
T(330) = getPowerDeriv(y(87),T(6),1);
T(331) = (-y(89))/(y(88)*y(88));
T(332) = getPowerDeriv(y(89)/y(88),T(3),1);
T(333) = getPowerDeriv(y(88),T(6),1);
T(334) = (-y(90))/(y(89)*y(89));
T(335) = getPowerDeriv(y(90)/y(89),T(3),1);
T(336) = getPowerDeriv(y(89),T(6),1);
T(337) = (-y(91))/(y(90)*y(90));
T(338) = getPowerDeriv(y(91)/y(90),T(3),1);
T(339) = getPowerDeriv(y(90),T(6),1);
T(340) = (-y(92))/(y(91)*y(91));
T(341) = getPowerDeriv(y(92)/y(91),T(3),1);
T(342) = getPowerDeriv(y(91),T(6),1);
T(343) = (-y(93))/(y(92)*y(92));
T(344) = getPowerDeriv(y(93)/y(92),T(3),1);
T(345) = getPowerDeriv(y(92),T(6),1);
T(346) = (-y(94))/(y(93)*y(93));
T(347) = getPowerDeriv(y(94)/y(93),T(3),1);
T(348) = getPowerDeriv(y(93),T(6),1);
T(349) = (-y(95))/(y(94)*y(94));
T(350) = getPowerDeriv(y(95)/y(94),T(3),1);
T(351) = getPowerDeriv(y(94),T(6),1);
T(352) = (-y(96))/(y(95)*y(95));
T(353) = getPowerDeriv(y(96)/y(95),T(3),1);
T(354) = getPowerDeriv(y(95),T(6),1);
T(355) = (-y(97))/(y(96)*y(96));
T(356) = getPowerDeriv(y(97)/y(96),T(3),1);
T(357) = getPowerDeriv(y(96),T(6),1);
T(358) = (-y(98))/(y(97)*y(97));
T(359) = getPowerDeriv(y(98)/y(97),T(3),1);
T(360) = getPowerDeriv(y(97),T(6),1);
T(361) = (-y(99))/(y(98)*y(98));
T(362) = getPowerDeriv(y(99)/y(98),T(3),1);
T(363) = getPowerDeriv(y(98),T(6),1);
T(364) = (-y(100))/(y(99)*y(99));
T(365) = getPowerDeriv(y(100)/y(99),T(3),1);
T(366) = getPowerDeriv(y(99),T(6),1);
T(367) = (-y(101))/(y(100)*y(100));
T(368) = getPowerDeriv(y(101)/y(100),T(3),1);
T(369) = getPowerDeriv(y(100),T(6),1);
T(370) = (-y(102))/(y(101)*y(101));
T(371) = getPowerDeriv(y(102)/y(101),T(3),1);
T(372) = getPowerDeriv(y(101),T(6),1);
T(373) = (-y(103))/(y(102)*y(102));
T(374) = getPowerDeriv(y(103)/y(102),T(3),1);
T(375) = getPowerDeriv(y(102),T(6),1);
T(376) = (-y(104))/(y(103)*y(103));
T(377) = getPowerDeriv(y(104)/y(103),T(3),1);
T(378) = getPowerDeriv(y(103),T(6),1);
T(379) = (-y(105))/(y(104)*y(104));
T(380) = getPowerDeriv(y(105)/y(104),T(3),1);
T(381) = getPowerDeriv(y(104),T(6),1);
T(382) = (-y(106))/(y(105)*y(105));
T(383) = getPowerDeriv(y(106)/y(105),T(3),1);
T(384) = getPowerDeriv(y(105),T(6),1);
T(385) = (-y(107))/(y(106)*y(106));
T(386) = getPowerDeriv(y(107)/y(106),T(3),1);
T(387) = getPowerDeriv(y(106),T(6),1);
T(388) = (-y(108))/(y(107)*y(107));
T(389) = getPowerDeriv(y(108)/y(107),T(3),1);
T(390) = getPowerDeriv(y(107),T(6),1);
T(391) = (-y(109))/(y(108)*y(108));
T(392) = getPowerDeriv(y(109)/y(108),T(3),1);
T(393) = getPowerDeriv(y(108),T(6),1);
T(394) = (-y(110))/(y(109)*y(109));
T(395) = getPowerDeriv(y(110)/y(109),T(3),1);
T(396) = getPowerDeriv(y(109),T(6),1);
T(397) = (-y(111))/(y(110)*y(110));
T(398) = getPowerDeriv(y(111)/y(110),T(3),1);
T(399) = getPowerDeriv(y(110),T(6),1);
T(400) = (-y(112))/(y(111)*y(111));
T(401) = getPowerDeriv(y(112)/y(111),T(3),1);
T(402) = getPowerDeriv(y(111),T(6),1);
T(403) = (y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(273)+y(276)+y(270)*y(278)+y(278)*params(3));
T(404) = y(276)*y(276);
T(405) = (1-params(1))*params(223)*getPowerDeriv(y(276)*params(223),T(212),1);
T(406) = getPowerDeriv(T(213),1/(params(5)-1),1);
T(407) = y(271)*T(405)*T(406);
T(408) = getPowerDeriv(y(276),T(217),1);
T(409) = getPowerDeriv(T(213),params(5)/(params(5)-1),1);
T(410) = params(1)*params(224)*getPowerDeriv(params(224)*y(278),T(212),1);
T(411) = y(271)*T(406)*T(410);
T(412) = getPowerDeriv(y(278),T(217),1);

end
