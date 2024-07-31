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

assert(length(T) >= 413);

T = dynare_ss_2015.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(239) = (-y(58))/(y(57)*y(57));
T(240) = getPowerDeriv(y(58)/y(57),T(3),1);
T(241) = getPowerDeriv(y(57),T(6),1);
T(242) = (-y(59))/(y(58)*y(58));
T(243) = getPowerDeriv(y(59)/y(58),T(3),1);
T(244) = getPowerDeriv(y(58),T(6),1);
T(245) = (-y(60))/(y(59)*y(59));
T(246) = getPowerDeriv(y(60)/y(59),T(3),1);
T(247) = getPowerDeriv(y(59),T(6),1);
T(248) = (-y(61))/(y(60)*y(60));
T(249) = getPowerDeriv(y(61)/y(60),T(3),1);
T(250) = getPowerDeriv(y(60),T(6),1);
T(251) = (-y(62))/(y(61)*y(61));
T(252) = getPowerDeriv(y(62)/y(61),T(3),1);
T(253) = getPowerDeriv(y(61),T(6),1);
T(254) = (-y(63))/(y(62)*y(62));
T(255) = getPowerDeriv(y(63)/y(62),T(3),1);
T(256) = getPowerDeriv(y(62),T(6),1);
T(257) = (-y(64))/(y(63)*y(63));
T(258) = getPowerDeriv(y(64)/y(63),T(3),1);
T(259) = getPowerDeriv(y(63),T(6),1);
T(260) = (-y(65))/(y(64)*y(64));
T(261) = getPowerDeriv(y(65)/y(64),T(3),1);
T(262) = getPowerDeriv(y(64),T(6),1);
T(263) = (-y(66))/(y(65)*y(65));
T(264) = getPowerDeriv(y(66)/y(65),T(3),1);
T(265) = getPowerDeriv(y(65),T(6),1);
T(266) = (-y(67))/(y(66)*y(66));
T(267) = getPowerDeriv(y(67)/y(66),T(3),1);
T(268) = getPowerDeriv(y(66),T(6),1);
T(269) = (-y(68))/(y(67)*y(67));
T(270) = getPowerDeriv(y(68)/y(67),T(3),1);
T(271) = getPowerDeriv(y(67),T(6),1);
T(272) = (-y(69))/(y(68)*y(68));
T(273) = getPowerDeriv(y(69)/y(68),T(3),1);
T(274) = getPowerDeriv(y(68),T(6),1);
T(275) = (-y(70))/(y(69)*y(69));
T(276) = getPowerDeriv(y(70)/y(69),T(3),1);
T(277) = getPowerDeriv(y(69),T(6),1);
T(278) = (-y(71))/(y(70)*y(70));
T(279) = getPowerDeriv(y(71)/y(70),T(3),1);
T(280) = getPowerDeriv(y(70),T(6),1);
T(281) = (-y(72))/(y(71)*y(71));
T(282) = getPowerDeriv(y(72)/y(71),T(3),1);
T(283) = getPowerDeriv(y(71),T(6),1);
T(284) = (-y(73))/(y(72)*y(72));
T(285) = getPowerDeriv(y(73)/y(72),T(3),1);
T(286) = getPowerDeriv(y(72),T(6),1);
T(287) = (-y(74))/(y(73)*y(73));
T(288) = getPowerDeriv(y(74)/y(73),T(3),1);
T(289) = getPowerDeriv(y(73),T(6),1);
T(290) = (-y(75))/(y(74)*y(74));
T(291) = getPowerDeriv(y(75)/y(74),T(3),1);
T(292) = getPowerDeriv(y(74),T(6),1);
T(293) = (-y(76))/(y(75)*y(75));
T(294) = getPowerDeriv(y(76)/y(75),T(3),1);
T(295) = getPowerDeriv(y(75),T(6),1);
T(296) = (-y(77))/(y(76)*y(76));
T(297) = getPowerDeriv(y(77)/y(76),T(3),1);
T(298) = getPowerDeriv(y(76),T(6),1);
T(299) = (-y(78))/(y(77)*y(77));
T(300) = getPowerDeriv(y(78)/y(77),T(3),1);
T(301) = getPowerDeriv(y(77),T(6),1);
T(302) = (-y(79))/(y(78)*y(78));
T(303) = getPowerDeriv(y(79)/y(78),T(3),1);
T(304) = getPowerDeriv(y(78),T(6),1);
T(305) = (-y(80))/(y(79)*y(79));
T(306) = getPowerDeriv(y(80)/y(79),T(3),1);
T(307) = getPowerDeriv(y(79),T(6),1);
T(308) = (-y(81))/(y(80)*y(80));
T(309) = getPowerDeriv(y(81)/y(80),T(3),1);
T(310) = getPowerDeriv(y(80),T(6),1);
T(311) = (-y(82))/(y(81)*y(81));
T(312) = getPowerDeriv(y(82)/y(81),T(3),1);
T(313) = getPowerDeriv(y(81),T(6),1);
T(314) = (-y(83))/(y(82)*y(82));
T(315) = getPowerDeriv(y(83)/y(82),T(3),1);
T(316) = getPowerDeriv(y(82),T(6),1);
T(317) = (-y(84))/(y(83)*y(83));
T(318) = getPowerDeriv(y(84)/y(83),T(3),1);
T(319) = getPowerDeriv(y(83),T(6),1);
T(320) = (-y(85))/(y(84)*y(84));
T(321) = getPowerDeriv(y(85)/y(84),T(3),1);
T(322) = getPowerDeriv(y(84),T(6),1);
T(323) = (-y(86))/(y(85)*y(85));
T(324) = getPowerDeriv(y(86)/y(85),T(3),1);
T(325) = getPowerDeriv(y(85),T(6),1);
T(326) = (-y(87))/(y(86)*y(86));
T(327) = getPowerDeriv(y(87)/y(86),T(3),1);
T(328) = getPowerDeriv(y(86),T(6),1);
T(329) = (-y(88))/(y(87)*y(87));
T(330) = getPowerDeriv(y(88)/y(87),T(3),1);
T(331) = getPowerDeriv(y(87),T(6),1);
T(332) = (-y(89))/(y(88)*y(88));
T(333) = getPowerDeriv(y(89)/y(88),T(3),1);
T(334) = getPowerDeriv(y(88),T(6),1);
T(335) = (-y(90))/(y(89)*y(89));
T(336) = getPowerDeriv(y(90)/y(89),T(3),1);
T(337) = getPowerDeriv(y(89),T(6),1);
T(338) = (-y(91))/(y(90)*y(90));
T(339) = getPowerDeriv(y(91)/y(90),T(3),1);
T(340) = getPowerDeriv(y(90),T(6),1);
T(341) = (-y(92))/(y(91)*y(91));
T(342) = getPowerDeriv(y(92)/y(91),T(3),1);
T(343) = getPowerDeriv(y(91),T(6),1);
T(344) = (-y(93))/(y(92)*y(92));
T(345) = getPowerDeriv(y(93)/y(92),T(3),1);
T(346) = getPowerDeriv(y(92),T(6),1);
T(347) = (-y(94))/(y(93)*y(93));
T(348) = getPowerDeriv(y(94)/y(93),T(3),1);
T(349) = getPowerDeriv(y(93),T(6),1);
T(350) = (-y(95))/(y(94)*y(94));
T(351) = getPowerDeriv(y(95)/y(94),T(3),1);
T(352) = getPowerDeriv(y(94),T(6),1);
T(353) = (-y(96))/(y(95)*y(95));
T(354) = getPowerDeriv(y(96)/y(95),T(3),1);
T(355) = getPowerDeriv(y(95),T(6),1);
T(356) = (-y(97))/(y(96)*y(96));
T(357) = getPowerDeriv(y(97)/y(96),T(3),1);
T(358) = getPowerDeriv(y(96),T(6),1);
T(359) = (-y(98))/(y(97)*y(97));
T(360) = getPowerDeriv(y(98)/y(97),T(3),1);
T(361) = getPowerDeriv(y(97),T(6),1);
T(362) = (-y(99))/(y(98)*y(98));
T(363) = getPowerDeriv(y(99)/y(98),T(3),1);
T(364) = getPowerDeriv(y(98),T(6),1);
T(365) = (-y(100))/(y(99)*y(99));
T(366) = getPowerDeriv(y(100)/y(99),T(3),1);
T(367) = getPowerDeriv(y(99),T(6),1);
T(368) = (-y(101))/(y(100)*y(100));
T(369) = getPowerDeriv(y(101)/y(100),T(3),1);
T(370) = getPowerDeriv(y(100),T(6),1);
T(371) = (-y(102))/(y(101)*y(101));
T(372) = getPowerDeriv(y(102)/y(101),T(3),1);
T(373) = getPowerDeriv(y(101),T(6),1);
T(374) = (-y(103))/(y(102)*y(102));
T(375) = getPowerDeriv(y(103)/y(102),T(3),1);
T(376) = getPowerDeriv(y(102),T(6),1);
T(377) = (-y(104))/(y(103)*y(103));
T(378) = getPowerDeriv(y(104)/y(103),T(3),1);
T(379) = getPowerDeriv(y(103),T(6),1);
T(380) = (-y(105))/(y(104)*y(104));
T(381) = getPowerDeriv(y(105)/y(104),T(3),1);
T(382) = getPowerDeriv(y(104),T(6),1);
T(383) = (-y(106))/(y(105)*y(105));
T(384) = getPowerDeriv(y(106)/y(105),T(3),1);
T(385) = getPowerDeriv(y(105),T(6),1);
T(386) = (-y(107))/(y(106)*y(106));
T(387) = getPowerDeriv(y(107)/y(106),T(3),1);
T(388) = getPowerDeriv(y(106),T(6),1);
T(389) = (-y(108))/(y(107)*y(107));
T(390) = getPowerDeriv(y(108)/y(107),T(3),1);
T(391) = getPowerDeriv(y(107),T(6),1);
T(392) = (-y(109))/(y(108)*y(108));
T(393) = getPowerDeriv(y(109)/y(108),T(3),1);
T(394) = getPowerDeriv(y(108),T(6),1);
T(395) = (-y(110))/(y(109)*y(109));
T(396) = getPowerDeriv(y(110)/y(109),T(3),1);
T(397) = getPowerDeriv(y(109),T(6),1);
T(398) = (-y(111))/(y(110)*y(110));
T(399) = getPowerDeriv(y(111)/y(110),T(3),1);
T(400) = getPowerDeriv(y(110),T(6),1);
T(401) = (-y(112))/(y(111)*y(111));
T(402) = getPowerDeriv(y(112)/y(111),T(3),1);
T(403) = getPowerDeriv(y(111),T(6),1);
T(404) = (y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(273)+y(276)+y(270)*y(278)+y(278)*params(3));
T(405) = y(276)*y(276);
T(406) = (1-params(1))*params(223)*getPowerDeriv(y(276)*params(223),T(214),1);
T(407) = getPowerDeriv(T(215),1/(params(5)-1),1);
T(408) = y(271)*T(406)*T(407);
T(409) = getPowerDeriv(y(276),T(219),1);
T(410) = getPowerDeriv(T(215),params(5)/(params(5)-1),1);
T(411) = params(1)*params(224)*getPowerDeriv(params(224)*y(278),T(214),1);
T(412) = y(271)*T(407)*T(411);
T(413) = getPowerDeriv(y(278),T(219),1);

end
