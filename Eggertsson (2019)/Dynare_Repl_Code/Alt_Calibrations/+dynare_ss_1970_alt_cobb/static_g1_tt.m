function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 292);

T = dynare_ss_1970_alt_cobb.static_resid_tt(T, y, x, params);

T(234) = getPowerDeriv(y(58)/y(57),T(3),1);
T(235) = getPowerDeriv(y(59)/y(58),T(3),1);
T(236) = getPowerDeriv(y(60)/y(59),T(3),1);
T(237) = getPowerDeriv(y(61)/y(60),T(3),1);
T(238) = getPowerDeriv(y(62)/y(61),T(3),1);
T(239) = getPowerDeriv(y(63)/y(62),T(3),1);
T(240) = getPowerDeriv(y(64)/y(63),T(3),1);
T(241) = getPowerDeriv(y(65)/y(64),T(3),1);
T(242) = getPowerDeriv(y(66)/y(65),T(3),1);
T(243) = getPowerDeriv(y(67)/y(66),T(3),1);
T(244) = getPowerDeriv(y(68)/y(67),T(3),1);
T(245) = getPowerDeriv(y(69)/y(68),T(3),1);
T(246) = getPowerDeriv(y(70)/y(69),T(3),1);
T(247) = getPowerDeriv(y(71)/y(70),T(3),1);
T(248) = getPowerDeriv(y(72)/y(71),T(3),1);
T(249) = getPowerDeriv(y(73)/y(72),T(3),1);
T(250) = getPowerDeriv(y(74)/y(73),T(3),1);
T(251) = getPowerDeriv(y(75)/y(74),T(3),1);
T(252) = getPowerDeriv(y(76)/y(75),T(3),1);
T(253) = getPowerDeriv(y(77)/y(76),T(3),1);
T(254) = getPowerDeriv(y(78)/y(77),T(3),1);
T(255) = getPowerDeriv(y(79)/y(78),T(3),1);
T(256) = getPowerDeriv(y(80)/y(79),T(3),1);
T(257) = getPowerDeriv(y(81)/y(80),T(3),1);
T(258) = getPowerDeriv(y(82)/y(81),T(3),1);
T(259) = getPowerDeriv(y(83)/y(82),T(3),1);
T(260) = getPowerDeriv(y(84)/y(83),T(3),1);
T(261) = getPowerDeriv(y(85)/y(84),T(3),1);
T(262) = getPowerDeriv(y(86)/y(85),T(3),1);
T(263) = getPowerDeriv(y(87)/y(86),T(3),1);
T(264) = getPowerDeriv(y(88)/y(87),T(3),1);
T(265) = getPowerDeriv(y(89)/y(88),T(3),1);
T(266) = getPowerDeriv(y(90)/y(89),T(3),1);
T(267) = getPowerDeriv(y(91)/y(90),T(3),1);
T(268) = getPowerDeriv(y(92)/y(91),T(3),1);
T(269) = getPowerDeriv(y(93)/y(92),T(3),1);
T(270) = getPowerDeriv(y(94)/y(93),T(3),1);
T(271) = getPowerDeriv(y(95)/y(94),T(3),1);
T(272) = getPowerDeriv(y(96)/y(95),T(3),1);
T(273) = getPowerDeriv(y(97)/y(96),T(3),1);
T(274) = getPowerDeriv(y(98)/y(97),T(3),1);
T(275) = getPowerDeriv(y(99)/y(98),T(3),1);
T(276) = getPowerDeriv(y(100)/y(99),T(3),1);
T(277) = getPowerDeriv(y(101)/y(100),T(3),1);
T(278) = getPowerDeriv(y(102)/y(101),T(3),1);
T(279) = getPowerDeriv(y(103)/y(102),T(3),1);
T(280) = getPowerDeriv(y(104)/y(103),T(3),1);
T(281) = getPowerDeriv(y(105)/y(104),T(3),1);
T(282) = getPowerDeriv(y(106)/y(105),T(3),1);
T(283) = getPowerDeriv(y(107)/y(106),T(3),1);
T(284) = getPowerDeriv(y(108)/y(107),T(3),1);
T(285) = getPowerDeriv(y(109)/y(108),T(3),1);
T(286) = getPowerDeriv(y(110)/y(109),T(3),1);
T(287) = getPowerDeriv(y(111)/y(110),T(3),1);
T(288) = getPowerDeriv(y(112)/y(111),T(3),1);
T(289) = (y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(273)+y(276)+y(270)*y(278)+y(278)*params(3));
T(290) = y(276)*y(276);
T(291) = getPowerDeriv(params(224)*y(278)/(y(276)*params(223)),params(1)-1,1);
T(292) = params(224)*getPowerDeriv(params(224)*y(278),params(1),1);

end
