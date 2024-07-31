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

assert(length(T) >= 243);

T = dynare_transition_alt.static_resid_tt(T, y, x, params);

T(183) = getPowerDeriv(y(58)/y(57),T(2),1);
T(184) = getPowerDeriv(y(59)/y(58),T(2),1);
T(185) = getPowerDeriv(y(60)/y(59),T(2),1);
T(186) = getPowerDeriv(y(61)/y(60),T(2),1);
T(187) = getPowerDeriv(y(62)/y(61),T(2),1);
T(188) = getPowerDeriv(y(63)/y(62),T(2),1);
T(189) = getPowerDeriv(y(64)/y(63),T(2),1);
T(190) = getPowerDeriv(y(65)/y(64),T(2),1);
T(191) = getPowerDeriv(y(66)/y(65),T(2),1);
T(192) = getPowerDeriv(y(67)/y(66),T(2),1);
T(193) = getPowerDeriv(y(68)/y(67),T(2),1);
T(194) = getPowerDeriv(y(69)/y(68),T(2),1);
T(195) = getPowerDeriv(y(70)/y(69),T(2),1);
T(196) = getPowerDeriv(y(71)/y(70),T(2),1);
T(197) = getPowerDeriv(y(72)/y(71),T(2),1);
T(198) = getPowerDeriv(y(73)/y(72),T(2),1);
T(199) = getPowerDeriv(y(74)/y(73),T(2),1);
T(200) = getPowerDeriv(y(75)/y(74),T(2),1);
T(201) = getPowerDeriv(y(76)/y(75),T(2),1);
T(202) = getPowerDeriv(y(77)/y(76),T(2),1);
T(203) = getPowerDeriv(y(78)/y(77),T(2),1);
T(204) = getPowerDeriv(y(79)/y(78),T(2),1);
T(205) = getPowerDeriv(y(80)/y(79),T(2),1);
T(206) = getPowerDeriv(y(81)/y(80),T(2),1);
T(207) = getPowerDeriv(y(82)/y(81),T(2),1);
T(208) = getPowerDeriv(y(83)/y(82),T(2),1);
T(209) = getPowerDeriv(y(84)/y(83),T(2),1);
T(210) = getPowerDeriv(y(85)/y(84),T(2),1);
T(211) = getPowerDeriv(y(86)/y(85),T(2),1);
T(212) = getPowerDeriv(y(87)/y(86),T(2),1);
T(213) = getPowerDeriv(y(88)/y(87),T(2),1);
T(214) = getPowerDeriv(y(89)/y(88),T(2),1);
T(215) = getPowerDeriv(y(90)/y(89),T(2),1);
T(216) = getPowerDeriv(y(91)/y(90),T(2),1);
T(217) = getPowerDeriv(y(92)/y(91),T(2),1);
T(218) = getPowerDeriv(y(93)/y(92),T(2),1);
T(219) = getPowerDeriv(y(94)/y(93),T(2),1);
T(220) = getPowerDeriv(y(95)/y(94),T(2),1);
T(221) = getPowerDeriv(y(96)/y(95),T(2),1);
T(222) = getPowerDeriv(y(97)/y(96),T(2),1);
T(223) = getPowerDeriv(y(98)/y(97),T(2),1);
T(224) = getPowerDeriv(y(99)/y(98),T(2),1);
T(225) = getPowerDeriv(y(100)/y(99),T(2),1);
T(226) = getPowerDeriv(y(101)/y(100),T(2),1);
T(227) = getPowerDeriv(y(102)/y(101),T(2),1);
T(228) = getPowerDeriv(y(103)/y(102),T(2),1);
T(229) = getPowerDeriv(y(104)/y(103),T(2),1);
T(230) = getPowerDeriv(y(105)/y(104),T(2),1);
T(231) = getPowerDeriv(y(106)/y(105),T(2),1);
T(232) = getPowerDeriv(y(107)/y(106),T(2),1);
T(233) = getPowerDeriv(y(108)/y(107),T(2),1);
T(234) = getPowerDeriv(y(109)/y(108),T(2),1);
T(235) = getPowerDeriv(y(110)/y(109),T(2),1);
T(236) = getPowerDeriv(y(111)/y(110),T(2),1);
T(237) = getPowerDeriv(y(112)/y(111),T(2),1);
T(238) = (1-params(1))*x(5)*getPowerDeriv(y(276)*x(5),T(171),1);
T(239) = getPowerDeriv(T(172),1/(params(5)-1),1);
T(240) = y(271)*T(238)*T(239);
T(241) = getPowerDeriv(T(172),params(5)/(params(5)-1),1);
T(242) = params(1)*params(48)*getPowerDeriv(params(48)*y(278),T(171),1);
T(243) = y(271)*T(239)*T(242);

end
