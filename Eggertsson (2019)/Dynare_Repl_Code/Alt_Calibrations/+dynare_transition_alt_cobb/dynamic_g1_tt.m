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

assert(length(T) >= 292);

T = dynare_transition_alt_cobb.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(234) = getPowerDeriv(y(707)/y(288),T(2),1);
T(235) = getPowerDeriv(y(708)/y(289),T(2),1);
T(236) = getPowerDeriv(y(709)/y(290),T(2),1);
T(237) = getPowerDeriv(y(710)/y(291),T(2),1);
T(238) = getPowerDeriv(y(711)/y(292),T(2),1);
T(239) = getPowerDeriv(y(712)/y(293),T(2),1);
T(240) = getPowerDeriv(y(713)/y(294),T(2),1);
T(241) = getPowerDeriv(y(714)/y(295),T(2),1);
T(242) = getPowerDeriv(y(715)/y(296),T(2),1);
T(243) = getPowerDeriv(y(716)/y(297),T(2),1);
T(244) = getPowerDeriv(y(717)/y(298),T(2),1);
T(245) = getPowerDeriv(y(718)/y(299),T(2),1);
T(246) = getPowerDeriv(y(719)/y(300),T(2),1);
T(247) = getPowerDeriv(y(720)/y(301),T(2),1);
T(248) = getPowerDeriv(y(721)/y(302),T(2),1);
T(249) = getPowerDeriv(y(722)/y(303),T(2),1);
T(250) = getPowerDeriv(y(723)/y(304),T(2),1);
T(251) = getPowerDeriv(y(724)/y(305),T(2),1);
T(252) = getPowerDeriv(y(725)/y(306),T(2),1);
T(253) = getPowerDeriv(y(726)/y(307),T(2),1);
T(254) = getPowerDeriv(y(727)/y(308),T(2),1);
T(255) = getPowerDeriv(y(728)/y(309),T(2),1);
T(256) = getPowerDeriv(y(729)/y(310),T(2),1);
T(257) = getPowerDeriv(y(730)/y(311),T(2),1);
T(258) = getPowerDeriv(y(731)/y(312),T(2),1);
T(259) = getPowerDeriv(y(732)/y(313),T(2),1);
T(260) = getPowerDeriv(y(733)/y(314),T(2),1);
T(261) = getPowerDeriv(y(734)/y(315),T(2),1);
T(262) = getPowerDeriv(y(735)/y(316),T(2),1);
T(263) = getPowerDeriv(y(736)/y(317),T(2),1);
T(264) = getPowerDeriv(y(737)/y(318),T(2),1);
T(265) = getPowerDeriv(y(738)/y(319),T(2),1);
T(266) = getPowerDeriv(y(739)/y(320),T(2),1);
T(267) = getPowerDeriv(y(740)/y(321),T(2),1);
T(268) = getPowerDeriv(y(741)/y(322),T(2),1);
T(269) = getPowerDeriv(y(742)/y(323),T(2),1);
T(270) = getPowerDeriv(y(743)/y(324),T(2),1);
T(271) = getPowerDeriv(y(744)/y(325),T(2),1);
T(272) = getPowerDeriv(y(745)/y(326),T(2),1);
T(273) = getPowerDeriv(y(746)/y(327),T(2),1);
T(274) = getPowerDeriv(y(747)/y(328),T(2),1);
T(275) = getPowerDeriv(y(748)/y(329),T(2),1);
T(276) = getPowerDeriv(y(749)/y(330),T(2),1);
T(277) = getPowerDeriv(y(750)/y(331),T(2),1);
T(278) = getPowerDeriv(y(751)/y(332),T(2),1);
T(279) = getPowerDeriv(y(752)/y(333),T(2),1);
T(280) = getPowerDeriv(y(753)/y(334),T(2),1);
T(281) = getPowerDeriv(y(754)/y(335),T(2),1);
T(282) = getPowerDeriv(y(755)/y(336),T(2),1);
T(283) = getPowerDeriv(y(756)/y(337),T(2),1);
T(284) = getPowerDeriv(y(757)/y(338),T(2),1);
T(285) = getPowerDeriv(y(758)/y(339),T(2),1);
T(286) = getPowerDeriv(y(759)/y(340),T(2),1);
T(287) = getPowerDeriv(y(760)/y(341),T(2),1);
T(288) = getPowerDeriv(y(761)/y(342),T(2),1);
T(289) = getPowerDeriv(T(226),params(1),1);
T(290) = getPowerDeriv(T(226),params(1)-1,1);
T(291) = getPowerDeriv(y(507)*x(it_, 5),1-params(1),1);
T(292) = (-T(223));

end
