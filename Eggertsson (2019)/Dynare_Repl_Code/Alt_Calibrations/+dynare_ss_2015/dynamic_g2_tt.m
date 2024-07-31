function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 477);

T = dynare_ss_2015.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(414) = getPowerDeriv(y(58)/y(57),T(3),2);
T(415) = getPowerDeriv(y(59)/y(58),T(3),2);
T(416) = getPowerDeriv(y(60)/y(59),T(3),2);
T(417) = getPowerDeriv(y(61)/y(60),T(3),2);
T(418) = getPowerDeriv(y(62)/y(61),T(3),2);
T(419) = getPowerDeriv(y(63)/y(62),T(3),2);
T(420) = getPowerDeriv(y(64)/y(63),T(3),2);
T(421) = getPowerDeriv(y(65)/y(64),T(3),2);
T(422) = getPowerDeriv(y(66)/y(65),T(3),2);
T(423) = getPowerDeriv(y(67)/y(66),T(3),2);
T(424) = getPowerDeriv(y(68)/y(67),T(3),2);
T(425) = getPowerDeriv(y(69)/y(68),T(3),2);
T(426) = getPowerDeriv(y(70)/y(69),T(3),2);
T(427) = getPowerDeriv(y(71)/y(70),T(3),2);
T(428) = getPowerDeriv(y(72)/y(71),T(3),2);
T(429) = getPowerDeriv(y(73)/y(72),T(3),2);
T(430) = getPowerDeriv(y(74)/y(73),T(3),2);
T(431) = getPowerDeriv(y(75)/y(74),T(3),2);
T(432) = getPowerDeriv(y(76)/y(75),T(3),2);
T(433) = getPowerDeriv(y(77)/y(76),T(3),2);
T(434) = getPowerDeriv(y(78)/y(77),T(3),2);
T(435) = getPowerDeriv(y(79)/y(78),T(3),2);
T(436) = getPowerDeriv(y(80)/y(79),T(3),2);
T(437) = getPowerDeriv(y(81)/y(80),T(3),2);
T(438) = getPowerDeriv(y(82)/y(81),T(3),2);
T(439) = getPowerDeriv(y(83)/y(82),T(3),2);
T(440) = getPowerDeriv(y(84)/y(83),T(3),2);
T(441) = getPowerDeriv(y(85)/y(84),T(3),2);
T(442) = getPowerDeriv(y(86)/y(85),T(3),2);
T(443) = getPowerDeriv(y(87)/y(86),T(3),2);
T(444) = getPowerDeriv(y(88)/y(87),T(3),2);
T(445) = getPowerDeriv(y(89)/y(88),T(3),2);
T(446) = getPowerDeriv(y(90)/y(89),T(3),2);
T(447) = getPowerDeriv(y(91)/y(90),T(3),2);
T(448) = getPowerDeriv(y(92)/y(91),T(3),2);
T(449) = getPowerDeriv(y(93)/y(92),T(3),2);
T(450) = getPowerDeriv(y(94)/y(93),T(3),2);
T(451) = getPowerDeriv(y(95)/y(94),T(3),2);
T(452) = getPowerDeriv(y(96)/y(95),T(3),2);
T(453) = getPowerDeriv(y(97)/y(96),T(3),2);
T(454) = getPowerDeriv(y(98)/y(97),T(3),2);
T(455) = getPowerDeriv(y(99)/y(98),T(3),2);
T(456) = getPowerDeriv(y(100)/y(99),T(3),2);
T(457) = getPowerDeriv(y(101)/y(100),T(3),2);
T(458) = getPowerDeriv(y(102)/y(101),T(3),2);
T(459) = getPowerDeriv(y(103)/y(102),T(3),2);
T(460) = getPowerDeriv(y(104)/y(103),T(3),2);
T(461) = getPowerDeriv(y(105)/y(104),T(3),2);
T(462) = getPowerDeriv(y(106)/y(105),T(3),2);
T(463) = getPowerDeriv(y(107)/y(106),T(3),2);
T(464) = getPowerDeriv(y(108)/y(107),T(3),2);
T(465) = getPowerDeriv(y(109)/y(108),T(3),2);
T(466) = getPowerDeriv(y(110)/y(109),T(3),2);
T(467) = getPowerDeriv(y(111)/y(110),T(3),2);
T(468) = getPowerDeriv(y(112)/y(111),T(3),2);
T(469) = y(276)+y(276);
T(470) = T(405)*T(405);
T(471) = (1-params(1))*params(223)*params(223)*getPowerDeriv(y(276)*params(223),T(214),2);
T(472) = getPowerDeriv(T(215),1/(params(5)-1),2);
T(473) = y(271)*(T(407)*T(471)+T(406)*T(406)*T(472));
T(474) = params(1)*params(224)*params(224)*getPowerDeriv(params(224)*y(278),T(214),2);
T(475) = y(271)*(T(411)*T(411)*T(472)+T(407)*T(474));
T(476) = getPowerDeriv(T(215),params(5)/(params(5)-1),2);
T(477) = (y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(270)+params(3))+(y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(270)+params(3));

end
