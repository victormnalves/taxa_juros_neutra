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

assert(length(T) >= 590);

T = dynare_ss_2015_alt.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(527) = getPowerDeriv(y(58)/y(57),T(3),2);
T(528) = getPowerDeriv(y(59)/y(58),T(3),2);
T(529) = getPowerDeriv(y(60)/y(59),T(3),2);
T(530) = getPowerDeriv(y(61)/y(60),T(3),2);
T(531) = getPowerDeriv(y(62)/y(61),T(3),2);
T(532) = getPowerDeriv(y(63)/y(62),T(3),2);
T(533) = getPowerDeriv(y(64)/y(63),T(3),2);
T(534) = getPowerDeriv(y(65)/y(64),T(3),2);
T(535) = getPowerDeriv(y(66)/y(65),T(3),2);
T(536) = getPowerDeriv(y(67)/y(66),T(3),2);
T(537) = getPowerDeriv(y(68)/y(67),T(3),2);
T(538) = getPowerDeriv(y(69)/y(68),T(3),2);
T(539) = getPowerDeriv(y(70)/y(69),T(3),2);
T(540) = getPowerDeriv(y(71)/y(70),T(3),2);
T(541) = getPowerDeriv(y(72)/y(71),T(3),2);
T(542) = getPowerDeriv(y(73)/y(72),T(3),2);
T(543) = getPowerDeriv(y(74)/y(73),T(3),2);
T(544) = getPowerDeriv(y(75)/y(74),T(3),2);
T(545) = getPowerDeriv(y(76)/y(75),T(3),2);
T(546) = getPowerDeriv(y(77)/y(76),T(3),2);
T(547) = getPowerDeriv(y(78)/y(77),T(3),2);
T(548) = getPowerDeriv(y(79)/y(78),T(3),2);
T(549) = getPowerDeriv(y(80)/y(79),T(3),2);
T(550) = getPowerDeriv(y(81)/y(80),T(3),2);
T(551) = getPowerDeriv(y(82)/y(81),T(3),2);
T(552) = getPowerDeriv(y(83)/y(82),T(3),2);
T(553) = getPowerDeriv(y(84)/y(83),T(3),2);
T(554) = getPowerDeriv(y(85)/y(84),T(3),2);
T(555) = getPowerDeriv(y(86)/y(85),T(3),2);
T(556) = getPowerDeriv(y(87)/y(86),T(3),2);
T(557) = getPowerDeriv(y(88)/y(87),T(3),2);
T(558) = getPowerDeriv(y(89)/y(88),T(3),2);
T(559) = getPowerDeriv(y(90)/y(89),T(3),2);
T(560) = getPowerDeriv(y(91)/y(90),T(3),2);
T(561) = getPowerDeriv(y(92)/y(91),T(3),2);
T(562) = getPowerDeriv(y(93)/y(92),T(3),2);
T(563) = getPowerDeriv(y(94)/y(93),T(3),2);
T(564) = getPowerDeriv(y(95)/y(94),T(3),2);
T(565) = getPowerDeriv(y(96)/y(95),T(3),2);
T(566) = getPowerDeriv(y(97)/y(96),T(3),2);
T(567) = getPowerDeriv(y(98)/y(97),T(3),2);
T(568) = getPowerDeriv(y(99)/y(98),T(3),2);
T(569) = getPowerDeriv(y(100)/y(99),T(3),2);
T(570) = getPowerDeriv(y(101)/y(100),T(3),2);
T(571) = getPowerDeriv(y(102)/y(101),T(3),2);
T(572) = getPowerDeriv(y(103)/y(102),T(3),2);
T(573) = getPowerDeriv(y(104)/y(103),T(3),2);
T(574) = getPowerDeriv(y(105)/y(104),T(3),2);
T(575) = getPowerDeriv(y(106)/y(105),T(3),2);
T(576) = getPowerDeriv(y(107)/y(106),T(3),2);
T(577) = getPowerDeriv(y(108)/y(107),T(3),2);
T(578) = getPowerDeriv(y(109)/y(108),T(3),2);
T(579) = getPowerDeriv(y(110)/y(109),T(3),2);
T(580) = getPowerDeriv(y(111)/y(110),T(3),2);
T(581) = getPowerDeriv(y(112)/y(111),T(3),2);
T(582) = y(276)+y(276);
T(583) = T(517)*T(517);
T(584) = (1-params(1))*params(223)*params(223)*getPowerDeriv(y(276)*params(223),T(214),2);
T(585) = getPowerDeriv(T(215),1/(params(5)-1),2);
T(586) = y(271)*(T(519)*T(584)+T(518)*T(518)*T(585));
T(587) = params(1)*params(224)*params(224)*getPowerDeriv(params(224)*y(278),T(214),2);
T(588) = y(271)*(T(523)*T(523)*T(585)+T(519)*T(587));
T(589) = getPowerDeriv(T(215),params(5)/(params(5)-1),2);
T(590) = (y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(270)+params(3))+(y(273)+y(276)+y(270)*y(278)+y(278)*params(3))*(y(270)+params(3));

end
