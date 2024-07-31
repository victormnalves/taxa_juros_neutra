function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = dynare_transition_alt.static_resid_tt(T, y, x, params);
end
residual = zeros(475, 1);
lhs = y(1);
rhs = x(1)*y(25)*1/x(142);
residual(1) = lhs - rhs;
lhs = y(2);
rhs = y(1)*x(6);
residual(2) = lhs - rhs;
lhs = y(3);
rhs = y(2)*x(7);
residual(3) = lhs - rhs;
lhs = y(4);
rhs = y(3)*x(8);
residual(4) = lhs - rhs;
lhs = y(5);
rhs = y(4)*x(9);
residual(5) = lhs - rhs;
lhs = y(6);
rhs = y(5)*x(10);
residual(6) = lhs - rhs;
lhs = y(7);
rhs = y(6)*x(11);
residual(7) = lhs - rhs;
lhs = y(8);
rhs = y(7)*x(12);
residual(8) = lhs - rhs;
lhs = y(9);
rhs = y(8)*x(13);
residual(9) = lhs - rhs;
lhs = y(10);
rhs = y(9)*x(14);
residual(10) = lhs - rhs;
lhs = y(11);
rhs = y(10)*x(15);
residual(11) = lhs - rhs;
lhs = y(12);
rhs = y(11)*x(16);
residual(12) = lhs - rhs;
lhs = y(13);
rhs = y(12)*x(17);
residual(13) = lhs - rhs;
lhs = y(14);
rhs = y(13)*x(18);
residual(14) = lhs - rhs;
lhs = y(15);
rhs = y(14)*x(19);
residual(15) = lhs - rhs;
lhs = y(16);
rhs = y(15)*x(20);
residual(16) = lhs - rhs;
lhs = y(17);
rhs = y(16)*x(21);
residual(17) = lhs - rhs;
lhs = y(18);
rhs = y(17)*x(22);
residual(18) = lhs - rhs;
lhs = y(19);
rhs = y(18)*x(23);
residual(19) = lhs - rhs;
lhs = y(20);
rhs = y(19)*x(24);
residual(20) = lhs - rhs;
lhs = y(21);
rhs = y(20)*x(25);
residual(21) = lhs - rhs;
lhs = y(22);
rhs = y(21)*x(26);
residual(22) = lhs - rhs;
lhs = y(23);
rhs = y(22)*x(27);
residual(23) = lhs - rhs;
lhs = y(24);
rhs = y(23)*x(28);
residual(24) = lhs - rhs;
lhs = y(25);
rhs = y(24)*x(29);
residual(25) = lhs - rhs;
lhs = y(26);
rhs = y(25)*x(30);
residual(26) = lhs - rhs;
lhs = y(27);
rhs = y(26)*x(31);
residual(27) = lhs - rhs;
lhs = y(28);
rhs = y(27)*x(32);
residual(28) = lhs - rhs;
lhs = y(29);
rhs = y(28)*x(33);
residual(29) = lhs - rhs;
lhs = y(30);
rhs = y(29)*x(34);
residual(30) = lhs - rhs;
lhs = y(31);
rhs = y(30)*x(35);
residual(31) = lhs - rhs;
lhs = y(32);
rhs = y(31)*x(36);
residual(32) = lhs - rhs;
lhs = y(33);
rhs = y(32)*x(37);
residual(33) = lhs - rhs;
lhs = y(34);
rhs = y(33)*x(38);
residual(34) = lhs - rhs;
lhs = y(35);
rhs = y(34)*x(39);
residual(35) = lhs - rhs;
lhs = y(36);
rhs = y(35)*x(40);
residual(36) = lhs - rhs;
lhs = y(37);
rhs = y(36)*x(41);
residual(37) = lhs - rhs;
lhs = y(38);
rhs = y(37)*x(42);
residual(38) = lhs - rhs;
lhs = y(39);
rhs = y(38)*x(43);
residual(39) = lhs - rhs;
lhs = y(40);
rhs = y(39)*x(44);
residual(40) = lhs - rhs;
lhs = y(41);
rhs = y(40)*x(45);
residual(41) = lhs - rhs;
lhs = y(42);
rhs = y(41)*x(46);
residual(42) = lhs - rhs;
lhs = y(43);
rhs = y(42)*x(47);
residual(43) = lhs - rhs;
lhs = y(44);
rhs = y(43)*x(48);
residual(44) = lhs - rhs;
lhs = y(45);
rhs = y(44)*x(49);
residual(45) = lhs - rhs;
lhs = y(46);
rhs = y(45)*x(50);
residual(46) = lhs - rhs;
lhs = y(47);
rhs = y(46)*x(51);
residual(47) = lhs - rhs;
lhs = y(48);
rhs = y(47)*x(52);
residual(48) = lhs - rhs;
lhs = y(49);
rhs = y(48)*x(53);
residual(49) = lhs - rhs;
lhs = y(50);
rhs = y(49)*x(54);
residual(50) = lhs - rhs;
lhs = y(51);
rhs = y(50)*x(55);
residual(51) = lhs - rhs;
lhs = y(52);
rhs = y(51)*x(56);
residual(52) = lhs - rhs;
lhs = y(53);
rhs = y(52)*x(57);
residual(53) = lhs - rhs;
lhs = y(54);
rhs = y(53)*x(58);
residual(54) = lhs - rhs;
lhs = y(55);
rhs = y(54)*x(59);
residual(55) = lhs - rhs;
lhs = y(56);
rhs = y(55)*x(60);
residual(56) = lhs - rhs;
lhs = T(1);
rhs = T(3)*T(4)+y(171)*T(6)*x(63)/(params(2)*x(118));
residual(57) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(7)+y(172)*T(8)*x(64)/T(9);
residual(58) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(10)+y(173)*T(11)*x(65)/T(12);
residual(59) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(13)+y(174)*T(14)*x(66)/T(15);
residual(60) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(16)+y(175)*T(17)*x(67)/T(18);
residual(61) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(19)+y(176)*T(20)*x(68)/T(21);
residual(62) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(22)+y(177)*T(23)*x(69)/T(24);
residual(63) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(25)+y(178)*T(26)*x(70)/T(27);
residual(64) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(28)+y(179)*T(29)*x(71)/T(30);
residual(65) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(31)+y(180)*T(32)*x(72)/T(33);
residual(66) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(34)+y(181)*T(35)*x(73)/T(36);
residual(67) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(37)+y(182)*T(38)*x(74)/T(39);
residual(68) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(40)+y(183)*T(41)*x(75)/T(42);
residual(69) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(43)+y(184)*T(44)*x(76)/T(45);
residual(70) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(46)+y(185)*T(47)*x(77)/T(48);
residual(71) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(49)+y(186)*T(50)*x(78)/T(51);
residual(72) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(52)+y(187)*T(53)*x(79)/T(54);
residual(73) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(55)+y(188)*T(56)*x(80)/T(57);
residual(74) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(58)+y(189)*T(59)*x(81)/T(60);
residual(75) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(61)+y(190)*T(62)*x(82)/T(63);
residual(76) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(64)+y(191)*T(65)*x(83)/T(66);
residual(77) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(67)+y(192)*T(68)*x(84)/T(69);
residual(78) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(70)+y(193)*T(71)*x(85)/T(72);
residual(79) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(73)+y(194)*T(74)*x(86)/T(75);
residual(80) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(76)+y(195)*T(77)*x(87)/T(78);
residual(81) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(79)+y(196)*T(80)*x(88)/T(81);
residual(82) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(82)+y(197)*T(83)*x(89)/T(84);
residual(83) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(85)+y(198)*T(86)*x(90)/T(87);
residual(84) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(88)+y(199)*T(89)*x(91)/T(90);
residual(85) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(91)+y(200)*T(92)*x(92)/T(93);
residual(86) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(94)+y(201)*T(95)*x(93)/T(96);
residual(87) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(97)+y(202)*T(98)*x(94)/T(99);
residual(88) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(100)+y(203)*T(101)*x(95)/T(102);
residual(89) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(103)+y(204)*T(104)*x(96)/T(105);
residual(90) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(106)+y(205)*T(107)*x(97)/T(108);
residual(91) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(109)+y(206)*T(110)*x(98)/T(111);
residual(92) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(112)+y(207)*T(113)*x(99)/T(114);
residual(93) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(115)+y(208)*T(116)*x(100)/T(117);
residual(94) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(118)+y(209)*T(119)*x(101)/T(120);
residual(95) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(121)+y(210)*T(122)*x(102)/T(123);
residual(96) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(124)+y(211)*T(125)*x(103)/T(126);
residual(97) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(127)+y(212)*T(128)*x(104)/T(129);
residual(98) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(130)+y(213)*T(131)*x(105)/T(132);
residual(99) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(133)+y(214)*T(134)*x(106)/T(135);
residual(100) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(136)+y(215)*T(137)*x(107)/T(138);
residual(101) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(139)+y(216)*T(140)*x(108)/T(141);
residual(102) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(142)+y(217)*T(143)*x(109)/T(144);
residual(103) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(145)+y(218)*T(146)*x(110)/T(147);
residual(104) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(148)+y(219)*T(149)*x(111)/T(150);
residual(105) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(151)+y(220)*T(152)*x(112)/T(153);
residual(106) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(154)+y(221)*T(155)*x(113)/T(156);
residual(107) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(157)+y(222)*T(158)*x(114)/T(159);
residual(108) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(160)+y(223)*T(161)*x(115)/T(162);
residual(109) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(163)+y(224)*T(164)*x(116)/T(165);
residual(110) = lhs - rhs;
lhs = T(1);
rhs = T(4)*T(166)+y(225)*T(167)*x(117)/T(168);
residual(111) = lhs - rhs;
lhs = y(266);
rhs = y(112)*T(181);
residual(112) = lhs - rhs;
residual(113) = y(113);
lhs = y(114);
rhs = y(226)+T(170)*params(8)-y(57);
residual(114) = lhs - rhs;
lhs = y(115);
rhs = y(114)*T(4)/x(63)+y(227)+T(170)*params(9)-y(58);
residual(115) = lhs - rhs;
lhs = y(116);
rhs = y(115)*T(4)/x(64)+y(228)+T(170)*params(10)-y(59);
residual(116) = lhs - rhs;
lhs = y(117);
rhs = y(116)*T(4)/x(65)+y(229)+T(170)*params(11)-y(60);
residual(117) = lhs - rhs;
lhs = y(118);
rhs = y(117)*T(4)/x(66)+y(230)+T(170)*params(12)-y(61);
residual(118) = lhs - rhs;
lhs = y(119);
rhs = y(118)*T(4)/x(67)+y(231)+T(170)*params(13)-y(62);
residual(119) = lhs - rhs;
lhs = y(120);
rhs = y(119)*T(4)/x(68)+y(232)+T(170)*params(14)-y(63);
residual(120) = lhs - rhs;
lhs = y(121);
rhs = y(120)*T(4)/x(69)+y(233)+T(170)*params(15)-y(64);
residual(121) = lhs - rhs;
lhs = y(122);
rhs = y(121)*T(4)/x(70)+y(234)+T(170)*params(16)-y(65);
residual(122) = lhs - rhs;
lhs = y(123);
rhs = y(122)*T(4)/x(71)+y(235)+T(170)*params(17)-y(66);
residual(123) = lhs - rhs;
lhs = y(124);
rhs = y(123)*T(4)/x(72)+y(236)+T(170)*params(18)-y(67);
residual(124) = lhs - rhs;
lhs = y(125);
rhs = y(124)*T(4)/x(73)+y(237)+T(170)*params(19)-y(68);
residual(125) = lhs - rhs;
lhs = y(126);
rhs = y(125)*T(4)/x(74)+y(238)+T(170)*params(20)-y(69);
residual(126) = lhs - rhs;
lhs = y(127);
rhs = y(126)*T(4)/x(75)+y(239)+T(170)*params(21)-y(70);
residual(127) = lhs - rhs;
lhs = y(128);
rhs = y(127)*T(4)/x(76)+y(240)+T(170)*params(22)-y(71);
residual(128) = lhs - rhs;
lhs = y(129);
rhs = y(128)*T(4)/x(77)+y(241)+T(170)*params(23)-y(72);
residual(129) = lhs - rhs;
lhs = y(130);
rhs = y(129)*T(4)/x(78)+y(242)+T(170)*params(24)-y(73);
residual(130) = lhs - rhs;
lhs = y(131);
rhs = y(130)*T(4)/x(79)+y(243)+T(170)*params(25)-y(74);
residual(131) = lhs - rhs;
lhs = y(132);
rhs = y(131)*T(4)/x(80)+y(244)+T(170)*params(26)-y(75);
residual(132) = lhs - rhs;
lhs = y(133);
rhs = y(132)*T(4)/x(81)+y(245)+T(170)*params(27)-y(76);
residual(133) = lhs - rhs;
lhs = y(134);
rhs = y(133)*T(4)/x(82)+y(246)+T(170)*params(28)-y(77);
residual(134) = lhs - rhs;
lhs = y(135);
rhs = y(134)*T(4)/x(83)+y(247)+T(170)*params(29)-y(78);
residual(135) = lhs - rhs;
lhs = y(136);
rhs = y(135)*T(4)/x(84)+y(248)+T(170)*params(30)-y(79);
residual(136) = lhs - rhs;
lhs = y(137);
rhs = y(136)*T(4)/x(85)+y(249)+T(170)*params(31)-y(80);
residual(137) = lhs - rhs;
lhs = y(138);
rhs = y(137)*T(4)/x(86)+y(250)+T(170)*params(32)-y(81);
residual(138) = lhs - rhs;
lhs = y(139);
rhs = y(138)*T(4)/x(87)+y(251)+T(170)*params(33)-y(82);
residual(139) = lhs - rhs;
lhs = y(140);
rhs = y(139)*T(4)/x(88)+y(252)+T(170)*params(34)-y(83);
residual(140) = lhs - rhs;
lhs = y(141);
rhs = y(140)*T(4)/x(89)+y(253)+T(170)*params(35)-y(84);
residual(141) = lhs - rhs;
lhs = y(142);
rhs = y(141)*T(4)/x(90)+y(254)+T(170)*params(36)-y(85);
residual(142) = lhs - rhs;
lhs = y(143);
rhs = y(142)*T(4)/x(91)+y(255)+T(170)*params(37)-y(86);
residual(143) = lhs - rhs;
lhs = y(144);
rhs = y(143)*T(4)/x(92)+y(256)+T(170)*params(38)+y(267)-y(87);
residual(144) = lhs - rhs;
lhs = y(267);
rhs = y(56)*x(1)*y(266)/y(32);
residual(145) = lhs - rhs;
lhs = y(145);
rhs = y(144)*T(4)/x(93)+y(257)+T(170)*params(39)-y(88);
residual(146) = lhs - rhs;
lhs = y(146);
rhs = y(145)*T(4)/x(94)+y(258)+T(170)*params(40)-y(89);
residual(147) = lhs - rhs;
lhs = y(147);
rhs = y(146)*T(4)/x(95)+y(259)+T(170)*params(41)-y(90);
residual(148) = lhs - rhs;
lhs = y(148);
rhs = y(147)*T(4)/x(96)+y(260)+T(170)*params(42)-y(91);
residual(149) = lhs - rhs;
lhs = y(149);
rhs = y(148)*T(4)/x(97)+y(261)+T(170)*params(43)-y(92);
residual(150) = lhs - rhs;
lhs = y(150);
rhs = y(149)*T(4)/x(98)+y(262)+T(170)*params(44)-y(93);
residual(151) = lhs - rhs;
lhs = y(151);
rhs = y(150)*T(4)/x(99)+y(263)+T(170)*params(45)-y(94);
residual(152) = lhs - rhs;
lhs = y(152);
rhs = y(151)*T(4)/x(100)+y(264)+T(170)*params(46)-y(95);
residual(153) = lhs - rhs;
lhs = y(153);
rhs = y(152)*T(4)/x(101)+y(265)+T(170)*params(47)-y(96);
residual(154) = lhs - rhs;
lhs = y(154);
rhs = y(153)*T(4)/x(102)-y(97);
residual(155) = lhs - rhs;
lhs = y(155);
rhs = y(154)*T(4)/x(103)-y(98);
residual(156) = lhs - rhs;
lhs = y(156);
rhs = y(155)*T(4)/x(104)-y(99);
residual(157) = lhs - rhs;
lhs = y(157);
rhs = y(156)*T(4)/x(105)-y(100);
residual(158) = lhs - rhs;
lhs = y(158);
rhs = y(157)*T(4)/x(106)-y(101);
residual(159) = lhs - rhs;
lhs = y(159);
rhs = y(158)*T(4)/x(107)-y(102);
residual(160) = lhs - rhs;
lhs = y(160);
rhs = y(159)*T(4)/x(108)-y(103);
residual(161) = lhs - rhs;
lhs = y(161);
rhs = y(160)*T(4)/x(109)-y(104);
residual(162) = lhs - rhs;
lhs = y(162);
rhs = y(161)*T(4)/x(110)-y(105);
residual(163) = lhs - rhs;
lhs = y(163);
rhs = y(162)*T(4)/x(111)-y(106);
residual(164) = lhs - rhs;
lhs = y(164);
rhs = y(163)*T(4)/x(112)-y(107);
residual(165) = lhs - rhs;
lhs = y(165);
rhs = y(164)*T(4)/x(113)-y(108);
residual(166) = lhs - rhs;
lhs = y(166);
rhs = y(165)*T(4)/x(114)-y(109);
residual(167) = lhs - rhs;
lhs = y(167);
rhs = y(166)*T(4)/x(115)-y(110);
residual(168) = lhs - rhs;
lhs = y(168);
rhs = y(167)*T(4)/x(116)-y(111);
residual(169) = lhs - rhs;
lhs = 0;
rhs = y(168)*T(4)/x(117)-x(1)*y(266)-y(112);
residual(170) = lhs - rhs;
residual(171) = y(169);
residual(172) = min(y(170),y(282)*x(174));
residual(173) = min(y(171),y(114)+y(283)*x(175));
residual(174) = min(y(172),y(115)+y(284)*x(176));
residual(175) = min(y(173),y(116)+y(285)*x(177));
residual(176) = min(y(174),y(117)+y(286)*x(178));
residual(177) = min(y(175),y(118)+y(287)*x(179));
residual(178) = min(y(176),y(119)+y(288)*x(180));
residual(179) = min(y(177),y(120)+y(289)*x(181));
residual(180) = min(y(178),y(121)+y(290)*x(182));
residual(181) = min(y(179),y(122)+y(291)*x(183));
residual(182) = min(y(180),y(123)+y(292)*x(184));
residual(183) = min(y(181),y(124)+y(293)*x(185));
residual(184) = min(y(182),y(125)+y(294)*x(186));
residual(185) = min(y(183),y(126)+y(295)*x(187));
residual(186) = min(y(184),y(127)+y(296)*x(188));
residual(187) = min(y(185),y(128)+y(297)*x(189));
residual(188) = min(y(186),y(129)+y(298)*x(190));
residual(189) = min(y(187),y(130)+y(299)*x(191));
residual(190) = min(y(188),y(131)+y(300)*x(192));
residual(191) = min(y(189),y(132)+y(301)*x(193));
residual(192) = min(y(190),y(133)+y(302)*x(194));
residual(193) = min(y(191),y(134)+y(303)*x(195));
residual(194) = min(y(192),y(135)+y(304)*x(196));
residual(195) = min(y(193),y(136)+y(305)*x(197));
residual(196) = min(y(194),y(137)+y(306)*x(198));
residual(197) = min(y(195),y(138)+y(307)*x(199));
residual(198) = min(y(196),y(139)+y(308)*x(200));
residual(199) = min(y(197),y(140)+y(309)*x(201));
residual(200) = min(y(198),y(141)+y(310)*x(202));
residual(201) = min(y(199),y(142)+y(311)*x(203));
residual(202) = min(y(200),y(143)+y(312)*x(204));
residual(203) = min(y(201),y(144)+y(313)*x(205));
residual(204) = min(y(202),y(145)+y(314)*x(206));
residual(205) = min(y(203),y(146)+y(315)*x(207));
residual(206) = min(y(204),y(147)+y(316)*x(208));
residual(207) = min(y(205),y(148)+y(317)*x(209));
residual(208) = min(y(206),y(149)+y(318)*x(210));
residual(209) = min(y(207),y(150)+y(319)*x(211));
residual(210) = min(y(208),y(151)+y(320)*x(212));
residual(211) = min(y(209),y(152)+y(321)*x(213));
residual(212) = min(y(210),y(153));
residual(213) = min(y(211),y(154));
residual(214) = min(y(212),y(155));
residual(215) = min(y(213),y(156));
residual(216) = min(y(214),y(157));
residual(217) = min(y(215),y(158));
residual(218) = min(y(216),y(159));
residual(219) = min(y(217),y(160));
residual(220) = min(y(218),y(161));
residual(221) = min(y(219),y(162));
residual(222) = min(y(220),y(163));
residual(223) = min(y(221),y(164));
residual(224) = min(y(222),y(165));
residual(225) = min(y(223),y(166));
residual(226) = min(y(224),y(167));
residual(227) = min(y(225),y(168));
lhs = y(282);
rhs = y(268)*params(8);
residual(228) = lhs - rhs;
lhs = y(283);
rhs = y(268)*params(9);
residual(229) = lhs - rhs;
lhs = y(284);
rhs = y(268)*params(10);
residual(230) = lhs - rhs;
lhs = y(285);
rhs = y(268)*params(11);
residual(231) = lhs - rhs;
lhs = y(286);
rhs = y(268)*params(12);
residual(232) = lhs - rhs;
lhs = y(287);
rhs = y(268)*params(13);
residual(233) = lhs - rhs;
lhs = y(288);
rhs = y(268)*params(14);
residual(234) = lhs - rhs;
lhs = y(289);
rhs = y(268)*params(15);
residual(235) = lhs - rhs;
lhs = y(290);
rhs = y(268)*params(16);
residual(236) = lhs - rhs;
lhs = y(291);
rhs = y(268)*params(17);
residual(237) = lhs - rhs;
lhs = y(292);
rhs = y(268)*params(18);
residual(238) = lhs - rhs;
lhs = y(293);
rhs = y(268)*params(19);
residual(239) = lhs - rhs;
lhs = y(294);
rhs = y(268)*params(20);
residual(240) = lhs - rhs;
lhs = y(295);
rhs = y(268)*params(21);
residual(241) = lhs - rhs;
lhs = y(296);
rhs = y(268)*params(22);
residual(242) = lhs - rhs;
lhs = y(297);
rhs = y(268)*params(23);
residual(243) = lhs - rhs;
lhs = y(298);
rhs = y(268)*params(24);
residual(244) = lhs - rhs;
lhs = y(299);
rhs = y(268)*params(25);
residual(245) = lhs - rhs;
lhs = y(300);
rhs = y(268)*params(26);
residual(246) = lhs - rhs;
lhs = y(301);
rhs = y(268)*params(27);
residual(247) = lhs - rhs;
lhs = y(302);
rhs = y(268)*params(28);
residual(248) = lhs - rhs;
lhs = y(303);
rhs = y(268)*params(29);
residual(249) = lhs - rhs;
lhs = y(304);
rhs = y(268)*params(30);
residual(250) = lhs - rhs;
lhs = y(305);
rhs = y(268)*params(31);
residual(251) = lhs - rhs;
lhs = y(306);
rhs = y(268)*params(32);
residual(252) = lhs - rhs;
lhs = y(307);
rhs = y(268)*params(33);
residual(253) = lhs - rhs;
lhs = y(308);
rhs = y(268)*params(34);
residual(254) = lhs - rhs;
lhs = y(309);
rhs = y(268)*params(35);
residual(255) = lhs - rhs;
lhs = y(310);
rhs = y(268)*params(36);
residual(256) = lhs - rhs;
lhs = y(311);
rhs = y(268)*params(37);
residual(257) = lhs - rhs;
lhs = y(312);
rhs = y(268)*params(38);
residual(258) = lhs - rhs;
lhs = y(313);
rhs = y(268)*params(39);
residual(259) = lhs - rhs;
lhs = y(314);
rhs = y(268)*params(40);
residual(260) = lhs - rhs;
lhs = y(315);
rhs = y(268)*params(41);
residual(261) = lhs - rhs;
lhs = y(316);
rhs = y(268)*params(42);
residual(262) = lhs - rhs;
lhs = y(317);
rhs = y(268)*params(43);
residual(263) = lhs - rhs;
lhs = y(318);
rhs = y(268)*params(44);
residual(264) = lhs - rhs;
lhs = y(319);
rhs = y(268)*params(45);
residual(265) = lhs - rhs;
lhs = y(320);
rhs = y(268)*params(46);
residual(266) = lhs - rhs;
lhs = y(321);
rhs = y(268)*params(47);
residual(267) = lhs - rhs;
lhs = y(226);
rhs = params(8)*y(273)/y(276);
residual(268) = lhs - rhs;
lhs = y(227);
rhs = params(9)*y(273)/y(276);
residual(269) = lhs - rhs;
lhs = y(228);
rhs = params(10)*y(273)/y(276);
residual(270) = lhs - rhs;
lhs = y(229);
rhs = params(11)*y(273)/y(276);
residual(271) = lhs - rhs;
lhs = y(230);
rhs = params(12)*y(273)/y(276);
residual(272) = lhs - rhs;
lhs = y(231);
rhs = params(13)*y(273)/y(276);
residual(273) = lhs - rhs;
lhs = y(232);
rhs = params(14)*y(273)/y(276);
residual(274) = lhs - rhs;
lhs = y(233);
rhs = params(15)*y(273)/y(276);
residual(275) = lhs - rhs;
lhs = y(234);
rhs = params(16)*y(273)/y(276);
residual(276) = lhs - rhs;
lhs = y(235);
rhs = params(17)*y(273)/y(276);
residual(277) = lhs - rhs;
lhs = y(236);
rhs = params(18)*y(273)/y(276);
residual(278) = lhs - rhs;
lhs = y(237);
rhs = params(19)*y(273)/y(276);
residual(279) = lhs - rhs;
lhs = y(238);
rhs = params(20)*y(273)/y(276);
residual(280) = lhs - rhs;
lhs = y(239);
rhs = params(21)*y(273)/y(276);
residual(281) = lhs - rhs;
lhs = y(240);
rhs = params(22)*y(273)/y(276);
residual(282) = lhs - rhs;
lhs = y(241);
rhs = params(23)*y(273)/y(276);
residual(283) = lhs - rhs;
lhs = y(242);
rhs = params(24)*y(273)/y(276);
residual(284) = lhs - rhs;
lhs = y(243);
rhs = params(25)*y(273)/y(276);
residual(285) = lhs - rhs;
lhs = y(244);
rhs = params(26)*y(273)/y(276);
residual(286) = lhs - rhs;
lhs = y(245);
rhs = params(27)*y(273)/y(276);
residual(287) = lhs - rhs;
lhs = y(246);
rhs = params(28)*y(273)/y(276);
residual(288) = lhs - rhs;
lhs = y(247);
rhs = params(29)*y(273)/y(276);
residual(289) = lhs - rhs;
lhs = y(248);
rhs = params(30)*y(273)/y(276);
residual(290) = lhs - rhs;
lhs = y(249);
rhs = params(31)*y(273)/y(276);
residual(291) = lhs - rhs;
lhs = y(250);
rhs = params(32)*y(273)/y(276);
residual(292) = lhs - rhs;
lhs = y(251);
rhs = params(33)*y(273)/y(276);
residual(293) = lhs - rhs;
lhs = y(252);
rhs = params(34)*y(273)/y(276);
residual(294) = lhs - rhs;
lhs = y(253);
rhs = params(35)*y(273)/y(276);
residual(295) = lhs - rhs;
lhs = y(254);
rhs = params(36)*y(273)/y(276);
residual(296) = lhs - rhs;
lhs = y(255);
rhs = params(37)*y(273)/y(276);
residual(297) = lhs - rhs;
lhs = y(256);
rhs = params(38)*y(273)/y(276);
residual(298) = lhs - rhs;
lhs = y(257);
rhs = params(39)*y(273)/y(276);
residual(299) = lhs - rhs;
lhs = y(258);
rhs = params(40)*y(273)/y(276);
residual(300) = lhs - rhs;
lhs = y(259);
rhs = params(41)*y(273)/y(276);
residual(301) = lhs - rhs;
lhs = y(260);
rhs = params(42)*y(273)/y(276);
residual(302) = lhs - rhs;
lhs = y(261);
rhs = params(43)*y(273)/y(276);
residual(303) = lhs - rhs;
lhs = y(262);
rhs = params(44)*y(273)/y(276);
residual(304) = lhs - rhs;
lhs = y(263);
rhs = params(45)*y(273)/y(276);
residual(305) = lhs - rhs;
lhs = y(264);
rhs = params(46)*y(273)/y(276);
residual(306) = lhs - rhs;
lhs = y(265);
rhs = params(47)*y(273)/y(276);
residual(307) = lhs - rhs;
lhs = y(271);
rhs = (x(4)-1)/x(4);
residual(308) = lhs - rhs;
lhs = y(273);
rhs = y(274)/x(4);
residual(309) = lhs - rhs;
lhs = y(268);
rhs = (1-params(1))*T(174)*T(175)*T(177)/params(49);
residual(310) = lhs - rhs;
lhs = y(269);
rhs = params(1)*T(174)*T(178)*T(179)/params(49);
residual(311) = lhs - rhs;
lhs = y(274);
rhs = T(172)^(params(5)/(params(5)-1))/params(49);
residual(312) = lhs - rhs;
lhs = y(270);
rhs = 1-params(3)+y(269)/x(2)-1;
residual(313) = lhs - rhs;
lhs = y(278)*y(281);
rhs = T(4)*y(278)*y(281)+y(274)*params(7)-y(279)*(1-y(280));
residual(314) = lhs - rhs;
lhs = y(280);
rhs = (y(274)*x(3)-y(278)*y(281))/y(279);
residual(315) = lhs - rhs;
lhs = y(279);
rhs = y(274)*params(7)+y(278)*y(270)*y(281);
residual(316) = lhs - rhs;
lhs = y(276)*y(272)*y(268);
rhs = y(279)*(1-y(280));
residual(317) = lhs - rhs;
lhs = y(275);
rhs = y(56)+y(55)+y(54)+y(53)+y(52)+y(51)+y(50)+y(49)+y(48)+y(47)+y(46)+y(45)+y(44)+y(43)+y(42)+y(41)+y(40)+y(39)+y(38)+y(37)+y(36)+y(35)+y(34)+y(33)+y(32)+y(31)+y(30)+y(29)+y(28)+y(27)+y(26)+y(25)+y(24)+y(23)+y(22)+y(21)+y(20)+y(19)+y(18)+y(17)+y(16)+y(15)+y(14)+y(13)+y(12)+y(11)+y(10)+y(9)+y(8)+y(7)+y(6)+y(5)+y(4)+y(3)+y(1)+y(2);
residual(318) = lhs - rhs;
lhs = y(276);
rhs = y(1)*params(8)+y(2)*params(9)+y(3)*params(10)+y(4)*params(11)+y(5)*params(12)+y(6)*params(13)+y(7)*params(14)+y(8)*params(15)+y(9)*params(16)+y(10)*params(17)+y(11)*params(18)+y(12)*params(19)+y(13)*params(20)+y(14)*params(21)+y(15)*params(22)+y(16)*params(23)+y(17)*params(24)+y(18)*params(25)+y(19)*params(26)+y(20)*params(27)+y(21)*params(28)+y(22)*params(29)+y(23)*params(30)+y(24)*params(31)+y(25)*params(32)+y(26)*params(33)+y(27)*params(34)+y(28)*params(35)+y(29)*params(36)+y(30)*params(37)+y(31)*params(38)+y(32)*params(39)+y(33)*params(40)+y(34)*params(41)+y(35)*params(42)+y(36)*params(43)+y(37)*params(44)+y(38)*params(45)+y(39)*params(46)+y(40)*params(47);
residual(319) = lhs - rhs;
lhs = y(277);
rhs = y(1)*y(57)+y(2)*y(58)+y(3)*y(59)+y(4)*y(60)+y(5)*y(61)+y(6)*y(62)+y(7)*y(63)+y(8)*y(64)+y(9)*y(65)+y(10)*y(66)+y(11)*y(67)+y(12)*y(68)+y(13)*y(69)+y(14)*y(70)+y(15)*y(71)+y(16)*y(72)+y(17)*y(73)+y(18)*y(74)+y(19)*y(75)+y(20)*y(76)+y(21)*y(77)+y(22)*y(78)+y(23)*y(79)+y(24)*y(80)+y(25)*y(81)+y(26)*y(82)+y(27)*y(83)+y(28)*y(84)+y(29)*y(85)+y(30)*y(86)+y(31)*y(87)+y(32)*y(88)+y(33)*y(89)+y(34)*y(90)+y(35)*y(91)+y(36)*y(92)+y(37)*y(93)+y(38)*y(94)+y(39)*y(95)+y(40)*y(96)+y(41)*y(97)+y(42)*y(98)+y(43)*y(99)+y(44)*y(100)+y(45)*y(101)+y(46)*y(102)+y(47)*y(103)+y(48)*y(104)+y(49)*y(105)+y(50)*y(106)+y(51)*y(107)+y(52)*y(108)+y(53)*y(109)+y(54)*y(110)+y(55)*y(111)+y(56)*y(112);
residual(320) = lhs - rhs;
lhs = y(278);
rhs = T(180)/T(182);
residual(321) = lhs - rhs;
lhs = y(322);
rhs = x(174);
residual(322) = lhs - rhs;
lhs = y(323);
rhs = x(175);
residual(323) = lhs - rhs;
lhs = y(324);
rhs = x(176);
residual(324) = lhs - rhs;
lhs = y(325);
rhs = x(177);
residual(325) = lhs - rhs;
lhs = y(326);
rhs = x(178);
residual(326) = lhs - rhs;
lhs = y(327);
rhs = x(179);
residual(327) = lhs - rhs;
lhs = y(328);
rhs = x(180);
residual(328) = lhs - rhs;
lhs = y(329);
rhs = x(181);
residual(329) = lhs - rhs;
lhs = y(330);
rhs = x(182);
residual(330) = lhs - rhs;
lhs = y(331);
rhs = x(183);
residual(331) = lhs - rhs;
lhs = y(332);
rhs = x(184);
residual(332) = lhs - rhs;
lhs = y(333);
rhs = x(185);
residual(333) = lhs - rhs;
lhs = y(334);
rhs = x(186);
residual(334) = lhs - rhs;
lhs = y(335);
rhs = x(187);
residual(335) = lhs - rhs;
lhs = y(336);
rhs = x(188);
residual(336) = lhs - rhs;
lhs = y(337);
rhs = x(189);
residual(337) = lhs - rhs;
lhs = y(338);
rhs = x(190);
residual(338) = lhs - rhs;
lhs = y(339);
rhs = x(191);
residual(339) = lhs - rhs;
lhs = y(340);
rhs = x(192);
residual(340) = lhs - rhs;
lhs = y(341);
rhs = x(193);
residual(341) = lhs - rhs;
lhs = y(342);
rhs = x(194);
residual(342) = lhs - rhs;
lhs = y(343);
rhs = x(195);
residual(343) = lhs - rhs;
lhs = y(344);
rhs = x(196);
residual(344) = lhs - rhs;
lhs = y(345);
rhs = x(197);
residual(345) = lhs - rhs;
lhs = y(346);
rhs = x(198);
residual(346) = lhs - rhs;
lhs = y(347);
rhs = x(199);
residual(347) = lhs - rhs;
lhs = y(348);
rhs = x(200);
residual(348) = lhs - rhs;
lhs = y(349);
rhs = x(201);
residual(349) = lhs - rhs;
lhs = y(350);
rhs = x(202);
residual(350) = lhs - rhs;
lhs = y(351);
rhs = x(203);
residual(351) = lhs - rhs;
lhs = y(352);
rhs = x(204);
residual(352) = lhs - rhs;
lhs = y(353);
rhs = x(205);
residual(353) = lhs - rhs;
lhs = y(354);
rhs = x(206);
residual(354) = lhs - rhs;
lhs = y(355);
rhs = x(207);
residual(355) = lhs - rhs;
lhs = y(356);
rhs = x(208);
residual(356) = lhs - rhs;
lhs = y(357);
rhs = x(209);
residual(357) = lhs - rhs;
lhs = y(358);
rhs = x(210);
residual(358) = lhs - rhs;
lhs = y(359);
rhs = x(211);
residual(359) = lhs - rhs;
lhs = y(360);
rhs = x(212);
residual(360) = lhs - rhs;
lhs = y(361);
rhs = x(213);
residual(361) = lhs - rhs;
lhs = y(362);
rhs = x(3);
residual(362) = lhs - rhs;
lhs = y(363);
rhs = x(142);
residual(363) = lhs - rhs;
lhs = y(364);
rhs = x(6);
residual(364) = lhs - rhs;
lhs = y(365);
rhs = x(7);
residual(365) = lhs - rhs;
lhs = y(366);
rhs = x(8);
residual(366) = lhs - rhs;
lhs = y(367);
rhs = x(9);
residual(367) = lhs - rhs;
lhs = y(368);
rhs = x(10);
residual(368) = lhs - rhs;
lhs = y(369);
rhs = x(11);
residual(369) = lhs - rhs;
lhs = y(370);
rhs = x(12);
residual(370) = lhs - rhs;
lhs = y(371);
rhs = x(13);
residual(371) = lhs - rhs;
lhs = y(372);
rhs = x(14);
residual(372) = lhs - rhs;
lhs = y(373);
rhs = x(15);
residual(373) = lhs - rhs;
lhs = y(374);
rhs = x(16);
residual(374) = lhs - rhs;
lhs = y(375);
rhs = x(17);
residual(375) = lhs - rhs;
lhs = y(376);
rhs = x(18);
residual(376) = lhs - rhs;
lhs = y(377);
rhs = x(19);
residual(377) = lhs - rhs;
lhs = y(378);
rhs = x(20);
residual(378) = lhs - rhs;
lhs = y(379);
rhs = x(21);
residual(379) = lhs - rhs;
lhs = y(380);
rhs = x(22);
residual(380) = lhs - rhs;
lhs = y(381);
rhs = x(23);
residual(381) = lhs - rhs;
lhs = y(382);
rhs = x(24);
residual(382) = lhs - rhs;
lhs = y(383);
rhs = x(25);
residual(383) = lhs - rhs;
lhs = y(384);
rhs = x(26);
residual(384) = lhs - rhs;
lhs = y(385);
rhs = x(27);
residual(385) = lhs - rhs;
lhs = y(386);
rhs = x(28);
residual(386) = lhs - rhs;
lhs = y(387);
rhs = x(29);
residual(387) = lhs - rhs;
lhs = y(388);
rhs = x(30);
residual(388) = lhs - rhs;
lhs = y(389);
rhs = x(31);
residual(389) = lhs - rhs;
lhs = y(390);
rhs = x(32);
residual(390) = lhs - rhs;
lhs = y(391);
rhs = x(33);
residual(391) = lhs - rhs;
lhs = y(392);
rhs = x(34);
residual(392) = lhs - rhs;
lhs = y(393);
rhs = x(35);
residual(393) = lhs - rhs;
lhs = y(394);
rhs = x(36);
residual(394) = lhs - rhs;
lhs = y(395);
rhs = x(37);
residual(395) = lhs - rhs;
lhs = y(396);
rhs = x(38);
residual(396) = lhs - rhs;
lhs = y(397);
rhs = x(39);
residual(397) = lhs - rhs;
lhs = y(398);
rhs = x(40);
residual(398) = lhs - rhs;
lhs = y(399);
rhs = x(41);
residual(399) = lhs - rhs;
lhs = y(400);
rhs = x(42);
residual(400) = lhs - rhs;
lhs = y(401);
rhs = x(43);
residual(401) = lhs - rhs;
lhs = y(402);
rhs = x(44);
residual(402) = lhs - rhs;
lhs = y(403);
rhs = x(45);
residual(403) = lhs - rhs;
lhs = y(404);
rhs = x(46);
residual(404) = lhs - rhs;
lhs = y(405);
rhs = x(47);
residual(405) = lhs - rhs;
lhs = y(406);
rhs = x(48);
residual(406) = lhs - rhs;
lhs = y(407);
rhs = x(49);
residual(407) = lhs - rhs;
lhs = y(408);
rhs = x(50);
residual(408) = lhs - rhs;
lhs = y(409);
rhs = x(51);
residual(409) = lhs - rhs;
lhs = y(410);
rhs = x(52);
residual(410) = lhs - rhs;
lhs = y(411);
rhs = x(53);
residual(411) = lhs - rhs;
lhs = y(412);
rhs = x(54);
residual(412) = lhs - rhs;
lhs = y(413);
rhs = x(55);
residual(413) = lhs - rhs;
lhs = y(414);
rhs = x(56);
residual(414) = lhs - rhs;
lhs = y(415);
rhs = x(57);
residual(415) = lhs - rhs;
lhs = y(416);
rhs = x(58);
residual(416) = lhs - rhs;
lhs = y(417);
rhs = x(59);
residual(417) = lhs - rhs;
lhs = y(418);
rhs = x(60);
residual(418) = lhs - rhs;
lhs = y(419);
rhs = x(1);
residual(419) = lhs - rhs;
lhs = y(420);
rhs = x(1);
residual(420) = lhs - rhs;
lhs = y(421);
rhs = x(1);
residual(421) = lhs - rhs;
lhs = y(422);
rhs = x(1);
residual(422) = lhs - rhs;
lhs = y(423);
rhs = x(1);
residual(423) = lhs - rhs;
lhs = y(424);
rhs = x(1);
residual(424) = lhs - rhs;
lhs = y(425);
rhs = x(1);
residual(425) = lhs - rhs;
lhs = y(426);
rhs = x(1);
residual(426) = lhs - rhs;
lhs = y(427);
rhs = x(1);
residual(427) = lhs - rhs;
lhs = y(428);
rhs = x(1);
residual(428) = lhs - rhs;
lhs = y(429);
rhs = x(1);
residual(429) = lhs - rhs;
lhs = y(430);
rhs = x(1);
residual(430) = lhs - rhs;
lhs = y(431);
rhs = x(1);
residual(431) = lhs - rhs;
lhs = y(432);
rhs = x(1);
residual(432) = lhs - rhs;
lhs = y(433);
rhs = x(1);
residual(433) = lhs - rhs;
lhs = y(434);
rhs = x(1);
residual(434) = lhs - rhs;
lhs = y(435);
rhs = x(1);
residual(435) = lhs - rhs;
lhs = y(436);
rhs = x(1);
residual(436) = lhs - rhs;
lhs = y(437);
rhs = x(1);
residual(437) = lhs - rhs;
lhs = y(438);
rhs = x(1);
residual(438) = lhs - rhs;
lhs = y(439);
rhs = x(1);
residual(439) = lhs - rhs;
lhs = y(440);
rhs = x(1);
residual(440) = lhs - rhs;
lhs = y(441);
rhs = x(1);
residual(441) = lhs - rhs;
lhs = y(442);
rhs = x(1);
residual(442) = lhs - rhs;
lhs = y(443);
rhs = x(1);
residual(443) = lhs - rhs;
lhs = y(444);
rhs = x(1);
residual(444) = lhs - rhs;
lhs = y(445);
rhs = x(1);
residual(445) = lhs - rhs;
lhs = y(446);
rhs = x(1);
residual(446) = lhs - rhs;
lhs = y(447);
rhs = x(1);
residual(447) = lhs - rhs;
lhs = y(448);
rhs = x(1);
residual(448) = lhs - rhs;
lhs = y(449);
rhs = x(1);
residual(449) = lhs - rhs;
lhs = y(450);
rhs = x(1);
residual(450) = lhs - rhs;
lhs = y(451);
rhs = x(1);
residual(451) = lhs - rhs;
lhs = y(452);
rhs = x(1);
residual(452) = lhs - rhs;
lhs = y(453);
rhs = x(1);
residual(453) = lhs - rhs;
lhs = y(454);
rhs = x(1);
residual(454) = lhs - rhs;
lhs = y(455);
rhs = x(1);
residual(455) = lhs - rhs;
lhs = y(456);
rhs = x(1);
residual(456) = lhs - rhs;
lhs = y(457);
rhs = x(1);
residual(457) = lhs - rhs;
lhs = y(458);
rhs = x(1);
residual(458) = lhs - rhs;
lhs = y(459);
rhs = x(1);
residual(459) = lhs - rhs;
lhs = y(460);
rhs = x(1);
residual(460) = lhs - rhs;
lhs = y(461);
rhs = x(1);
residual(461) = lhs - rhs;
lhs = y(462);
rhs = x(1);
residual(462) = lhs - rhs;
lhs = y(463);
rhs = x(1);
residual(463) = lhs - rhs;
lhs = y(464);
rhs = x(1);
residual(464) = lhs - rhs;
lhs = y(465);
rhs = x(1);
residual(465) = lhs - rhs;
lhs = y(466);
rhs = x(1);
residual(466) = lhs - rhs;
lhs = y(467);
rhs = x(1);
residual(467) = lhs - rhs;
lhs = y(468);
rhs = x(1);
residual(468) = lhs - rhs;
lhs = y(469);
rhs = x(1);
residual(469) = lhs - rhs;
lhs = y(470);
rhs = x(1);
residual(470) = lhs - rhs;
lhs = y(471);
rhs = x(1);
residual(471) = lhs - rhs;
lhs = y(472);
rhs = x(1);
residual(472) = lhs - rhs;
lhs = y(473);
rhs = x(1);
residual(473) = lhs - rhs;
lhs = y(474);
rhs = x(1);
residual(474) = lhs - rhs;
lhs = y(475);
rhs = x(2);
residual(475) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end