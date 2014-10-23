function tests = testMMkinetics
    tests = functiontests(localfunctions);
end

function testfunc1(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Elem,rt);
mathml = enzobj.toMathExpr;
calcvalue=mathml;
mathexp.mathexpr = 'Vm*S/(Km+S)';
mathexp.speciesnames = {'S'};
mathexp.rtconsts = {'Vm','Km'};
expvalue=mathexp;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc2(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Elem,rt);
ele2ems(enzobj,MMenKinetics.SMM);
calcvalue=enzobj.sumconsts;
enzobj.sumconsts.vm=0.6;
enzobj.sumconsts.km=0.5;
expvalue=enzobj.sumconsts;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc3(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Elem,rt);
ele2ems(enzobj,MMenKinetics.SMM);
calcvalue=num2str(simVel(enzobj,1));
expvalue=num2str(0.4);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc4(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'kfi',5,'kri',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.Noncomp,MMenKinetics.Elem,rt);
mathml = enzobj.toMathExpr;
calcvalue=mathml;
mathexp.sbmlmath = 'Vm/((1+Ci/Ki)*(1+Km/S))';
mathexp.speciesnames = {'S','Ci'};
mathexp.rtconsts = {'Vm','Km','Ki'};
expvalue=mathexp;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc5(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'kfi',5,'kri',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.Noncomp,MMenKinetics.Elem,rt);
ele2ems(enzobj,MMenKinetics.Noncomp);
calcvalue=enzobj.sumconsts;
enzobj.sumconsts.vm=0.6;
enzobj.sumconsts.km=0.5;
enzobj.sumconsts.ki=0.4;
expvalue=enzobj.sumconsts;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc6(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'kfi',5,'kri',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.Noncomp,MMenKinetics.Elem,rt);
ele2ems(enzobj,MMenKinetics.Noncomp);
calcvalue=num2str(simVel(enzobj,0.5,0.4));
expvalue=num2str(0.15);
verifyEqual(testCase,calcvalue,expvalue);
end


function testfunc7(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'kfi',5,'kri',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.Uncomp,MMenKinetics.Elem,rt);
mathml = enzobj.toMathExpr;
calcvalue=mathml;
mathexp.sbmlmath = '(Vm/(1+Ci/Ki))*S/(Km/(1+Ci/Ki)+S)';
mathexp.speciesnames = {'S','Ci'};
mathexp.rtconsts = {'Vm','Km','Ki'};
expvalue=mathexp;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc8(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'kfi',5,'kri',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.Uncomp,MMenKinetics.Elem,rt);
ele2ems(enzobj,MMenKinetics.Uncomp);
calcvalue=enzobj.sumconsts;
enzobj.sumconsts.vm=0.6;
enzobj.sumconsts.km=0.5;
enzobj.sumconsts.ki=0.4;
expvalue=enzobj.sumconsts;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc9(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'kfi',5,'kri',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.Uncomp,MMenKinetics.Elem,rt);
ele2ems(enzobj,MMenKinetics.Uncomp);
calcvalue=num2str(simVel(enzobj,0.5,0.4));
expvalue=num2str(0.2);
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc10(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('vm',3,'km',6,'ki',2);
enzobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt);
mathml = enzobj.toMathExpr;
calcvalue=mathml;
mathexp.sbmlmath = 'Vm*S/(Km*(1+Ci/Ki)+S)';
mathexp.speciesnames = {'S','Ci'};
mathexp.rtconsts = {'Vm','Km','Ki'};
expvalue=mathexp;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc11(testCase)
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('vm',3,'km',6,'ki',2);
ad     = struct('ki2',3,'ki3',1,'ki4',2);
enzobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt,ad);
mathml = enzobj.toMathExpr;
calcvalue=mathml;
mathexp.mathexpr = 'Vm*S/(S+Km(1+Ci1/Ki1+Ci2/Ki2+Ci3/Ki3+Ci4/Ki4))';
mathexp.speciesnames = {'S'  'Ci1'  'Ci2'  'Ci3'  'Ci4'};
mathexp.rtconsts = {'Vm'  'Km'  'Ki1'  'Ki2'  'Ki3'  'Ki4'};
expvalue=mathexp;
verifyEqual(testCase,calcvalue,expvalue);
end















