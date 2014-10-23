function [ enzobj, mathml] = testMMkinetics()
%TESTMMKINETICS Summary of this function goes here
%   Detailed explanation goes here
rt     = struct('kf',10,'kr',5,'kcat',2,'enzconc',0.3);
enzobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Elem,rt);
mathml = enzobj.toMathExpr;
end

