function tests = EleRxnKineticstest
    tests = functiontests(localfunctions);
end

function testfunc1(testCase)
speciesR1 = CellArrayList;
speciesP1 = CellArrayList;
speciesR1.add(struct('speciesid','A1','coef',2))
speciesR1.add(struct('speciesid','A2','coef',0.3))
speciesR1.add(struct('speciesid','A3','coef',0.2))
speciesP1.add(struct('speciesid','B1','coef',0.2))
speciesP1.add(struct('speciesid','B2','coef',0.4))
rateconsts = struct('kf',5,'kr',2);
obj1       = EleRxnKinetics(speciesR1,speciesP1,rateconsts);
calcvalue  = obj1.kineticsconst;
expvalue = rateconsts;
verifyEqual(testCase,calcvalue,expvalue);
end

function testfunc2(testCase)
speciesR2 = CellArrayList;
speciesP2 = CellArrayList;
speciesR2.add(struct('speciesid','A1','coef',2))
speciesR2.add(struct('speciesid','A2','coef',0.3))
speciesR2.add(struct('speciesid','A3','coef',0.2))
speciesP2.add(struct('speciesid','B1','coef',0.2))
speciesP2.add(struct('speciesid','B2','coef',0.4))
rateconsts = struct('kf',5,'kr',2);
obj2       = EleRxnKinetics(speciesR2,speciesP2,rateconsts);
setMathFormula(obj2,'rxn1')
calcvalue  = obj2.kineticlaws.mathformula;
expvalue = 'kf_rxn1*(A1^a1_rxn1)*(A2^a2_rxn1)*(A3^a3_rxn1)-kr_rxn1*(B1^b1_rxn1)*(B2^b2_rxn1)';
verifyEqual(testCase,calcvalue,expvalue);
end