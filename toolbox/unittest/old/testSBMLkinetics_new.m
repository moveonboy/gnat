% test example 
clc;clear

theCompt  = Compt('testsbml');

%% create an array of species
m5species   = GlycanSpecies(glycanMLread('M5.glycoct_xml'),theCompt);
m5gnspecies = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),theCompt);
allglycans  = CellArrayList;
allglycans.add(m5species);
allglycans.add(m5gnspecies);

%% create an array of reactions
allrxns  = CellArrayList;

% load enzyme database
enzdbmatfilename   = 'glyenzDB.mat';
enzdb  = enzdbmatLoad(enzdbmatfilename);
mgat1  = enzdb('mgat1');

rt     = struct('kf',10,'kr',5,'kcat',2,'enzconc',0.3);
enzkineticsobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Elem,rt);
mgat1.enzkinetics = enzkineticsobj;

%set up rxn
m5gnrxn = Rxn(m5species,m5gnspecies,mgat1);
% m5gnrxn.setrxnkinetics;
allrxns.add(m5gnrxn);

%set up enzyme pool
allenzs = CellArrayList;
allenzs.add(mgat1);

testpathway = Pathway;
testpathway.addGlycansByStruct(allglycans);
testpathway.addRxnsByStruct(allrxns)
testpathway.addEnzs(allenzs);

% set ID for species, rxn and kinetics for SBML output
testpathway.setEnzID;
testpathway.setSpeciesID;
testpathway.setRxnsID;
testpathway.setRxnKinetics;

modelname ='testSBML';
testmodel = GlycanNetModel(theCompt,testpathway,modelname);
testsbmlmodel = testmodel.toSBMLStruct;
glycanNetViewer(testmodel);
OutputSBML(testsbmlmodel);

