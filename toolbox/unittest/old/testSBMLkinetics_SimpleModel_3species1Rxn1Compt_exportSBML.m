% test example 
clc;clear

thecisCompt  = Compt('cis-golgi');
theCompt = CellArrayList;
theCompt.add(thecisCompt);

%% create an array of species
m5gnspecies    = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),thecisCompt);
m4gnspecies    = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),thecisCompt);
m3gnspecies    = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),thecisCompt);

allglycans  = CellArrayList;
allglycans.add(m5gnspecies);
allglycans.add(m4gnspecies);
allglycans.add(m3gnspecies);

%% create an array of reactions
allrxns  = CellArrayList;

% load enzyme database
enzdbmatfilename   = 'glyenzDB.mat';
enzdb  = enzdbmatLoad(enzdbmatfilename);
manii  = enzdb('manii');

rt8            = struct('Vm',300,'Km',100);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt8);
specistruct    = CellArrayList;
specistruct.add(m5gnspecies.glycanStruct);
specifiratio   = CellArrayList;
ratio1         = struct('Vm',1,'Km',2);
specifiratio.add(ratio1);
enzkineticsobj.setSubstrateSpecificityK(specistruct,specifiratio,1);
manii.enzkinetics = enzkineticsobj;

m4gnrxn = Rxn(m5gnspecies,m4gnspecies,manii);
allrxns.add(m4gnrxn);

m3gnrxn = Rxn(m4gnspecies,m3gnspecies,manii); 
allrxns.add(m3gnrxn);

%set up enzyme pool
allenzs = CellArrayList;
allenzs.add(manii);

testpathway = Pathway;
testpathway.addGlycans(allglycans);
testpathway.addRxns(allrxns)
testpathway.addEnzs(allenzs);
testpathway.setSpeciesID;
testpathway.setRxnsID;
testpathway.setEnzID;
testpathway.setRxnsKinetics;

testpathway.setComptReactorType(Reactor.CSTR);
qflow  = 10;
k      = 2;
testpathway.setComptsQflow(qflow);
testpathway.setComptskt(k);

modelname ='testSBML';
testmodel = GlycanNetModel(theCompt,testpathway,modelname);
testsbmlstruct= testmodel.toSBMLStruct;
glycanNetViewer(testmodel);
% OutputSBML(testsbmlstruct);

