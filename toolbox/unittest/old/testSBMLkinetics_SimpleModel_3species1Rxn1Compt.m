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

% rt1     = struct('Vm',450,'Km',260);
% enzkineticsobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt1);
% mgat1.enzkinetics = enzkineticsobj;
% 
% rt2     = struct('Vm',140,'Km',190);
% enzkineticsobj= MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt2);
% mgat2.enzkinetics = enzkineticsobj;
% 
% rt3     = struct('Vm',4000,'Km',190);
% enzkineticsobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt3);
% specistruct= CellArrayList;
% specistruct.add(m5gnspecies.glycanStruct);
% specistruct.add(m4gnspecies.glycanStruct);
% specistruct.add(m3gnspecies.glycanStruct);
% specifiratio = CellArrayList;
% ratio1 = struct('Vm',1,'Km',21.053);
% ratio2 = struct('Vm',1,'Km',21.053);
% ratio3 = struct('Vm',1,'Km',21.053);
% specifiratio.add(ratio1);
% specifiratio.add(ratio2);
% specifiratio.add(ratio3);
% setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
% mgat3.enzkinetics = enzkineticsobj;
% 
% rt4     = struct('Vm',10,'Km',3400);
% enzkineticsobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt4);
% mgat4.enzkinetics = enzkineticsobj;
% 
% rt5     = struct('Vm',10,'Km',130);
% enzkineticsobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt5);
% specistruct= CellArrayList;
% specistruct.add(m3gn3species.glycanStruct);
% specifiratio = CellArrayList;
% ratio1 = struct('Vm',1,'Km',0.692);
% specifiratio.add(ratio1);
% setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
% mgat5.enzkinetics = enzkineticsobj;
% 
% rt6     = struct('Vm',300,'Km',100);
% enzkineticsobj6 = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt6);
% mani.enzkinetics = enzkineticsobj6;
% 
% rt7     = struct('Vm',580,'Km',130);
% enzkineticsobj = MMenKinetics(MMenKinetics.SMM,MMenKinetics.Ensm,rt7);
% specistruct= CellArrayList;
% specistruct.add(m5gnspecies.glycanStruct);
% specistruct.add(m4gnspecies.glycanStruct);
% specistruct.add(m3gnspecies.glycanStruct);
% specistruct.add(m3gn3bspecies.glycanStruct);
% specistruct.add(m3gn3species.glycanStruct);
% specistruct.add(m3gn4species.glycanStruct);
% specistruct.add(m5gngnbspecies.glycanStruct);
% specistruct.add(m4gngnbspecies.glycanStruct);
% specistruct.add(m3gngnbspecies.glycanStruct);
% specistruct.add(m3gn2gnbspecies.glycanStruct);
% specistruct.add(m3gn3bgnbspecies.glycanStruct);
% specistruct.add(m3gn3gnbspecies.glycanStruct);
% specistruct.add(m3gn4gnbspecies.glycanStruct);
% specifiratio = CellArrayList;
% ratio1 = struct('Vm',1,'Km',30.769);
% ratio2 = struct('Vm',1,'Km',30.769);
% ratio3 = struct('Vm',1,'Km',30.769);
% ratio4 = struct('Vm',1,'Km',0.538);
% ratio5 = struct('Vm',1,'Km',0.385);
% ratio6 = struct('Vm',1,'Km',0.308);
% ratio7 = struct('Vm',1,'Km',30.769);
% ratio8 = struct('Vm',1,'Km',30.769);
% ratio9 = struct('Vm',1,'Km',30.769);
% ratio10 = struct('Vm',1,'Km',3.846);
% ratio11 = struct('Vm',1,'Km',1.692);
% ratio12 = struct('Vm',1,'Km',1.538);
% ratio13 = struct('Vm',1,'Km',1.077);
% specifiratio.add(ratio1);
% specifiratio.add(ratio2);
% specifiratio.add(ratio3);
% specifiratio.add(ratio4);
% specifiratio.add(ratio5);
% specifiratio.add(ratio6);
% specifiratio.add(ratio7);
% specifiratio.add(ratio8);
% specifiratio.add(ratio9);
% specifiratio.add(ratio10);
% specifiratio.add(ratio11);
% specifiratio.add(ratio12);
% specifiratio.add(ratio13);
% setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
% galt.enzkinetics = enzkineticsobj;

rt8            = struct('Vm',300,'Km',100);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt8);
specistruct    = CellArrayList;
specistruct.add(m5gnspecies.glycanStruct);
specifiratio   = CellArrayList;
ratio1         = struct('Vm',1,'Km',2);
specifiratio.add(ratio1);
enzkineticsobj.setSubstrateSpecificityK(specistruct,specifiratio,1);
manii.enzkinetics = enzkineticsobj;

%set up rxn
% m8rxn   = Rxn(m9species,m8species,mani);
% %m8rxn.setrxnkinetics;
% allrxns.add(m8rxn);
% 
% m7rxn   = Rxn(m8species,m7species,mani);
% %m7rxn.setrxnkinetics;
% allrxns.add(m7rxn);
% 
% m6rxn   = Rxn(m7species,m6species,mani);
% %m6rxn.setrxnkinetics;
% allrxns.add(m6rxn);
% 
% m5rxn   = Rxn(m6species,m5species,mani);
% %m5rxn.setrxnkinetics;
% allrxns.add(m5rxn);
% 
% m5gnrxn = Rxn(m5species,m5gnspecies,mgat1);
%m5gnrxn.setrxnkinetics;
% allrxns.add(m5gnrxn);

m4gnrxn = Rxn(m5gnspecies,m4gnspecies,manii);
%m4gnrxn.setrxnkinetics;
allrxns.add(m4gnrxn);

m3gnrxn = Rxn(m4gnspecies,m3gnspecies,manii); 
%m3gnrxn.setrxnkinetics;
allrxns.add(m3gnrxn);

% m3gn2rxn = Rxn(m3gnspecies,m3gn2species,mgat2);
% %m3gn2rxn.setrxnkinetics;
% allrxns.add(m3gn2rxn);
% 
% m3gn3brxn = Rxn(m3gn2species,m3gn3bspecies,mgat5);
% %m3gn3brxn.setrxnkinetics;
% allrxns.add(m3gn3brxn);
% 
% m3gn3rxn  = Rxn(m3gn2species,m3gn3species,mgat4);
% %m3gn3rxn.setrxnkinetics;
% allrxns.add(m3gn3rxn);
% 
% m3gn4arxn = Rxn(m3gn3bspecies,m3gn4species,mgat4);
% %m3gn4arxn.setrxnkinetics;
% allrxns.add(m3gn4arxn);
% 
% m3gn4brxn = Rxn(m3gn3species,m3gn4species,mgat5);
% %m3gn4brxn.setrxnkinetics;
% allrxns.add(m3gn4brxn);
% 
% m5gngrxn = Rxn(m5gnspecies,m5gngspecies,galt);
% %m5gngrxn.setrxnkinetics;
% allrxns.add(m5gngrxn);
% 
% m4gngrxn = Rxn(m4gnspecies,m4gngspecies,galt);
% %m4gngrxn.setrxnkinetics;
% allrxns.add(m4gngrxn);
% 
% m3gngrxn = Rxn(m3gnspecies,m3gngspecies,galt);
% %m3gngrxn.setrxnkinetics;
% allrxns.add(m3gngrxn);
% 
% m3gn2grxn = Rxn(m3gn2species,m3gn2gspecies,galt);
% %m3gn2grxn.setrxnkinetics;
% allrxns.add(m3gn2grxn);
% 
% m3gn3bgrxn = Rxn(m3gn3bspecies,m3gn3bgspecies,galt);
% %m3gn3bgrxn.setrxnkinetics;
% allrxns.add(m3gn3bgrxn);
% 
% m3gn3grxn = Rxn(m3gn3species,m3gn3gspecies,galt);
% %m3gn3grxn.setrxnkinetics;
% allrxns.add(m3gn3grxn);
% 
% m3gn4grxn = Rxn(m3gn4species,m3gn4gspecies,galt);
% %m3gn4grxn.setrxnkinetics;
% allrxns.add(m3gn4grxn);
% 
% m5gngnbrxn = Rxn(m5gnspecies,m5gngnbspecies,mgat3);
% %m5gngnbrxn.setrxnkinetics;
% allrxns.add(m5gngnbrxn);
% 
% m4gngnbrxn = Rxn(m4gnspecies,m4gngnbspecies,mgat3);
% %m4gngnbrxn.setrxnkinetics;
% allrxns.add(m4gngnbrxn);
% 
% m3gngnbrxn = Rxn(m3gnspecies,m3gngnbspecies,mgat3);
% %m3gngnbrxn.setrxnkinetics;
% allrxns.add(m3gngnbrxn);
% 
% m3gn2gnbrxn = Rxn(m3gn2species,m3gn2gnbspecies,mgat3);
% %m3gn2gnbrxn.setrxnkinetics;
% allrxns.add(m3gn2gnbrxn);
% 
% m3gn3bgnbrxn = Rxn(m3gn3bspecies,m3gn3bgnbspecies,mgat3);
% %m3gn3bgnbrxn.setrxnkinetics;
% allrxns.add(m3gn3bgnbrxn);
% 
% m3gn3gnbrxn = Rxn(m3gn3species,m3gn3gnbspecies,mgat3);
% %m3gn3gnbrxn.setrxnkinetics;
% allrxns.add(m3gn3gnbrxn);
% 
% m3gn4gnbrxn = Rxn(m3gn4species,m3gn4gnbspecies,mgat3);
% %m3gn4gnbrxn.setrxnkinetics;
% allrxns.add(m3gn4gnbrxn);
% 
% m5gngnbgrxn = Rxn(m5gngnbspecies,m5gngnbgspecies,galt);
% %m5gngnbgrxn.setrxnkinetics;
% allrxns.add(m5gngnbgrxn);
% 
% m4gngnbgrxn = Rxn(m4gngnbspecies,m4gngnbgspecies,galt);
% %m4gngnbgrxn.setrxnkinetics;
% allrxns.add(m4gngnbgrxn);
% 
% m3gngnbgrxn = Rxn(m3gngnbspecies,m3gngnbgspecies,galt);
% %m3gngnbgrxn.setrxnkinetics;
% allrxns.add(m3gngnbgrxn);
% 
% m3gn2gnbgrxn = Rxn(m3gn2gnbspecies,m3gn2gnbgspecies,galt);
% %m3gn2gnbgrxn.setrxnkinetics;
% allrxns.add(m3gn2gnbgrxn);
% 
% m3gn3bgnbgrxn = Rxn(m3gn3bgnbspecies,m3gn3bgnbgspecies,galt);
% %m3gn3bgnbgrxn.setrxnkinetics;
% allrxns.add(m3gn3bgnbgrxn);
% 
% m3gn3gnbgrxn = Rxn(m3gn3gnbspecies,m3gn3gnbgspecies,galt);
% %m3gn3gnbgrxn.setrxnkinetics;
% allrxns.add(m3gn3gnbgrxn);
% 
% m3gn4gnbgrxn = Rxn(m3gn4gnbspecies,m3gn4gnbgspecies,galt);
% %m3gn4gnbgrxn.setrxnkinetics;
% allrxns.add(m3gn4gnbgrxn);

%set up enzyme pool
allenzs = CellArrayList;
% allenzs.add(mgat1);
% allenzs.add(mgat2);
% allenzs.add(mgat3);
% allenzs.add(mgat4);
% allenzs.add(mgat5);
% allenzs.add(mani);
allenzs.add(manii);
% allenzs.add(galt);

testpathway = Pathway;
testpathway.addGlycans(allglycans);
testpathway.addRxns(allrxns)
testpathway.addEnzs(allenzs);
testpathway.setSpeciesID;
testpathway.setRxnsID;
testpathway.setEnzID;
testpathway.setRxnsKinetics;

% testpathway.setComptReactorType(Reactor.BatchR);
testpathway.setComptReactorType(Reactor.CSTR);

% set compartment parameter

modelname ='testSBML';
testmodel = GlycanNetModel(theCompt,testpathway,modelname);
% testsbmlstruct= testmodel.toSBMLStruct;
glycanNetViewer(testmodel);
% OutputSBML(testsbmlstruct);

