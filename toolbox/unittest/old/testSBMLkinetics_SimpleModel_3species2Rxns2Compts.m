% test example 
clc;clear

theCompt                  = CellArrayList;
cisCompt                  = Compt('cis-golgi');
mediaCompt                = Compt('media-golgi');
cisCompt.posteriorcompt   = mediaCompt;
mediaCompt.priorcompt     = cisCompt;
theCompt.add(cisCompt);
theCompt.add(mediaCompt);

%% create an array of species
m5gnspecies_media    = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),mediaCompt);
m4gnspecies_media    = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),mediaCompt);
m3gnspecies_media    = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),mediaCompt);

m5gnspecies_cis      = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),cisCompt);
m4gnspecies_cis      = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),cisCompt);
m3gnspecies_cis      = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),cisCompt);

allglycans  = CellArrayList;
allglycans.add(m5gnspecies_cis);
allglycans.add(m4gnspecies_cis);
allglycans.add(m3gnspecies_cis);
allglycans.add(m5gnspecies_media);
allglycans.add(m4gnspecies_media);
allglycans.add(m3gnspecies_media);

%% create an array of reactions
allrxns  = CellArrayList;

% load enzyme database
enzdbmatfilename   = 'glyenzDB.mat';
enzdb              = enzdbmatLoad(enzdbmatfilename);
manii              = enzdb('manii');

rt8                = struct('Vm',300,'Km',100);
enzkineticsobj     = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt8);
specistruct        = CellArrayList;
specistruct.add(m5gnspecies_media.glycanStruct);
specifiratio       = CellArrayList;
ratio1         = struct('Vm',1,'Km',2);
specifiratio.add(ratio1);
enzkineticsobj.setSubstrateSpecificityK(specistruct,specifiratio,1);
manii.enzkinetics = enzkineticsobj;

%set up rxn
m4gnrxn_media = Rxn(m5gnspecies_media,m4gnspecies_media,manii);
m4gnrxn_cis = Rxn(m5gnspecies_cis,m4gnspecies_cis,manii);
%m4gnrxn.setrxnkinetics;
allrxns.add(m4gnrxn_media);
allrxns.add(m4gnrxn_cis);

m3gnrxn_cis = Rxn(m4gnspecies_cis,m3gnspecies_cis,manii); 
m3gnrxn_media = Rxn(m4gnspecies_media,m3gnspecies_media,manii); 

%m3gnrxn.setrxnkinetics;
allrxns.add(m3gnrxn_cis);
allrxns.add(m3gnrxn_media);

% set up enzyme pool
allenzs = CellArrayList;
allenzs.add(manii);

testpathway = Pathway;
testpathway.addGlycans(allglycans);
testpathway.addRxns(allrxns)
testpathway.addEnzs(allenzs);

testpathway.setSpeciesID;
testpathway.setRxnsID;
testpathway.setEnzID;
testpathway.setEnzRxnsKinetics;

testpathway.setComptReactorType(Reactor.CSTR);  
% set qflow and k in each compts
qflow  = 10;
k      = 2;
testpathway.setComptsQflow(qflow);
testpathway.setComptskt(k);

%first compt
speciesInitConc = CellArrayList;
speciesInitConc.add(struct('species',m5gnspecies_cis,'initconc',0.8));
speciesInitConc.add(struct('species',m4gnspecies_cis,'initconc',0.1));
speciesInitConc.add(struct('species',m3gnspecies_cis,'initconc',0.1));
cisCompt.reactortype.setSpeciesInitDB(speciesInitConc)

testpathway.setTransporRxnKinetics;

modelname ='testSBML';
testmodel = GlycanNetModel(theCompt,testpathway,modelname);
testsbmlstruct= testmodel.toSBMLStruct;
glycanNetViewer(testmodel);
OutputSBML(testsbmlstruct);

