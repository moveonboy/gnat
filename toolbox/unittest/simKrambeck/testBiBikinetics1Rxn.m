% testBiBiKinetics for a single rxn
function testBiBikinetics1Rxn()
cisCompt          =  Compt('cis_golgi');
enzdbmatfilename  =  'glyenzDB.mat';
enzdb             =  enzdbmatLoad(enzdbmatfilename);
mgat1             =  enzdb('mgat1');

listofstructs     =  setGlyStruct;
m5species         =  GlycanSpecies(listofstructs('M5'),cisCompt);
m5gnspecies       =  GlycanSpecies(listofstructs('M5gn'),cisCompt);
m5species.id      = 'm5';
m5gnspecies.id    = 'm5gn';
m5gnrxn           =  Rxn(m5species,m5gnspecies,mgat1);

rt1               =  struct('Vm',1,'Kmd_G1E',2,'Kmd_SNE',3,'Keq',4);
enzkineticsobj    =  BiBiKinetics(mgat1,BiBiKinetics.bibiSeqOrderwSubstInhib,rt1);
mgat1.enzkinetics =  enzkineticsobj;

allrxns = CellArrayList;
allrxns.add(m5gnrxn);

allenzs = CellArrayList;
allenzs.add(mgat1);

allspecies = CellArrayList;
allspecies.add(m5species);
allspecies.add(m5gnspecies);

testpathway = Pathway;
testpathway.addGlycans(allspecies);
testpathway.addRxns(allrxns)
testpathway.addEnzs(allenzs);
testpathway.setRxnsID;
testpathway.setEnzID;

sugarnucldb = setsugarstruct();
testpathway.setEnzRxnsKinetics(sugarnucldb);  % BiBiRxnKinetics need a database for sugar nucleotide database
theCompt = CellArrayList;
theCompt.add(cisCompt);

modelname       = 'testBiBi';
testmodel       = GlycanNetModel(theCompt,testpathway,modelname);
testsbmlstruct  = testmodel.toSBMLStruct;
OutputSBML(testsbmlstruct);
disp('end');
end

function sugarnucldb = setsugarstruct()
sugarnucldb = containers.Map;
sugarnucldb('UDP_GlcNAc')=1;
sugarnucldb('UDP')=2;
end

function listofstructs = setGlyStruct()
listofstructs = containers.Map;
listofstructs('M9')=glycanMLread('M9.glycoct_xml');
listofstructs('M8')=glycanMLread('M8.glycoct_xml');
listofstructs('M7')=glycanMLread('M7.glycoct_xml');
listofstructs('M6')=glycanMLread('M6.glycoct_xml');
listofstructs('M5')=glycanMLread('M5.glycoct_xml');
listofstructs('M5gn')=glycanMLread('M5gn.glycoct_xml');
listofstructs('M4gn')=glycanMLread('M4gn.glycoct_xml');
listofstructs('M3gn')=glycanMLread('M3gn.glycoct_xml');
listofstructs('M3gn2')=glycanMLread('M3gn2.glycoct_xml');
listofstructs('M3gn3')=glycanMLread('M3gn3.glycoct_xml');
listofstructs('M3gn3b')=glycanMLread('M3gn3b.glycoct_xml');
listofstructs('M3gn4')=glycanMLread('M3gn4.glycoct_xml');
listofstructs('M5gng')=glycanMLread('M5gng.glycoct_xml');
listofstructs('M4gng')=glycanMLread('M4gng.glycoct_xml');
listofstructs('M3gng')=glycanMLread('M3gng.glycoct_xml');
listofstructs('M3gn2g')=glycanMLread('M3gn2g.glycoct_xml');
listofstructs('M3gn3g')=glycanMLread('M3gn3g.glycoct_xml');
listofstructs('M3gn3bg')=glycanMLread('M3gn3bg.glycoct_xml');
listofstructs('M3gn4g')=glycanMLread('M3gn4g.glycoct_xml');
listofstructs('M5gngnb')=glycanMLread('M5gngnb.glycoct_xml');
listofstructs('M4gngnb')=glycanMLread('M4gngnb.glycoct_xml');
listofstructs('M3gngnb')=glycanMLread('M3gngnb.glycoct_xml');
listofstructs('M3gn2gnb')=glycanMLread('M3gn2gnb.glycoct_xml');
listofstructs('M3gn3gnb')=glycanMLread('M3gn3gnb.glycoct_xml');
listofstructs('M3gn3bgnb')=glycanMLread('M3gn3bgnb.glycoct_xml');
listofstructs('M3gn4gnb')=glycanMLread('M3gn4gnb.glycoct_xml');
listofstructs('M5gngnbg')=glycanMLread('M5gngnbg.glycoct_xml');
listofstructs('M4gngnbg')=glycanMLread('M4gngnbg.glycoct_xml');
listofstructs('M3gngnbg')=glycanMLread('M3gngnbg.glycoct_xml');
listofstructs('M3gn2gnbg')=glycanMLread('M3gn2gnbg.glycoct_xml');
listofstructs('M3gn3gnbg')=glycanMLread('M3gn3gnbg.glycoct_xml');
listofstructs('M3gn3bgnbg')=glycanMLread('M3gn3bgnbg.glycoct_xml');
listofstructs('M3gn4gnbg')=glycanMLread('M3gn4gnbg.glycoct_xml');
end
