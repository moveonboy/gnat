residueMap=load('residueTypes.mat');

iGnt                     = GTEnz([2;4;1;149]);
iGnt.resfuncgroup        = residueMap.allresidues('GlcNAc');
galResType               = residueMap.allresidues('Gal');
iGnt.resAtt2FG           = galResType;
manBond                  = GlycanBond('4','1');
iGnt. linkresAtt2FG      = struct('bond',manBond,'anomer','b');
glcnacbond               = GlycanBond('3','1');
iGnt.linkFG              = struct('anomer','b','bond',glcnacbond);
iGnt.substMinStruct      = glycanMLread('igntminstruct.glycoct_xml');
iGnt.targetbranchcontain = glycanMLread('ignttargetbranchcontain.glycoct_xml');

a1  = GlycanSpecies(glycanMLread('2489.25bi.glycoct_xml'));

[numSubstr,substrSpecies] = inferGlySubstr(a1,iGnt);