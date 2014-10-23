residueMap=load('residueTypes.mat');

%mgat1
mgat1                    = GTEnz([2;4;1;101]);
mgat1.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('3','1');
mgat1.resAtt2FG          = manResType;
mgat1. linkresAtt2FG     = struct('bond',manBond,'anomer','a');
glcnacbond               = GlycanBond('2','1');
mgat1.linkFG             = struct('anomer','b','bond',glcnacbond);
a2                       = glycanMLread('1579.78.glycoct_xml');
mgat1.substMinStruct     = a2;
mgat1.substMaxStruct     = a2;
mgat1.targetBranch       = glycanMLread('mgat1targetbranch.glycoct_xml');

%mgat2
mgat2                    = GTEnz([2;4;1;143]);
mgat2.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('6','1');
mgat2.resAtt2FG          = manResType;
mgat2. linkresAtt2FG     = struct('bond',manBond,'anomer','a');
glcnacbond               = GlycanBond('2','1');
mgat2.linkFG             = struct('anomer','b','bond',glcnacbond);
gn2m3gn                  = glycanMLread('gnm3gn.glycoct_xml');
mgat2.substMinStruct     = gn2m3gn;
mgat2.substNABranch      = glycanMLread('mgat2NAbranch2.glycoct_xml');
mgat2.substNAResidue     = residueMap.allresidues('Gal');

%mgat3
mgat3                    = GTEnz([2;4;1;143]);
mgat3.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('4','1');
mgat3.resAtt2FG          = manResType;
mgat3.linkresAtt2FG      = struct('bond', manBond,'anomer','b');
glcnacbond               = GlycanBond('4','1');
mgat3.linkFG             = struct('anomer','b','bond',glcnacbond);
m3gn                     = glycanMLread('m3gn.glycoct_xml');
mgat3.substMinStruct     = m3gn;
mgat3.substNABranch      = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
mgat3.targetBranch       = glycanMLread('mgat3targetbranch.glycoct_xml');
mgat3.substNAResidue     = residueMap.allresidues('Gal');

%mgat4
mgat4                    = GTEnz([2;4;1;145]);
mgat4.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('3','1');
mgat4.resAtt2FG          = manResType;
mgat4.linkresAtt2FG      = struct('bond', manBond,'anomer','a');
glcnacbond               = GlycanBond('4','1');
mgat4.linkFG             = struct('anomer','b','bond',glcnacbond);
gn2m3gn                  = glycanMLread('gnm3gn.glycoct_xml');
mgat4.substMinStruct     = gn2m3gn;
mgat4.substNABranch      = CellArrayList;
mgat4.substNABranch.add(glycanMLread('mgat4NAbranch1.glycoct_xml'));
mgat4.substNABranch.add(glycanMLread('mgat4NAbranch2.glycoct_xml'));

%mgat5
mgat5                    = GTEnz([2;4;1;155]);
mgat5.resfuncgroup       = residueMap.allresidues('GlcNAc');
manResType               = residueMap.allresidues('Man');
manBond                  = GlycanBond('6','1');
mgat5.resAtt2FG          = manResType;
mgat5.linkresAtt2FG      = struct('bond', manBond,'anomer','a');
glcnacbond               = GlycanBond('6','1');
mgat5.linkFG             = struct('anomer','b','bond',glcnacbond);
gn2m3gn                  = glycanMLread('mgat5minimal.glycoct_xml');
mgat5.substMinStruct     = gn2m3gn;
mgat5.substNABranch      =CellArrayList;
mgat5.substNABranch.add(glycanMLread('mgat5NAbranch1.glycoct_xml'));
mgat5.substNABranch.add(glycanMLread('mgat5NAbranch2.glycoct_xml'));

%manii
manResType               = residueMap.allresidues('Man');
manii                    = GHEnz([3;2;1;114]);
manii.resfuncgroup       = manResType;
manBond(1,1)             = GlycanBond('3','1');
manBond(2,1)             = GlycanBond('6','1');
manii.linkFG             = struct('bond',manBond ,'anomer','a');
manii.resAtt2FG          = manResType;
resbond                  = GlycanBond('6','1');
manii.linkresAtt2FG      = struct('bond', resbond,'anomer','a');
manii.substNABranch      = glycanMLread('1008.51.glycoct_xml');
manii.substMaxStruct     = glycanMLread('1824.91.glycoct_xml');

%mani
mani                     = GTEnz([3;2;1;113]);
mani.resfuncgroup        = residueMap.allresidues('Man');
mani.resAtt2FG           = residueMap.allresidues('Man');
manBond                  = GlycanBond('2','1');
mani.linkFG              = struct('anomer','a','bond',manBond);
manunknownbond           = GlycanBond('?','1');
mani.linkresAtt2FG       = struct('bond', manunknownbond,'anomer','a');
mani.substMinStruct      = glycanMLread('maniminstruct.glycoct_xml');
mani.substMaxStruct      = glycanMLread('manimaxstruct.glycoct_xml');

%ST6GalI
ST6GalI                  = GTEnz([2;4;99;6]);
ST6GalI.isTerminalTarget = true;
ST6GalI.resfuncgroup     = residueMap.allresidues('NeuAc');
siaTbond                 = GlycanBond('3','2');
ST6GalI.linkFG           = struct('anomer','a','bond',siaTbond);
galResType               = residueMap.allresidues('Gal');
galBond                  = GlycanBond('4','1');
ST6GalI.resAtt2FG        = galResType;
ST6GalI.linkresAtt2FG    = struct('bond', galBond,'anomer','b');
ST6GalI.substNABranch    = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
ST6GalI.targetbranchcontain = glycanMLread('gntetargetbranch.glycoct_xml');

%Fut8
Fut8                     = GTEnz([2;4;1;68]);
Fut8.resfuncgroup        = residueMap.allresidues('Fuc');
Fut8.resAtt2FG           = residueMap.allresidues('GlcNAc');
fucbond                  = GlycanBond('6','1');
Fut8.linkFG              = struct('anomer','a','bond',fucbond);
gnbond                   = GlycanBond('?','?');
Fut8.linkresAtt2FG       = struct('bond',gnbond,'anomer','b');
gnm3                     = glycanMLread('fucminstruct.glycoct_xml');
Fut8.substMinStruct      = gnm3;
Fut8.targetBranch        = glycanMLread('fuctargetbranch.glycoct_xml');

%b4GalI
b4GalI                   = GTEnz([2;4;1;38]);
b4GalI.resfuncgroup      = residueMap.allresidues('Gal');
b4GalI.resAtt2FG         = residueMap.allresidues('GlcNAc');
fucbond                  = GlycanBond('4','1');
b4GalI.linkFG            = struct('anomer','b','bond',fucbond);
gnbond                   = GlycanBond('?','1');
b4GalI.linkresAtt2FG     = struct('bond',gnbond,'anomer','b');
gn2m3gn                  = glycanMLread('b4galiminstruct.glycoct_xml');
b4GalI.substMinStruct    = gn2m3gn;
b4GalI.targetNABranch    = glycanMLread('b4galitargetNAbranch.glycoct_xml');

%iGnt
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

enzArray=CellArrayList;
enzArray.add(mgat1);
enzArray.add(mgat2);
enzArray.add(mgat3);
enzArray.add(mgat4);
enzArray.add(mgat5);
enzArray.add(manii);
enzArray.add(mani);
enzArray.add(iGnt);
enzArray.add(b4GalI);
enzArray.add(Fut8);
enzArray.add(ST6GalI);

%define the end prod
a1  = GlycanSpecies(glycanMLread('M5.glycoct_xml'));
% a2  = GlycanSpecies(glycanMLread('3142.58bi.glycoct_xml'));
a3  = GlycanSpecies(glycanMLread('3864.92.glycoct_xml'));
% a4  = GlycanSpecies(glycanMLread('3142.58tri.glycoct_xml'));

glycanArray = CellArrayList;
glycanArray.add(a1);
% glycanArray.add(a2);
glycanArray.add(a3);
% glycanArray.add(a4);

%Perform reaction
displayOptions = displayset('showMass',true,'showLinkage',true,'showRedEnd',true);
fprintf(1,'Input of glycan product structure is \n');
glycanViewer(a1.glycanStruct,displayOptions);
% glycanViewer(a2.glycanStruct,displayOptions);
glycanViewer(a3.glycanStruct,displayOptions);
% glycanViewer(a4.glycanStruct,displayOptions);
[isPath,nlinkedpath]=inferGlyConnPath(glycanArray,enzArray);
fprintf(1,'Inferred network is shown below:\n');
glycanPathViewer(nlinkedpath);