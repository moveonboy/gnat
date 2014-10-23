residueMap=load('residueTypes.mat');

%mgat1
mgat1=GTEnz([2;4;1;101]);
mgat1.resfuncgroup=residueMap.allresidues('GlcNAc');
manResType=residueMap.allresidues('Man');
manBond=GlycanBond('3','1');
mgat1.resAtt2FG=manResType;
mgat1. linkresAtt2FG=struct('bond',manBond,'anomer','a');
glcnacbond=GlycanBond('2','1');
mgat1.linkFG=struct('anomer','b','bond',glcnacbond);
a2=glycanMLread('m5gn.glycoct_xml');
mgat1.substMinStruct=a2;
mgat1.substMaxStruct=a2;
mgat1.targetBranch=glycanMLread('mgat1targetbranch.glycoct_xml');

%mgat2
mgat2=GTEnz([2;4;1;143]);
mgat2.resfuncgroup=residueMap.allresidues('GlcNAc');
manResType=residueMap.allresidues('Man');
manBond=GlycanBond('6','1');
mgat2.resAtt2FG=manResType;
mgat2. linkresAtt2FG=struct('bond',manBond,'anomer','a');
glcnacbond=GlycanBond('2','1');
mgat2.linkFG=struct('anomer','b','bond',glcnacbond);
gn2m3gn=glycanMLread('gnm3gn.glycoct_xml');
mgat2.substMinStruct=gn2m3gn;
mgat2.substNABranch=CellArrayList;
mgat2.substNABranch.add(glycanMLread('mgat2NAbranch1.glycoct_xml'));
mgat2.substNABranch.add(glycanMLread('mgat2NAbranch2.glycoct_xml'));

%mgat3
 mgat3    = GTEnz([2;4;1;143]);
 mgat3.resfuncgroup  = residueMap.allresidues('GlcNAc');
 manResType   = residueMap.allresidues('Man');
 manBond   = GlycanBond('4','1');
 mgat3.resAtt2FG  = manResType;
 mgat3.linkresAtt2FG  = struct('bond', manBond,'anomer','b');
 glcnacbond   = GlycanBond('4','1');
 mgat3.linkFG   = struct('anomer','b','bond',glcnacbond);
 m3gn    = glycanMLread('m3gn.glycoct_xml');
 mgat3.substMinStruct= m3gn; 
 mgat3.substNABranch  = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
 mgat3.targetBranch = glycanMLread('mgat3targetbranch.glycoct_xml');
 mgat3.substNAResidue= residueMap.allresidues('Gal');
 
 %mgat4
 mgat4    = GTEnz([2;4;1;145]);
 mgat4.resfuncgroup  = residueMap.allresidues('GlcNAc');
 manResType   = residueMap.allresidues('Man');
 manBond   = GlycanBond('3','1');
 mgat4.resAtt2FG  = manResType;
 mgat4.linkresAtt2FG  = struct('bond', manBond,'anomer','a');
 glcnacbond   = GlycanBond('4','1');
 mgat4.linkFG   = struct('anomer','b','bond',glcnacbond);
 gn2m3gn    = glycanMLread('gnm3gn.glycoct_xml');
 mgat4.substMinStruct= gn2m3gn; 
 mgat4.substNABranch=CellArrayList;
 mgat4.substNABranch.add(glycanMLread('mgat4NAbranch1.glycoct_xml'));
 mgat4.substNABranch.add(glycanMLread('mgat4NAbranch2.glycoct_xml'));
 
 %mgat5
 mgat5    = GTEnz([2;4;1;155]);
 mgat5.resfuncgroup  = residueMap.allresidues('GlcNAc');
 manResType   = residueMap.allresidues('Man');
 manBond   = GlycanBond('6','1');
 mgat5.resAtt2FG  = manResType;
 mgat5.linkresAtt2FG  = struct('bond', manBond,'anomer','a');
 glcnacbond   = GlycanBond('6','1');
 mgat5.linkFG   = struct('anomer','b','bond',glcnacbond);
 gn2m3gn    = glycanMLread('mgat5minimal.glycoct_xml');
 mgat5.substMinStruct= gn2m3gn; 
 mgat5.substNABranch=CellArrayList;
 mgat5.substNABranch.add(glycanMLread('mgat5NAbranch1.glycoct_xml'));
 mgat5.substNABranch.add(glycanMLread('mgat5NAbranch2.glycoct_xml'));
 
 %manii
 manii = GHEnz([3;2;1;114]);
 manii.resfuncgroup =residueMap.allresidues('Man');
 manii.linkFG.anomer = 'a';
 manBond(1,1) = GlycanBond('3','1');
 manBond(2,1) = GlycanBond('6','1');
 manii.linkFG.bond = manBond;
 manii.resAtt2FG  = residueMap.allresidues('Man');
 manii.substMinStruct=glycanMLread('maniiminstruct.glycoct_xml');
 manii.substMaxStruct=glycanMLread('maniimaxstruct.glycoct_xml');
 
 %mani
 mani    = GTEnz([3;2;1;113]);
 mani.resfuncgroup  = residueMap.allresidues('Man');
 mani.resAtt2FG     = residueMap.allresidues('Man');
 manBond   = GlycanBond('2','1');
 mani.linkFG   = struct('anomer','a','bond',manBond);
 manunknownbond=GlycanBond('?','1');
 mani.linkresAtt2FG  = struct('bond', manunknownbond,'anomer','a');
 mani.substMinStruct=glycanMLread('maniminstruct.glycoct_xml');
 mani.substMaxStruct=glycanMLread('manimaxstruct.glycoct_xml');
 
 %ST6GalI
 ST6GalI=GTEnz([2;4;99;1]);
 ST6GalI.resfuncgroup=residueMap.allresidues('NeuAc');
 ST6GalI.resAtt2FG=residueMap.allresidues('Gal');
 fucbond=GlycanBond('6','2');
 ST6GalI.linkFG=struct('anomer','a','bond',fucbond);
 gnbond=GlycanBond('4','1');
 ST6GalI.linkresAtt2FG=struct('bond',gnbond,'anomer','b');
 ST6GalI.isTerminalTarget=true;
 lacnac2m3gn=glycanMLread('st6gal1minstruct.glycoct_xml');
 ST6GalI.substMinStruct=lacnac2m3gn;
 galglcnac=glycanMLread('galglcnac.glycoct_xml');
 ST6GalI.targetbranchcontain=galglcnac;
 ST6GalI.substNABranch=glycanMLread('st6gal1NAbranch.glycoct_xml');
 
 %Fut8
 Fut8=GTEnz([2;4;1;68]);
 Fut8.resfuncgroup=residueMap.allresidues('Fuc');
 Fut8.resAtt2FG=residueMap.allresidues('GlcNAc');
 fucbond=GlycanBond('6','1');
 Fut8.linkFG=struct('anomer','a','bond',fucbond);
 gnbond=GlycanBond('?','?');
 Fut8.linkresAtt2FG=struct('bond',gnbond,'anomer','b');
 gnm3=glycanMLread('fucminstruct.glycoct_xml');
 Fut8.substMinStruct=gnm3;
 Fut8.targetBranch=glycanMLread('fuctargetbranch.glycoct_xml');
 
 %b4GalI
 b4GalI=GTEnz([2;4;1;38]);
 b4GalI.resfuncgroup=residueMap.allresidues('Gal');
 b4GalI.resAtt2FG=residueMap.allresidues('GlcNAc');
 fucbond=GlycanBond('4','1');
 b4GalI.linkFG=struct('anomer','b','bond',fucbond);
 gnbond=GlycanBond('?','1');
 b4GalI.linkresAtt2FG=struct('bond',gnbond,'anomer','b');
 gn2m3gn=glycanMLread('b4galiminstruct.glycoct_xml');
 b4GalI.substMinStruct=gn2m3gn;
 b4GalI.targetNABranch=glycanMLread('b4galitargetNAbranch.glycoct_xml');
 
 %iGnt
 iGnt=GTEnz([2;4;1;149]);
 iGnt.resfuncgroup=residueMap.allresidues('GlcNAc');
 galResType=residueMap.allresidues('Gal');
 iGnt.resAtt2FG=galResType;
 manBond=GlycanBond('4','1');
 iGnt. linkresAtt2FG=struct('bond',manBond,'anomer','b');
 glcnacbond=GlycanBond('3','1');
 iGnt.linkFG=struct('anomer','b','bond',glcnacbond);
 iGnt.substMinStruct=glycanMLread('igntminstruct.glycoct_xml');
 iGnt.targetbranchcontain=glycanMLread('ignttargetbranchcontain.glycoct_xml');

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
 a1 = GlycanSpecies(glycanMLread('1579.8.glycoct_xml'));
 a2 = GlycanSpecies(glycanMLread('3142.5a.glycoct_xml'));
 a3 = GlycanSpecies(glycanMLread('3142.5b.glycoct_xml'));
 a4 = GlycanSpecies(glycanMLread('3142.5c.glycoct_xml'));
 a5 = GlycanSpecies(glycanMLread('3142.5d.glycoct_xml'));
 a6 = GlycanSpecies(glycanMLread('3142.5e.glycoct_xml'));
 a7 = GlycanSpecies(glycanMLread('3142.5f.glycoct_xml'));
 a8 = GlycanSpecies(glycanMLread('3142.5g.glycoct_xml'));
 a9 = GlycanSpecies(glycanMLread('3142.5h.glycoct_xml'));
 a10 = GlycanSpecies(glycanMLread('3142.5i.glycoct_xml'));
 a11 = GlycanSpecies(glycanMLread('3142.5j.glycoct_xml'));
 a12 = GlycanSpecies(glycanMLread('m3gngn.glycoct_xml'));
 
 glycanArray = CellArrayList;
%  glycanArray.add(a1);
%  glycanArray.add(a2);
%  glycanArray.add(a3);
 glycanArray.add(a4);
 glycanArray.add(a12);
%  glycanArray.add(a5);
%  glycanArray.add(a6);
%  glycanArray.add(a7);
%  glycanArray.add(a8);
%  glycanArray.add(a9);
%  glycanArray.add(a10);
%  glycanArray.add(a11);


 
 %Perform reaction
 displayOptions = displayset('showMass',true,'showLinkage',true,'showRedEnd',true);
 fprintf(1,'Input of glycan product structure is \n');
%  glycanViewer(a1.glycanStruct,displayOptions);
%  glycanViewer(a2.glycanStruct,displayOptions);
%  glycanViewer(a3.glycanStruct,displayOptions);
 glycanViewer(a4.glycanStruct,displayOptions);
 glycanViewer(a12.glycanStruct,displayOptions);
%  glycanViewer(a5.glycanStruct,displayOptions);
%  glycanViewer(a6.glycanStruct,displayOptions);
%  glycanViewer(a7.glycanStruct,displayOptions);
%  glycanViewer(a8.glycanStruct,displayOptions);
%  glycanViewer(a9.glycanStruct,displayOptions);
%  glycanViewer(a10.glycanStruct,displayOptions);
%  glycanViewer(a11.glycanStruct,displayOptions);
 [isPath,nlinkedpath]=inferGlyConnPath(glycanArray,enzArray);
 fprintf(1,'Inferred network is shown below:\n');
 glycanPathViewer(nlinkedpath);