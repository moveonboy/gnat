%% Enzyme definition 
%  MANI, II, MGAT 1,2,3,4,5, GalT,IGNT,FucT,SiaT enzymes
%   are defined. 
residueMap                      = load('residueTypes.mat');
manResType                      = residueMap.allresidues('Man');
m3gn                            = glycanMLread('m3gn.glycoct_xml');

%% define sia T
siaT                            = GTEnz([2;4;99;6]);
siaT.isTerminalTarget           = true;
siaT.resfuncgroup               = residueMap.allresidues('NeuAc');
siaTbond                        = GlycanBond('3','2');
siaT.linkFG                     = struct('anomer','a','bond',siaTbond);
galResType                      = residueMap.allresidues('Gal');
galBond                         = GlycanBond('4','1');
siaT.resAtt2FG                  = galResType;
siaT.linkresAtt2FG              = struct('bond', galBond,'anomer','b');
siaT.substNABranch              = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
siaT.targetbranchcontain        = glycanMLread('gntetargetbranch.glycoct_xml');
enzViewer(siaT);

%% Define Fuc T
fucT8                           = GTEnz([2;4;1;68]);
fucT8.isTerminalTarget          = false;
fucT8.resfuncgroup              = residueMap.allresidues('Fuc');
fuctbond                        = GlycanBond('6','1');
fucT8.linkFG                    = struct('anomer','a','bond',fuctbond);
glcnacResType                   = residueMap.allresidues('GlcNAc');
glcnacBond                      = GlycanBond('?','?');
fucT8.resAtt2FG                 = glcnacResType;
fucT8.linkresAtt2FG             = struct('bond', glcnacBond,'anomer','?');
fucT8.targetNABranch            = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
fucT8.substNAResidue            = residueMap.allresidues('Gal');
fucT8.substMinStruct            = m3gn;
%fucT8.targetBranch             = glycanMLread('fuctTargetbranch.glycoct_xml');
enzViewer(fucT8);

%% Define MGAT1 
mgat1                            = GTEnz([2;4;1;101]);
mgat1.resfuncgroup               = residueMap.allresidues('GlcNAc');
glcnacbond                       = GlycanBond('2','1');
mgat1.linkFG                     = struct('anomer','b','bond',glcnacbond);
manBond                          = GlycanBond('3','1');
mgat1.resAtt2FG                  = manResType;
mgat1.linkresAtt2FG              = struct('bond', manBond,'anomer','a');
mgat1.substMinStruct             = glycanMLread('M5.glycoct_xml');
mgat1.substMaxStruct             = glycanMLread('M5.glycoct_xml');
mgat1.targetBranch               = glycanMLread('M5_lowerbranch.glycoct_xml');
enzViewer(mgat1);

%% Define IGNT 
gnte                             = GTEnz([2;4;1;149]);
gnte.isTerminalTarget            = true;
gnte.resfuncgroup                = residueMap.allresidues('GlcNAc');
galtbond                         = GlycanBond('3','1');
gnte.linkFG                      = struct('anomer','b','bond',galtbond);
galResType                       = residueMap.allresidues('Gal');
galBond                          = GlycanBond('4','1');
gnte.resAtt2FG                   = galResType;
gnte.linkresAtt2FG               = struct('bond', galBond,'anomer','b');
gnte.targetNABranch              = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
% galt.substNABranch=glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
%gnte.targetbranchcontain      = glycanMLread('gntetargetbranch.glycoct_xml');
enzViewer(gnte);

%% Define MGAT2 
mgat2                             = GTEnz([2;4;1;143]);
residueMap                        = load('residueTypes.mat');
mgat2.resfuncgroup                = residueMap.allresidues('GlcNAc');
glcnacbond                        = GlycanBond('2','1');
mgat2.linkFG                      = struct('anomer','b','bond',glcnacbond);
manResType                        = residueMap.allresidues('Man');
manBond                           = GlycanBond('6','1');
mgat2.resAtt2FG                   = manResType;
mgat2.linkresAtt2FG               = struct('bond', manBond,'anomer','a');
m3gn                              = glycanMLread('m3gn.glycoct_xml');
mgat2.isTerminalTarget            = true;
mgat2.substMinStruct              = m3gn;
mgat2.targetBranch                = glycanMLread('mgat2actingbranch.glycoct_xml');
mgat2.substNABranch               = CellArrayList;
mgat2.substNABranch.add(glycanMLread('NGlycanBisectGlcNAc.glycoct_xml'));
mgat2.substNABranch.add(glycanMLread('mgat2substrateNAbranch.glycoct_xml'));
enzViewer(mgat2)

%% Define MANI 
mani                              = GHEnz([3;2;1;113]);
mani.resfuncgroup                 = manResType;
mani.linkFG.anomer                ='a';
manBond                           = GlycanBond('2','1');
mani.linkFG.bond                  = manBond;
mani.resAtt2FG                    = manResType;
manAttachBond                     = GlycanBond('?','?');
mani.linkresAtt2FG                = struct('bond', manAttachBond,'anomer','?');
mani.prodMinStruct                = glycanMLread('M5.glycoct_xml'); 
mani.substMaxStruct               = glycanMLread('M9.glycoct_xml');
mani.substNAResidue               = residueMap.allresidues('Gal');
enzViewer(mani);

%% Define MANII
manii                            = GHEnz([3;2;1;114]);
manii.resfuncgroup               = manResType;
manii.linkFG.anomer              = 'a'; 
manBond(1,1)                     = GlycanBond('3','1');
manBond(2,1)                     = GlycanBond('6','1');
manii.linkFG.bond                = manBond;
manii.resAtt2FG                  = manResType;
manii.prodMinStruct              = glycanMLread('m3gn.glycoct_xml'); % mannosidase II can not work on M3 glycan
manii.substNABranch              = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
manii.substMaxStruct             = glycanMLread('m5gn.glycoct_xml');
enzViewer(manii)

%% Definition and visualization of MGAT3 enzyme 
mgat3                            = GTEnz([2;4;1;143]);
mgat3.resfuncgroup               = residueMap.allresidues('GlcNAc');
manResType                       = residueMap.allresidues('Man');
manBond                          = GlycanBond('4','1');
mgat3.resAtt2FG                  = manResType;
mgat3.linkresAtt2FG              = struct('bond', manBond,'anomer','b');
glcnacbond                       = GlycanBond('4','1');
mgat3.linkFG                     = struct('anomer','b','bond',glcnacbond);
m3gn                             = glycanMLread('m3gn.glycoct_xml');
mgat3.substMinStruct             = m3gn;
mgat3.substNABranch              = glycanMLread('NGlycanBisectGlcNAc.glycoct_xml');
mgat3.targetBranch               = glycanMLread('mgat3targetbranch.glycoct_xml');
mgat3.substNAResidue             = residueMap.allresidues('Gal');
enzViewer(mgat3)

%% Definition and visualization of MGAT4 enzyme 
mgat4                             = GTEnz([2;4;1;145]);
mgat4.resfuncgroup                = residueMap.allresidues('GlcNAc');
manResType                        = residueMap.allresidues('Man');
manBond                           = GlycanBond('3','1');
mgat4.resAtt2FG                   = manResType;
mgat4.linkresAtt2FG               = struct('bond', manBond,'anomer','a');
glcnacbond                        = GlycanBond('4','1');
mgat4.linkFG                      = struct('anomer','b','bond',glcnacbond);
m3gn                              = glycanMLread('m3gn.glycoct_xml');
mgat4.substMinStruct              = m3gn;
mgat4.targetBranch                =  glycanMLread('mgat4targetbranch.glycoct_xml');
mgat4.substNABranch               = CellArrayList;
mgat4.substNABranch.add(glycanMLread('NGlycanBisectGlcNAc.glycoct_xml'));
mgat4.substNABranch.add(glycanMLread('mgat4subsNAbranch.glycoct_xml'));
mgat4.substNAResidue              = residueMap.allresidues('Gal');
%enzViewer(mgat4)

%% Definition and visualization of MGAT5 enzyme 
mgat5                             = GTEnz([2;4;1;155]);
mgat5.resfuncgroup                = residueMap.allresidues('GlcNAc');
manResType                        = residueMap.allresidues('Man');
manBond                           = GlycanBond('6','1');
mgat5.resAtt2FG                   = manResType;
mgat5.linkresAtt2FG               = struct('bond', manBond,'anomer','a');
glcnacbond                        = GlycanBond('6','1');
mgat5.linkFG                      = struct('anomer','b','bond',glcnacbond);
m3gn                              = glycanMLread('m3gn.glycoct_xml');
%mgat5.substMinStruct  = m3gn;
mgat5.targetBranch                = glycanMLread('mgat5targetbranch.glycoct_xml');
mgat5.substMinStruct              = glycanMLread('mgat5substMinStruct.glycoct_xml');
mgat5.substNABranch               = CellArrayList;
mgat5.substNABranch.add(glycanMLread('NGlycanBisectGlcNAc.glycoct_xml'));
mgat5.substNABranch.add(glycanMLread('mgat5substrateNAbranch.glycoct_xml'));
mgat5.substNAResidue              = residueMap.allresidues('Gal');

%% Define Beta4GalT 
beta4galt                              = GTEnz([2;4;1;38]);
beta4galt.isTerminalTarget             = true;
beta4galt.resfuncgroup                 = residueMap.allresidues('Gal');
glcnacResType                     = residueMap.allresidues('GlcNAc');
glcnacBond                        = GlycanBond('?','1');
beta4galt.resAtt2FG                    = glcnacResType;
beta4galt.linkresAtt2FG                = struct('bond', glcnacBond,'anomer','b');
galtbond                          = GlycanBond('4','1');
beta4galt.linkFG                       = struct('anomer','b','bond',galtbond);

%% Define IGNT
ignt                              = GTEnz([2;4;1;149]);
ignt.isTerminalTarget             = true;
ignt.resfuncgroup                 = residueMap.allresidues('GlcNAc');
galResType                        = residueMap.allresidues('Gal');
galBond                           = GlycanBond('4','1');
ignt.resAtt2FG                    = galResType;
ignt.linkresAtt2FG                = struct(...
    'bond', galBond,'anomer','b');
glcnacbond                        = GlycanBond('3','1');
ignt.linkFG                       = struct('anomer','b','bond',glcnacbond);

%% Define FUCT7
fucT7                      = GTEnz([2;4;1;152]);
fucT7.isTerminalTarget     = false;
fucT7.resfuncgroup         = residueMap.allresidues('Fuc');
fuctbond                   = GlycanBond('4','1');
fucT7.linkFG               = struct('anomer','a','bond',fuctbond);

glcnacResType              = residueMap.allresidues('GlcNAc');
glcnacBond                 = GlycanBond('?','1');
fucT7.resAtt2FG            =  glcnacResType;
fucT7.linkresAtt2FG        =  struct('bond', glcnacBond,'anomer','?');
fucT7.substNAResidue       = residueMap.allresidues('Fuc');
fucT7.substMinStruct       =  glycanMLread('fucT7substmin.glycoct_xml');

%% Define FUT2
FUT2                       = GTEnz([2;4;1;69]); %not sure if EC # is correct
FUT2.isTerminalTarget      = true;
FUT2.resfuncgroup          = residueMap.allresidues('Fuc');
FUT2bond                   = GlycanBond('2','1');
FUT2.linkFG                = struct('anomer','a','bond',FUT2bond);
galResType                 = residueMap.allresidues('Gal');
galBond                    = GlycanBond('4','1');
FUT2.resAtt2FG             = galResType;
FUT2.linkresAtt2FG         = struct('bond', galBond,'anomer','b');
% FUT2.targetNABranch        = glycanMLread('FUT2NAbranch.glycoct_xml');
% FUT2.substNAResidue      = residueMap.allresidues('Gal');
% FUT2.substMinStruct        = m3gn;
FUT2.targetbranchcontain      = glycanMLread('FUT2Targetbranch.glycoct_xml');
% enzViewer(FUT2);

%% Define A Transferase
TRAA                             = GTEnz([2;4;1;40]);
TRAA.isTerminalTarget            = false;
TRAA.resfuncgroup                = residueMap.allresidues('GalNAc');
galResType                       = residueMap.allresidues('Gal');
TRAABond                         = GlycanBond('3','1');
TRAA.resAtt2FG                   = galResType;
galBond                          = GlycanBond('4','1');
TRAA.linkresAtt2FG               = struct('bond',galBond,'anomer','b');
TRAA.linkFG                      = struct('anomer','a','bond',TRAABond);
TRAA.targetbranchcontain         = glycanMLread('FUT2Targetbranch.glycoct_xml');

%% Define B Transferase
TRAB                             = GTEnz([2;4;1;37]);
TRAB.isTerminalTarget            = false;
TRAB.resfuncgroup                = residueMap.allresidues('Gal');
galResType                       = residueMap.allresidues('GalNAc');
TRABBond                         = GlycanBond('3','1');
TRAB.resAtt2FG                   = galResType;
galBond                          = GlycanBond('4','1');
TRAB.linkresAtt2FG               = struct('bond', galBond,'anomer','b');
TRAB.linkFG                      = struct('anomer','a','bond',TRABBond);
TRAB.targetbranchcontain         = glycanMLread('FUT2Targetbranch.glycoct_xml');

%% PP-GalNAc-T
ppGalNAcT = GTEnz([2;4;1;41]);
ppGalNAcT.isTerminalTarget        = false;
ppGalNAcT.resfuncgroup            = residueMap.allresidues('GalNAc');
ppGalNAcT.resAtt2FG               = residueMap.allresidues('freeEnd');
ppGalNAcTBond                     = GlycanBond('?','1');
ppGalNAcT.linkFG                  = struct('anomer','a','bond',ppGalNAcTBond);

%% ST6GalNAc
st6galnac                  = GTEnz([2;4;99;3]);
st6galnac.resfuncgroup     = residueMap.allresidues('NeuAc');
galNAcResType              = residueMap.allresidues('GalNAc');
st6galnac.resAtt2FG        = galNAcResType;
galNAcBond                 = GlycanBond('?','1');
st6galnac.linkresAtt2FG    = struct('bond', galNAcBond,'anomer','?');
neuacgalnacbond            = GlycanBond('6','2');
st6galnac.linkFG           = struct('anomer','a','bond',neuacgalnacbond);

%% C1GalT-1 
c1galt1 = GTEnz([2;4;1;122]);
c1galt1.isTerminalTarget   = true;
c1galt1.resfuncgroup       = residueMap.allresidues('Gal');
c1galt1.resAtt2FG          = residueMap.allresidues('GalNAc');
galNAcBond                  = GlycanBond('?','1');
c1galt1.linkresAtt2FG      = struct('bond', galNAcBond,'anomer','?');
galbond                     = GlycanBond('3','1');
c1galt1.linkFG             = struct('anomer','b','bond',galbond);

%% ST3Gal1 
st3gal1                  = GTEnz([2;4;99;4]);
st3gal1.isTerminalTarget = true;
st3gal1.resfuncgroup     = residueMap.allresidues('NeuAc');
st3gal1.resAtt2FG        = residueMap.allresidues('Gal');
galbond                  = GlycanBond('3','1');
st3gal1.linkresAtt2FG    = struct('bond',galbond,'anomer','b');
neuacgalnacbond          = GlycanBond('3','2');
st3gal1.linkFG           = struct('anomer','a','bond',neuacgalnacbond);

%% ST3GalIV
st3galIV                   = GTEnz([2;4;99;6]);
st3galIV.isTerminalTarget  = true;
st3galIV.resfuncgroup      = residueMap.allresidues('NeuAc');
galResType                 = residueMap.allresidues('Gal');
galBond                    = GlycanBond('4','1');
st3bond                    = GlycanBond('3','2');
st3galIV.linkFG            = struct('anomer','a','bond',st3bond);
st3galIV.resAtt2FG         = galResType;
st3galIV.linkresAtt2FG     = struct('bond', galBond,'anomer','b');

%% Define b3GlcNAcT elongating beta3GnT1
b3glcnacT                   = GTEnz([2;4;1;146]);
b3glcnacT.isTerminalTarget  = true;
b3glcnacT.resfuncgroup      = residueMap.allresidues('GlcNAc');
b3glcnacT.resAtt2FG         = residueMap.allresidues('Gal');
galBond                     = GlycanBond('3','1');
glcnacbond                  = GlycanBond('3','1');
b3glcnacT.linkFG            = struct('anomer','b','bond',glcnacbond);
b3glcnacT.linkresAtt2FG     = struct('bond', galBond,'anomer','b');

%% Define C2GNT beta6 
c2gnt                       = GTEnz([2;4;1;102]);
c2gnt.isTerminalTarget      = false;
c2gnt.resfuncgroup          = residueMap.allresidues('GlcNAc');
c2gnt.resAtt2FG             = residueMap.allresidues('GalNAc');
glcnacbond                  = GlycanBond('6','1');
galnacbond                  = GlycanBond('?','1');
c2gnt.linkFG                = struct('anomer','b','bond',glcnacbond);
c2gnt.linkresAtt2FG         = struct('bond', galnacbond,'anomer','?');
ct2ngtnastruct              = CellArrayList;
% ct2ngtnastruct.add(glycanMLread('core2_Oglycan.glycoct_xml'));
% ct2ngtnastruct.add(glycanMLread('sialylcore2_Oglycan.glycoct_xml'));
% ct2ngtnastruct.add(glycanMLread('sialylcore2_Oglycan.glycoct_xml'));
% ct2ngtnastruct.add(glycanMLread('sialylcore2_Oglycan.glycoct_xml'));
c2gnt.substNAStruct         =  glycanMLread('core2_Oglycan.glycoct_xml');

%% Define C3GNT 
c3gnt                      = GTEnz([2;4;1;102]);
c3gnt.isTerminalTarget     = true;
c3gnt.resfuncgroup         = residueMap.allresidues('GlcNAc');
c3gnt.resAtt2FG            = residueMap.allresidues('Gal');
glcnacBond                 = GlycanBond('3','1');
galnacbond                 = GlycanBond('?','1');
c3gnt.linkFG               = struct('anomer','b','bond',glcnacBond);
c3gnt.linkresAtt2FG        = struct('bond', galnacbond,'anomer','a');
c3gnt.substMinStruct       = glycanMLread('core2_Oglycan.glycoct_xml');

enzDB = containers.Map;
enzDB('mgat5')=mgat5;
enzDB('mgat4')=mgat4;
enzDB('mgat3')=mgat3;
enzDB('mgat2')=mgat2;
enzDB('mgat1')=mgat1;
enzDB('manii')=manii;
enzDB('mani') =mani;
enzDB('ignt')=ignt;

enzDB('fucT8')=fucT8;
enzDB('fucT7')=fucT7;
enzDB('siaT') =siaT;
enzDB('gnte') =gnte;
enzDB('galt') =beta4galt;
enzDB('st3galI')=st3gal1;
enzDB('st3galIV')=st3galIV;
enzDB('fucT7')=fucT7;
enzDB('FUT2')=FUT2;
enzDB('TRAA')=TRAA;
enzDB('TRAB')=TRAB;
enzDB('ppGalNAcT')=ppGalNAcT;
enzDB('TRAB')=TRAB;
enzDB('b3glcnacT')=b3glcnacT;
enzDB('c3gnt')=c3gnt;
enzDB('c2gnt')=c2gnt;
enzDB('c1galt1')=c1galt1;
enzDB('st6galnac')=st6galnac;
enzDB('ppGalNAcT')=ppGalNAcT;

save('glyBigenzDB.mat','enzDB');



