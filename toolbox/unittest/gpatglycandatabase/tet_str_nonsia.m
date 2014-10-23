%% Enzyme definition
% MAN II, MGAT 1,2,3,4,5,GalT, FucT, SiaT enzymes
% are loaded from a pre-defined local database.

enzdbmatfilename   = 'glyBigenzDB.mat';
enzdb = enzdbmatLoad(enzdbmatfilename);
mgat1 = enzdb('mgat1');
mgat2 = enzdb('mgat2');
mgat3 = enzdb('mgat3');
mgat4 = enzdb('mgat4');
mgat5 = enzdb('mgat5');
mani = enzdb('mani');
manii  = enzdb('manii');
galt     = enzdb('galt');
gnte    = enzdb('gnte');
siaT     = enzdb('siaT');
fucT8 = enzdb('fucT8');
FUT2 = enzdb('FUT2');
TRAA = enzdb('TRAA');
TRAB = enzdb('TRAB');

%% Definition of Glycan Structures
% Initial glycans are defined based on 20 peaks and stored in substrateArray variable.

gDispOption1    = displayset('showmass',true,'showLinkage',true,...
    'showRedEnd',true);

s1species    = GlycanSpecies(glycanMLread('M3gn.glycoct_xml')) ;
glycanViewer(s1species.glycanStruct,gDispOption1);

s2species    = GlycanSpecies(glycanMLread('tetra_AAAA.glycoct_xml')) ;
glycanViewer(s2species.glycanStruct,gDispOption1);
% 
% s3species    = GlycanSpecies(glycanMLread('tetra_AAAB.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s4species    = GlycanSpecies(glycanMLread('tetra_AAAH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s5species    = GlycanSpecies(glycanMLread('tetra_AABB.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s6species    = GlycanSpecies(glycanMLread('tetra_AABH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s7species    = GlycanSpecies(glycanMLread('tetra_AAHH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s8species    = GlycanSpecies(glycanMLread('tetra_ABBB.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s9species    = GlycanSpecies(glycanMLread('tetra_ABBH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s10species    = GlycanSpecies(glycanMLread('tetra_ABHH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s11species    = GlycanSpecies(glycanMLread('tetra_AHHH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s12species    = GlycanSpecies(glycanMLread('tetra_BBBB.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s13species    = GlycanSpecies(glycanMLread('tetra_BBBH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s14species    = GlycanSpecies(glycanMLread('tetra_BBHH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s15species    = GlycanSpecies(glycanMLread('tetra_BHHH.glycoct_xml')) ;
% % % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 
% s16species    = GlycanSpecies(glycanMLread('tetra_HHHH.glycoct_xml')) ;
% % % glycanViewer(s1species.glycanStruct,gDispOption1);
% % 


substrateArray = CellArrayList;
substrateArray.add(s1species);
substrateArray.add(s2species);
% substrateArray.add(s3species);
% substrateArray.add(s4species);
% substrateArray.add(s5species);
% substrateArray.add(s6species);
% substrateArray.add(s7species);
% substrateArray.add(s8species);
% substrateArray.add(s9species);
% substrateArray.add(s10species);
% substrateArray.add(s11species);
% substrateArray.add(s12species);
% substrateArray.add(s13species);
% substrateArray.add(s14species);
% substrateArray.add(s15species);
% substrateArray.add(s16species);

%% Store enzymes and glycan in CellArrayList variables
% enzArray are created to store  enzymes which might act on substrates.
enzArray = CellArrayList;
enzArray.add(mani);
enzArray.add(manii);
enzArray.add(mgat1);
enzArray.add(mgat2);
enzArray.add(mgat3);
enzArray.add(mgat4);
enzArray.add(mgat5);
enzArray.add(galt);
enzArray.add(fucT8);
enzArray.add(gnte);
enzArray.add(siaT);
enzArray.add(FUT2);
enzArray.add(TRAA);
enzArray.add(TRAB);

%% Pathway reconstruction using connection inferrence
%  inferGlyConnPath command is used to construct a
%    pathway and the pathway can be visualized using glycanPathViewer
%    command.  The reconstructed network has 288 reactions and 151 species.
%    This construction takes about 5 mins to finish on Core 7 computer.
%
[isPath,nglycanpath]=inferGlyConnPath(substrateArray, enzArray);

if(isPath)
    glycanPathViewer(nglycanpath);
    fprintf(1,'the generated network has: \n');
    fprintf(1,'    %i  reactions\n',nglycanpath.getNReactions);
    fprintf(1,'    %i  species\n',   nglycanpath.getNSpecies);
else
    fprintf(1,'no path is found\n');
end
glycanpath2 = removeIsomerStructk(nglycanpath);


