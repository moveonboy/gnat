%% GNAT Demo 5
% This demo shows how to construct an N-linked glycosylation pathway using glycans annoated from
%   mass spectra data.  The inputs are glycan structures annoated from 24 peaks and 11 enzymes.

%% Enzyme definition
% MAN II, MGAT 1,2,3,4,5,GalT, FucT, SiaT enzymes
% are loaded from a pre-defined local database.

enzdbmatfilename   = 'glyenzDB.mat';
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
enzViewer(fucT8);

%% Definition of Glycan Structures
% Initial glycans are defined based on 20 peaks and stored in substrateArray variable.

gDispOption1    = displayset('showmass',true,'showLinkage',true,...
    'showRedEnd',true);
% s2species    = GlycanSpecies(glycanMLread('Man9.glycoct_xml')) ;
% glycanViewer(s1species.glycanStruct,gDispOption1);
% s2species    = GlycanSpecies(glycanMLread('s2_1590.glycoct_xml')) ;
% glycanViewer(s2species.glycanStruct,gDispOption1);
s1species    = GlycanSpecies(glycanMLread('Man9.glycoct_xml')) ;
glycanViewer(s1species.glycanStruct,gDispOption1);

s2species    = GlycanSpecies(glycanMLread('tri_sia_fuc.glycoct_xml')) ;
glycanViewer(s2species.glycanStruct,gDispOption1);

s3species    = GlycanSpecies(glycanMLread('bi_sia.glycoct_xml')) ;
glycanViewer(s3species.glycanStruct,gDispOption1);

% s3v2species    = GlycanSpecies(glycanMLread('s3v2_1620.glycoct_xml')) ;
% glycanViewer(s3v2species.glycanStruct,gDispOption1);
% 
s4species    = GlycanSpecies(glycanMLread('bi_sia_fuc.glycoct_xml')) ;
glycanViewer(s4species.glycanStruct,gDispOption1);
% 
s5species    = GlycanSpecies(glycanMLread('tetra_sia.glycoct_xml')) ;
glycanViewer(s5species.glycanStruct,gDispOption1);
% 
s6species    = GlycanSpecies(glycanMLread('tri_sia.glycoct_xml')) ;
glycanViewer(s6species.glycanStruct,gDispOption1);

% 
s7species    = GlycanSpecies(glycanMLread('tetra_sia_fuc.glycoct_xml')) ;
glycanViewer(s7species.glycanStruct,gDispOption1);

substrateArray = CellArrayList;
substrateArray.add(s1species);
% substrateArray.add(s2species);
% substrateArray.add(s3species);
% % substrateArray.add(s3v2species);
% substrateArray.add(s4species);
% substrateArray.add(s5species);
% substrateArray.add(s6species);
% substrateArray.add(s6v2species);
substrateArray.add(s7species);
% substrateArray.add(s8v1species);
% substrateArray.add(s8v2species);
% substrateArray.add(s9species);
% substrateArray.add(s10species);
% % substrateArray.add(s11species);
% substrateArray.add(s12species);
% % substrateArray.add(s13species);
% substrateArray.add(s14v1species);
% substrateArray.add(s14v2species);
% 
% substrateArray.add(s15species);
% substrateArray.add(s16v1species);
% substrateArray.add(s16v2species);
% 
% substrateArray.add(s17v1species);
% substrateArray.add(s17v2species);
% 
% substrateArray.add(s18v1species);
% substrateArray.add(s18v2species);
% 
% substrateArray.add(s19v1species);
% substrateArray.add(s19v2species);
% substrateArray.add(s19v3species);
% substrateArray.add(s19v4species);
% substrateArray.add(s19v5species);
% substrateArray.add(s19v6species);
% 
% substrateArray.add(s20v1species);
% substrateArray.add(s20v2species);
% 
% substrateArray.add(s21species);
% substrateArray.add(s22v1species);
% substrateArray.add(s22v2species);
% substrateArray.add(s22v3species);
% substrateArray.add(s22v4species);
% substrateArray.add(s22v5species);
% substrateArray.add(s22v6species);
% 
% substrateArray.add(s23species);
% substrateArray.add(s24v1species);
% substrateArray.add(s24v2species);
% substrateArray.add(s24v3species);
% substrateArray.add(s24v4species);
% substrateArray.add(s24v5species);
% substrateArray.add(s24v6species);

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

%% Pathway reconstruction using connection inferrence
%  inferGlyConnPath command is used to construct a
%    pathway and the pathway can be visualized using glycanPathViewer
%    command.  The reconstructed network has 288 reactions and 151 species.
%    This construction takes about 5 mins to finish on Core 7 computer.

[isPath,nglycanpath]=inferGlyConnPath(substrateArray, enzArray,'iterativedisp',true);

if(isPath)
    glycanPathViewer(nglycanpath);
    fprintf(1,'the generated network has: \n');
    fprintf(1,'    %i  reactions\n',nglycanpath.getNReactions);
    fprintf(1,'    %i  species\n',   nglycanpath.getNSpecies);
else
    fprintf(1,'no path is found\n');
end
% glycanpath2 = removeIsomerStructk(nglycanpath);


