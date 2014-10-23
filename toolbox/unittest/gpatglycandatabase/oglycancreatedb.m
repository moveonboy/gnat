%% O-Glycan Database
% This demo shows how to construct an O-linked glycosylation pathway on
% PGSL-1 developed by Liu et al. (Bioinformatic, 2008).  The inputs are 12
% glycan structures and 5 enzymes.

%% Enzyme definition
% define beta1.4 GalTIV, beta 1,3 GlcNAcTgnte enzyme,
%  alpha2,3ST3Gal-I , alpha2,3ST3Gal-IV , and alpha1,3FT-VII
clc;clear;
enzdbmatfilename  = 'glyBigenzDB.mat';
enzDB             = enzdbmatLoad(enzdbmatfilename);
galtiv     = enzDB('galt');
ignt       = enzDB('ignt');
st3galI    = enzDB('st3galI');
st3galIV   = enzDB('st3galIV');
fucT7      = enzDB('fucT7');
b3GlcNAc3  = enzDB('b3glcnacT');
st6galnac  = enzDB('st6galnac');
c1galt1    = enzDB('c1galt1');
c2gnt      = enzDB('c2gnt');

enzArray  = CellArrayList;
enzArray.add(ignt);
enzArray.add(galtiv);
enzArray.add(st3galI);
enzArray.add(st3galIV);
enzArray.add(b3GlcNAc3);
enzArray.add(st6galnac);
enzArray.add(fucT7);
enzArray.add(c1galt1);
enzArray.add(c2gnt);


%%  Inputs of Glycan Structure
%  12 Glycan structure are constructed using Glycoworkbench toolbox and they
%   are stored as Glycoct xml format. These files can be imported into
%   GNAT.
gDispOption1    = displayset('showmass',true,'showLinkage',true,...
    'showRedEnd',true);
Tantigen        = GlycanSpecies(glycanMLread('tnAntigen_Oglycan.glycoct_xml'));
%Tantigen.glycanStruct.root.linkageChildren.child.anomer.symbol='a';

glycanViewer(Tantigen.glycanStruct,gDispOption1);
[numglycanprod,glycanprod]= inferGlyProd(Tantigen,st6galnac);
if(numglycanprod>0)
    for i = 1 : length(glycanprod)
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);
    end
end
sialylAntigen = glycanprod.get(1);

[numglycanprod,glycanprod]= inferGlyProd(Tantigen,c1galt1);
if(numglycanprod>0)
    for i = 1 : length(glycanprod)
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
    end
end
core1struct = glycanprod.get(1);

[numglycanprod,glycanprod]= inferGlyProd(core1struct,c2gnt);
if(numglycanprod>0)
    for i = 1 : length(glycanprod)        
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
    end
end
core2struct = glycanprod.get(1);

[numglycanprod,glycanprod]= inferGlyProd(core1struct,b3GlcNAc3);
if(numglycanprod>0)
    for i = 1 : length(glycanprod)        
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
    end
end
glcnacCore1=glycanprod.get(1);

[numglycanprod,glycanprod]= inferGlyProd(core1struct,st3galI);
if(numglycanprod>0)
    for i = 1 : length(glycanprod)        
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
    end
end
sialycore1glycan = glycanprod.get(1);

[numglycanprod,glycanprod]= inferGlyProd(sialycore1glycan,st6galnac);
if(numglycanprod>0)
    for i = 1 : length(glycanprod)        
        glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
    end
end
dualsialylcore1glycan=glycanprod.get(1);

% numsubstr=inferGlySubstr(dualsialylcore1glycan,st6galnac);

%% Test programme
% gDispOption1    = displayset('showmass',true,'showLinkage',true,...
%     'showRedEnd',true);
% Tantigen  = GlycanSpecies(glycanMLread('tnAntigen_Oglycan.glycoct_xml'));
% Tantigen.glycanStruct.root.linkageChildren.child.anomer.symbol='a';
% 
% glycanViewer(Tantigen.glycanStruct,gDispOption1);
% [numglycanprod,glycanprod]= inferGlyProd(Tantigen,st6galnac);
% if(numglycanprod>0)
%     for i = 1 : length(glycanprod)
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);
%     end
% end
% 
% [numglycanprod,glycanprod]= inferGlyProd(Tantigen,c1galt1);
% if(numglycanprod>0)
%     for i = 1 : length(glycanprod)
%         
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);
%         
%     end
% end
% 
% core1 = glycanprod.get(1);
% 
% [numglycanprod,glycanprod]= inferGlyProd(core1,c2gnt);
% if(numglycanprod>0)
%     for i = 1 : length(glycanprod)        
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
%     end
% end
% 
% [numglycanprod,glycanprod]= inferGlyProd(core1,b3GlcNAc3);
% if(numglycanprod>0)
%     for i = 1 : length(glycanprod)        
%         glycanViewer(glycanprod.get(i).glycanStruct,gDispOption1);        
%     end
% end

 s1species = GlycanSpecies(glycanMLread('OG1.glycoct_xml'));
 glycanViewer(s1species.glycanStruct,gDispOption1);

 s2species = GlycanSpecies(glycanMLread('OG2.glycoct_xml'));
 s4species = GlycanSpecies(glycanMLread('OG4.glycoct_xml'));
 s5species = GlycanSpecies(glycanMLread('OG5.glycoct_xml'));
s7species  = GlycanSpecies(glycanMLread('OG7.glycoct_xml'));
s10species = GlycanSpecies(glycanMLread('OG10.glycoct_xml'));
s11species = GlycanSpecies(glycanMLread('OG11.glycoct_xml'));
s13species = GlycanSpecies(glycanMLread('OG13.glycoct_xml'));
s14species = GlycanSpecies(glycanMLread('OG14.glycoct_xml'));
s17species = GlycanSpecies(glycanMLread('OG17.glycoct_xml'));
s18species = GlycanSpecies(glycanMLread('OG18.glycoct_xml'));
s19species = GlycanSpecies(glycanMLread('OG19.glycoct_xml'));
s20species = GlycanSpecies(glycanMLread('OG20.glycoct_xml'));
sialylewisspecies = GlycanSpecies(glycanMLread('sialylewisOglycan.glycoct_xml'));
glycanViewer(sialylewisspecies.glycanStruct,gDispOption1);

[numglycansubstr,glycansubst]= inferGlySubstr(sialylewisspecies,fucT7);
 if(numglycansubstr>0)
     for i = 1 : numglycansubstr
         glycanViewer(glycansubst.get(i).glycanStruct,gDispOption1);
     end
 end


%% Storage of enzymes and glycans in the CellArrayList variables
%  Two input variables are created as CellArrayList variables and store
%    a list of enzymes and glyan structure respectively.
%
glycanArray = CellArrayList;
glycanArray.add(s1species);
glycanArray.add(s2species);
glycanArray.add(s4species);
glycanArray.add(s5species);
glycanArray.add(s7species);
glycanArray.add(s10species);
glycanArray.add(s11species);
glycanArray.add(s13species);
glycanArray.add(s17species);
glycanArray.add(s14species);
glycanArray.add(s18species);
glycanArray.add(s19species);
glycanArray.add(s20species);
glycanArray.add(Tantigen);
glycanArray.add(sialylAntigen);
glycanArray.add(glcnacCore1);
glycanArray.add(sialycore1glycan);
glycanArray.add(dualsialylcore1glycan);
glycanArray.add(sialylewisspecies);

%% Network reconstruction from inputs of glycans and enzymes
%  An O-linked glycosylation network can be constructed from 12 defined
%   glycan structures and 5 enzymes using inferGlyConnPath command.
%
[isPath,oglycanpath]=inferGlyConnPath(glycanArray, enzArray);

%% Visualization of Reconstructed Network
if(isPath)
    glycanPathViewer(oglycanpath);
end
oglycanpath.toGlycanListinSmallGlyPep(1,'test.txt');