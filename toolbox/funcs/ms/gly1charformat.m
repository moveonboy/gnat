function glycan1letstring=gly1charformat(glycanresiduestring)
%gly3to1letter 
%Example 1:
%      glycan3letterstring='Hex5HexNAc2'; % Hex5HexNAc2
%      glycan1letterstring=gly1charformat(glycanstring);%
%      glycan3letstring = 'Hex5HexNAc2';
%      glycanresiduestring ='Hex25HexNAc25';


% Date Lastly Updated:
% Author: 
glycan1letstring =[];

% HexNAC -h
nHexNAC  =  regexp(glycanresiduestring,'(HexNAc\d+)','match');
numHexNAc     =  str2num(nHexNAC{1,1}(7:end));
for i = 1 : numHexNAc
    glycan1letstring = strcat(glycan1letstring,'n');
end

% Hex
nHex       =  regexp(glycanresiduestring,'(Hex\d+)','match');
numHex     =  str2num(nHex{1,1}(4:end));
for i = 1 : numHex
    glycan1letstring = strcat(glycan1letstring,'h');
end

% Sailic Acid NeuAc
% nNeuAc=findstr('s',gly);
nNeuAc             =  regexp(glycanresiduestring,'(NeuAc\d+)','match');

if(~isempty(nNeuAc))
numNeuAc       =  str2num(nNeuAc{1,1}(6:end));
for i = 1 : numNeuAc
    glycan1letstring = strcat(glycan1letstring,'s');
end

end

% Sailic Acid NeuGc
% nNeuAc=findstr('g',gly);
nNeuGc              =  regexp(glycanresiduestring,'(NeuGc\d+)','match') ;
if(~isempty(nNeuGc))
numNeuGc       =  str2num(nNeuGc{1,1}(6:end));
for i = 1 : numNeuGc
    glycan1letstring = strcat(glycan1letstring,'g');
end
end

% Fucose f
nFuc               =  regexp(glycanresiduestring,'(Fuc\d+)','match') ;  
if(~isempty(nFuc))

numFuc         =  str2num(nFuc{1,1}(4:end));
for i = 1 : numFuc
    glycan1letstring = strcat(glycan1letstring,'f');
end

end

% pentose
nPen           = regexp(glycanresiduestring,'(Pen\d+)','match') ;
if(~isempty(nPen))

numPen         =  str2num(nPen{1,1}(4:end)); 
for i = 1 : numPen
    glycan1letstring = strcat(glycan1letstring,'x');
end

end

end

