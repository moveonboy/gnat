function [glycanForm,varargout]= glycanFormula(gly,varargin)
%glycanFormula calculate glycan chemical formula
%
%   glycanForm = glycanFormula(glycanstring)  takes an input of the glycan 
%       string and outputs its chemical formula. Default options for the 
%      calculation of chemical formula are "methylated", "Na ionized", "showIsotopicDistPlot'. 
%
%   glycanForm = glycanFormula(glycanstring,options)  uses customized option
%       structure for the calcuation of chemical formula. The fields in the option
%      include 'methylation' ([true]/false),'ion'
%      ('none',['Na'],'K'),'showIsotopicDistPlot' (true/[false]).
%    
%    [glycanForm,md]=glycanFormula(glycanstring,options) returns isotopic
%    distribution md. 
%
% Example 1:
%      glycanstring='hhhhhnn'; % Hex5HexNAc2
%      glycanFormula(glycanstring);
%
% Example 2:
%      glycanstring='hhhhhnn'; % Hex5HexNAc2
%       [glycanForm,md]=glycanFormula(glycanstring);
%
% Example 3:
%      options.showIsotopicDistPlot =false;
%       glycanstring='hhhhhnn'; % Hex5HexNAc2
%       [glycanForm,md]=glycanFormula(glycanstring,options);
% 
%  See also glycanMolWt, mwconstant.

% Date Lastly Updated: 2/19/14
% Author: Gang Liu, 

% option setup 
if(length(varargin)==1)
    options = varargin{1};
else
    options = struct('methylation',true,'ion','Na');
end

if(isfield(options,'showIsotopicDistPlot'))
    showIsotopicDistPlot = options.showIsotopicDistPlot;
else
    showIsotopicDistPlot = false;
end

if(isfield(options,'methylation'))
    ismethyl = options.methylation;
else
    ismethyl = true;
    % error('options is not setup properly');
end

if(isfield(options,'ion'))
    ion = options.ion;
else
    ion = 'Na';    
end

% handle input argument
if(ismethyl)
    Hexoseunit   = mwconstant.HexMethyl;
    HexNACunit = mwconstant.HexNACMethyl;
    NeuAcunit   = mwconstant.NeuACMethyl;
    NeuGcunit = mwconstant.NeuGCMethyl;
    DeoHexunit =mwconstant.DeoHexMethyl;
    Pentoseunit = mwconstant.PentoseMethyl;
    Sulfunit =mwconstant.SulfMethyl;
    Phosunit  =mwconstant.PhosMethyl;
    KDNunit =mwconstant.KDNMethyl;
    HexAunit =mwconstant.HexAMethyl;
    redEndunit = mwconstant.redEndMethylForm;
else
    Hexoseunit   =mwconstant.Hex;
    HexNACunit =mwconstant.HexNAC;
    NeuAcunit   =mwconstant.NeuAC;
    NeuGcunit =mwconstant.NeuGC;
    DeoHexunit =mwconstant.DeoHex;
    Pentoseunit = mwconstant.Pentose;
    Sulfunit    =mwconstant.Sulf;
    Phosunit      =mwconstant.Phos;
    KDNunit     = mwconstant.KDN;
    HexAunit    = mwconstant.HexA;  
    redEndunit = mwconstant.redEndForm;
end

formulafieldnames = fieldnames(Hexoseunit);
for i=1:length(formulafieldnames)
    glycanForm.(formulafieldnames{i})= 0;    
end   

% add reducing end
glycanForm=formulaAdd(glycanForm,redEndunit);

%Hexose
nHex=findstr('h',gly);

glycanForm=formulaAdd(glycanForm,formulaMultiply(Hexoseunit,length(nHex)));

% HexNAC
nHexNAc=findstr('n',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(HexNACunit,length(nHexNAc)));

% Sailic Acid
nNeuAc=findstr('s',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(NeuAcunit,length(nNeuAc)));

% NEUGC
nNeuGc=findstr('g',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(NeuGcunit,length(nNeuGc)));

% fucose
nFuc=findstr('f',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(DeoHexunit,length(nFuc)));

%pentose
nXyl = findstr('x',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(Pentoseunit,length(nXyl)));

%sulfation
nSO3 = findstr('z',gly);                  % sulfation adds SO3 to MW
glycanForm=formulaAdd(glycanForm,formulaMultiply(Sulfunit,length(nSO3)));

%phosphoration
nPO3=findstr('p',gly);                  % phosphoration adds PO3H to MW
glycanForm=formulaAdd(glycanForm,formulaMultiply(Phosunit,length(nPO3)));

% KDN
nKDN=findstr('k',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(KDNunit,length(nKDN)));

%hexa
nHexA=findstr('u',gly);
glycanForm=formulaAdd(glycanForm,formulaMultiply(HexAunit,length(nHexA)));

%  add Na / Ka
if(strcmpi(ion,'none'))
    % do nothing
elseif(strcmpi(ion,'Na'))
    glycanForm.('Na')=1;
elseif(strcmpi(ion,'K'))
    glycanForm.('K')=1;
end

if(nargout==2)
    md=isotopicdist(glycanForm,'SHOWPLOT',showIsotopicDistPlot);
    varargout{1} =md;
elseif((nargout==1)||(nargout==0))
    if(showIsotopicDistPlot)
        isotopicdist(glycanForm);
    end
end

end

function aplusb=formulaAdd(a,b)
      formulafieldnames = fieldnames(a);
      if(isempty(b))
          aplusb = a;
          return
      end      
        for i=1:length(formulafieldnames)
              aplusb.(formulafieldnames{i})=...
                  a.(formulafieldnames{i})+....
                  b.(formulafieldnames{i});
        end
end
  
function amultiplyb=formulaMultiply(formula,numberStruct)
    if(numberStruct==0) 
        amultiplyb=[];
        return
    end    
    formulafieldnames = fieldnames(formula);
    for i=1:length(formulafieldnames)
         amultiplyb.(formulafieldnames{i}) = ...
            formula.(formulafieldnames{i}) * ....
            numberStruct;
    end
end
