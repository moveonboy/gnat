function glycanMW = glycanMolWt(gly,varargin)
%glycanMolWt calculate glycan mass
%
%    glycanMW = glycanMolWt(glycanstring)  takes an input of the glycan 
%       string and outputs its molecular weight. Default options for the 
%       calculation of molecular weight are methylated, Na ionized,
%       monoisotopic and free reducing end.
%
%   glycanMW = glycanMolWt(glycanstring,options)  uses customized options
%       structure for the calcuation of molecular weight. The fields in the
%       options include 'methylation' ([true]/false),'ion' ('none',['Na'],'K'), 
%       'mono' ([true]/false), 'redEnd' ([true]/false).
%   
%   Example1:
%        options.ion='none'; 
%        options.methylation=false;
%        glycanMolWt('shn',options); for 'NeuAcHexHexNAc'
%        ans=674.2381  -- this is the same values as GlycanMass
%        Linkage specificity mentioned in the text is ignored as are the brackets
%
%    Example 2:
%        glycanMolWt('sa2,3hb1,3{s2,6}n',options); for NeuAc2,3Hex2,3[NeuAc2,6]HexNAc
%        ans=965.3336
% 
%    Example 3:
%        glycanMolWt('h(s)n',options); Hex(NeuAc)HexNAc
%        ans=674.2381
%
%    Example 4:
%        slex='fshn'; forFucNeuAcHexHeN
%        glycanMolWt(strcat(slex,'hn'),options)
%        ans=1185.428231
%   
%   Example 5: 
%      glycanstring='hhhhhnn'; % Hex5HexNAc2
%      glycanMolWt(glycanstring);
%     
% See also mwconstant,glycanFormula.

%  Date: 8/2/12 3:32 PM Modified from Dr.Neelamegham's code
%  Copyright 2012 Neelamegham Lab

if(length(varargin)==1)
    options = varargin{1};
else
    options = struct('methylation',true,'ion','Na','mono',true,'redEnd',true);
end

if(isfield(options,'methylation'))
    ismethyl = options.methylation;
else
    ismethyl = true;
    % error('options is not setup properly');
end

if(isfield(options,'mono'))
    ismono = options.mono;
else
    ismono = true;
    % error('options is not setup properly');
end

if(isfield(options,'redEnd'))
    isredEnd = options.redEnd;
else
    isredEnd = true;
end

if(isfield(options,'ion'))
    ion = options.ion;
else
    ion = 'Na';
end

% handle input argument
% Note: % all masses are internal glycan fragments
if((ismethyl)&&(ismono))
    Hexosemwunit    = mwconstant.HexoseMonoMethyl;  
    HexNACmwunit  = mwconstant.HexNACMonoMethyl;  
    NeuAcmwunit    = mwconstant.NeuACMonoMethyl;  
    NeuGcmwunit    =  mwconstant.NeuGCMonoMethyl;
    DeoHexmwunit  = mwconstant.DeoHexMonoMethyl; 
    Pentosemwunit  = mwconstant.PentoseMonoMethyl;
    Sulfmwunit = mwconstant.SulfMonoMethyl;
    Phosmwunit = mwconstant.PhosMonoMethyl;
    KDNmwunit = mwconstant.KDNMonoMethyl;
    HexAmwunit  = mwconstant.HexAMonoMethyl;  
    
    
elseif((ismethyl)&&(~ismono))
    Hexosemwunit =mwconstant.HexoseAvgMethyl; 
    HexNACmwunit =mwconstant.HexNACAvgMethyl;  
    NeuAcmwunit =mwconstant.NeuACAvgMethyl;  
    NeuGcmwunit =mwconstant.NeuGCAvgMethyl;
    DeoHexmwunit =mwconstant.DeoHexAvgMethyl; 
    Pentosemwunit = mwconstant.PentoseAvgMethyl;
    Sulfmwunit  = mwconstant.SulfAvgMethyl;
    Phosmwunit  = mwconstant.PhosAvgMethyl;
    KDNmwunit  = mwconstant.KDNAvgMethyl;
    HexAmwunit =mwconstant.HexAAvgMethyl; 
    
elseif((~ismethyl)&&(ismono))
    Hexosemwunit =mwconstant.HexoseMonoUnde ;  
    HexNACmwunit =mwconstant.HexNACMonoUnde ; 
    NeuAcmwunit =mwconstant.NeuACMonoUnde ;  
    NeuGcmwunit =mwconstant.NeuGCMonoUnde ;
    DeoHexmwunit =mwconstant.DeoHexMonoUnde;  
    Pentosemwunit=mwconstant.PentoseMonoUnde ;
    Sulfmwunit=mwconstant.SulfMonoUnde ;
    Phosmwunit=mwconstant.PhosMonoUnde ;
    KDNmwunit=mwconstant.KDNMonoUnde ;
    HexAmwunit =mwconstant.HexAMonoUnde ;  
    
elseif((~ismethyl)&&(~ismono))
    Hexosemwunit =mwconstant.HexoseAvgUnde;  
    HexNACmwunit =mwconstant.HexNACAvgUnde; 
    NeuAcmwunit =mwconstant.NeuACAvgUnde; 
    NeuGcmwunit =mwconstant.NeuGCAvgUnde;
    DeoHexmwunit =mwconstant.DeoHexAvgUnde; 
    Pentosemwunit =mwconstant.PentoseAvgUnde;
    Sulfmwunit=mwconstant.SulfAvgUnde;
    Phosmwunit=mwconstant.PhosAvgUnde;
    KDNmwunit=mwconstant.KDNAvgUnde;
    HexAmwunit =mwconstant.HexAAvgUnde;  
end

if(isredEnd)
    if(ismethyl)
       glycanMW=mwconstant.redEndMethyl;
    else
        glycanMW=mwconstant.redEnd;
    end
else
    glycanMW=0;
end

% if ion options
if(strcmpi(ion,'none'))
    % do nothing
elseif(strcmpi(ion,'Na'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.NaionMono;
    else
      glycanMW= glycanMW+mwconstant.NaionAvg;
    end    
elseif(strcmpi(ion,'K'));    
    if(ismono)
       glycanMW=glycanMW+mwconstant.KionMono;
    else
      glycanMW= glycanMW+mwconstant.KionAvg;
    end        
end

%Hexose
nHex=findstr('h',gly);
glycanMW=glycanMW+length(nHex)*Hexosemwunit;

% HexNAC
nHexNAc=findstr('n',gly);
glycanMW=glycanMW+length(nHexNAc)*HexNACmwunit;

% Sailic Acid
nNeuAc=findstr('s',gly);
glycanMW=glycanMW+length(nNeuAc)*NeuAcmwunit;

%NEUGC
nNeuGc=findstr('g',gly);
glycanMW=glycanMW+length(nNeuGc)*NeuGcmwunit;

%DeoxyHex
nFuc=findstr('f',gly);
glycanMW=glycanMW+length(nFuc)*DeoHexmwunit;

%Xylose
nXyl=findstr('x',gly);
glycanMW=glycanMW+length(nXyl)*Pentosemwunit;

%Sulfation
nSO3=findstr('z',gly);                  % sulfation adds SO3 to MW
glycanMW=glycanMW+length(nSO3)*Sulfmwunit;

%Phosphorylation
nPO3=findstr('p',gly);                  % phosphoration adds PO3H to MW
glycanMW=glycanMW+length(nPO3)*Phosmwunit;

%KDN
nKDN=findstr('k',gly);
glycanMW=glycanMW+length(nKDN)*KDNmwunit;

%HEXA
nHexA=findstr('u',gly);
glycanMW=glycanMW+length(nHexA)*HexAmwunit;

end