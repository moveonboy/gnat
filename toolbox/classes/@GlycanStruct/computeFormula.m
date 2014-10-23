function glycanformula = computeFormula(obj,varargin)
%computeFormula compute the chemical formula of a glycan object
%  computeFormula(GLYCANSTRUCTobj) computes the chemical
%   formula
%
% See also GlycanStruct

if(length(varargin)==1)
    options = varargin{1};
else
    options = struct('methylation',true,'ion','Na');
end
glycanCompString       = glyCompString(obj);
gly1letstring          = gly1charformat(glycanCompString);
glycanformula          = glycanFormula(gly1letstring,options);
end

