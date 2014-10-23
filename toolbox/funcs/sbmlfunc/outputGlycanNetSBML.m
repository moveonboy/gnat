function  outputGlycanNetSBML(glycanNetModelobj,varargin)
%OUTPUTGLYCANNETSBML output SBML file to a local file 
%  Input: 
%    glycanNetModelobj - a GlycanNetModel object
%    fullfilename (optional) - SBML file storing glycosylation reaction networks
%  
%  Output:
%    an SBML file will be saved in the local directory. 
%
%See also GlycanNetModel.

% Author: Gang Liu
% Date Lastly Updated: 4/9/14

narginchk(1,2);

if(length(varargin)==1)
   if(~ischar(varargin{1}))
       error('MATLAB:GNAT:ERRORINPUT','WRONG TYPE OF INPUTS');    
   end
end

if(isempty(glycanNetModelobj.glycanNet_sbmlmodel))
    glycanNetModelobj.glycanNet_sbmlmodel = glycanNetModelobj.toSBMLStruct();
end

if(length(varargin)==1)
  outputSBML(glycanNetModelobj.glycanNet_sbmlmodel,varargin{1});
elseif(isempty(varargin))
  outputSBML(glycanNetModelobj.glycanNet_sbmlmodel);  
end

end

