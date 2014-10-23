function glycanNetReactorModel = pathway2ReactorModel(glycanpathway,reactortype,varargin)
%pathway2ReactorModel create a reactor model based on reactor type    
%   provided enzyme kinetics is defined in the model.
%
%  Input: reactortype, reactorpara (optional)
%      1. reactortype, the type of reactor by which the pathway is modelled
%
%      2. reactorpara, a MATLAB structure containing 3 fields: qflow (scalar), 
%           kt (scalar), inlet conc(a CellArrayList with elemetents bearing
%           two fields including "species" and "inletconc". 
%  
%  Output: glycanNetReactorModel
%       glycanNetReactorModel, a GlycanNetModel object.
%      
%See also GlycanNetModel.

% Author: Gang Liu
% Date Lastly Updated: 4/8/14

narginchk(2,3);

% set compartment reactor type
if(strcmpi(reactortype,'CSTR'))
  glycanpathway.setComptReactorType(TypeStringConst.CSTR);
elseif(strcmpi(reactortype,'BatchR'))
  glycanpathway.setComptReactorType(TypeStringConst.BatchR);
elseif(strcmpi(reactortype,'PFR'))
  glycanpathway.setComptReactorType(TypeStringConst.PFR);
else
  error('MATLAB:GNAT:ERRORINPUT','THE REACTOR TYPE IS NOT SUPPORTED');
end 
  
% set qflow and k in each compt
if(strcmpi(reactortype,'PFR')) ...
  ||(strcmpi(reactortype,'CSTR')) 
 if(length(varargin)==1)
   qflow  = varargin{1}.qflow;
   k      = varargin{1}.kt;
   speciesInletConc = varargin{1}.inletConc;
   
   % speciesInletConc must be a CellArrayList object
   if(~isa(speciesInletConc,'CellArrayList'))
       error('MATLAB:GNAT:ERRORINPUT',...
           'Inlet Concentration must be a CellArrayList object');
   end
   
   glycanpathway.setComptsQflow(qflow);
   glycanpathway.setComptskt(k);
    
   %find first compt and set inlet concentrations
   firstcompt = glycanpathway.identifyfirstcompt;
   firstcompt.reactortype.setSpeciesInletDB(speciesInletConc)   
  
   %set up transport kinetics
   glycanpathway.setTransportRxnKinetics;
  else
    error('MATLAB:GNAT:ERRORINPUT','THE INPUT FORMAT IS NOT CORRECT');  
 end
end

glycanNetReactorModel = GlycanNetModel(glycanpathway.compartment,...
                        glycanpathway,glycanpathway.name);
end

