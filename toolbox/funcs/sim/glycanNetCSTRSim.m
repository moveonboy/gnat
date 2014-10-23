function [glycanNetReactorModel,exitflag,speciesID,speciesConc,funceval] = ... 
   glycanNetCSTRSim(glyPathway,enzkineticsdb,substrComptOption,reactorpara,odesimoption)
%glycanNetCSTRSim simulate the glycosylation reaction network in the
%   specified reactor
%  
%  Syntax:  [tspan,specieids,speciesvalues] = glycanNetReactorSim(...
%                glyPathway,enzkineticsdb,substrComptOption,reactorpara) 
%
%  Input: 
%         glyPathway:    a Pathway Object
%         enzKineticsdb: a Containers.Map object for storing enzyme kinetics 
%         substrComptOption: option for substrate competition 
%         reactorpara (optional): a structure with 3 fields including qflow, kt, inletConcdb. 
%         odesimoption: option for ODE simulation
% 
%  Output:  
%         1) exitflag
%         2) speciesID
%         3) speciesConc 
%         4) funceval
%
%  Syntax: [] = glycanNetReactorSim()
%
%
%See also glycanDSim, glyanSSim.

% Author: Gang Liu  
% Date Lastly Updated: 4/6/14

% set up enzyme kinetics
glyPathway = pathwaySetKinetics(glyPathway,enzkineticsdb,substrComptOption);

% set up CSTR reactor
glycanNetReactorModel = pathway2ReactorModel(glyPathway,'CSTR',reactorpara);

% perform steady-state analysis
ODEmfileName                   = 'tempODE.m';
currentpath                    = pwd;
fullodefilename                = fullfile(currentpath,ODEmfileName);
rmfileexpr                     = ['delete ', fullodefilename];

if(exist(fullodefilename,'file')==2)
   eval(rmfileexpr);
end

outputGlycanNetModelODE(glycanNetReactorModel,fullodefilename);
[exitflag, speciesConc, speciesID, funceval] =  glycanSSim(...
    glycanNetReactorModel,fullodefilename,odesimoption);

rmfileexpr = ['delete ', fullodefilename];
eval(rmfileexpr);

end

