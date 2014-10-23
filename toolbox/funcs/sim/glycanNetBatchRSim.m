function [glycanNetReactorModel,exitflag,speciesConc,speciesID,funceval] = glycanNetBatchRSim(...
               glyPathway,enzkineticsdb,substrComptOption,odesimoption)
%glycanNetReactorSim simulate the glycosylation reaction network in the
% specified reactor
%  
%  Syntax:  [tspan,specieids,speciesvalues] = glycanNetReactorSim(...
%                glyPathway,enzkineticsdb,reactortype,varargin) 
%
%  Input: 
%         glyPathway :  Pathway Object
%         enzKineticsdb : a Containers.Map object for storing enzyme kinetics 
%         reactortype: a string representing type of reactor assumed
%         flowpara (optional): a structure with 3 fields including qflow, kt, inletConcdb.          
% 
%  Output: tspan 
%
%
%See also glycanDSim, glyanSSim.

% Author: Gang Liu  
% Date Lastly Updated: 4/6/14

% set up enzyme kinetics
glyPathway = pathwaySetKinetics(glyPathway,enzkineticsdb,substrComptOption);

% set up CSTR reactor
glycanNetReactorModel = pathway2ReactorModel(glyPathway,'BatchR');

% perform steady-state analysis
ODEmfileName       =  'tempODE.m';
currentpath        =  pwd;
fullodefilename    =  fullfile(currentpath,ODEmfileName);
rmfileexpr = ['delete ', fullodefilename];

if(exist(fullodefilename,'file')==2)
   eval(rmfileexpr);
end

[exitflag,speciesConc,speciesID,funceval] =  glycanSSim(glycanNetReactorModel,fullodefilename,...
   odesimoption);

rmfileexpr = ['delete ', fullodefilename];
eval(rmfileexpr);







