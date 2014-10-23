function [tspan,specieids,speciesvalues] = glycanNetPFRSim(...
                    glyPathway,enzkineticsdb,varargin)
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

% output ODE function for glycanNetModel object


