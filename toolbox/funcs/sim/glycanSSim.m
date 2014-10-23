function [exitflag, speciesConc,speciesID,funceval]=glycanSSim(glycanNMObj,fullODEmfileName,varargin)
%glycanSSim simulate glycosylation reaction network model under steady state 
%
% [exitflag, speciesConc,speciesName,funceval]=glycanSSim(glycanNMObj,ODEmfileName,options) 
%  takes GlycanNetModel object, ODEFunction file to simulate the system under steady states using 
%  fsolve function provided by MATLAB. Options can be provided using the optimset command. See 
%  optimset and fsolve for details. 
%
% glycanSSim(glycanNMObj,ODEmfileName) uses default options. The iterative steps will show in the 
%  command windows,and the algorithm used to solve equations is levenberg-marquardt methods. 
%  See optimset for other default options. 
%
% Example : 
%  testGlycanNetSBMLfileName ='gnat_test_ssim.xml';
%  testGlycanNetModel= glycanNetSBMLread(testGlycanNetSBMLfileName);
%  % write out ODE function for simulation
%  if(isempty(testGlycanNetModel.glycanNet_sbmlmodel.name))
%   exampleODEFileName = testGlycanNetModel.glycanNet_sbmlmodel.id ;
%  else
%   exampleODEFileName = testGlycanNetModel.glycanNet_sbmlmodel.name ;
%  end 
%  if(~(exist(strcat(exampleODEFileName,'.m'),'file')==2))
%   WriteODEFunction(testGlycanNetModel.glycanNet_sbmlmodel,exampleODEFileName);
%  end
%  [exitflag, speciesconc, speciesnames,funceval] = glycanSSim(testGlycanNetModel,exampleODEFileName);
%    %To see simulation result in plot, type:
%  bar(speciesconc);
%  set(gca,'Xtick',1:4,'XTickLabel',speciesnames); 
%
% See also GlycanNetModel,glycanDSim.

% Author: Gang Liu
% Copyright 2012 Neelamegham Lab

if(isempty(glycanNMObj.glycanNet_sbmlmodel))
    glycanNet_SBMLStruct = glycanNMObj.toSBMLStruct;
else
    glycanNet_SBMLStruct = glycanNMObj.glycanNet_sbmlmodel;
end

% retrieve species name for the list
nSpecies =  length(glycanNet_SBMLStruct.species);
speciesID = cell(nSpecies,1);
for i=1: nSpecies
   speciesID{i,1} = glycanNet_SBMLStruct.species(1,i).id;
end

% write M file for ODE simulation
% ODEmfileName = [ODEmfileName '.m'];
if(exist(fullODEmfileName,'file')~=2)
    outputGlycanNetModelODE(glycanNMObj,fullODEmfileName);   
end

% create an anoymous function for fsolve 
[path,ODEfuncname,ext]=fileparts(fullODEmfileName);
% ODEmfileHandle = str2func(strrep(fullODEmfileName,'.m',''));
ODEmfileHandle = str2func(ODEfuncname);
fsolvefun = @(x) ODEmfileHandle(1,x);
x0 = ODEmfileHandle();

if(length(varargin)==1)
    options = varargin{1};
    [speciesConc,funceval,exitflag]= fsolve(fsolvefun,x0,options);
else
    options=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.01}); 
     [speciesConc,funceval,exitflag]=fsolve(fsolvefun,x0,options);
end

end

