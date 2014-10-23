function glycanNetWriteODEFunction(varargin)
% glycanNetWriteODEFunction(SBMLModel, name(optional)) adpated from
% glycanNetWriteODEFunction 
%
% Takes 
% 
% 1. SBMLModel, an SBML Model structure
% 2. name, an optional string representing the name of the ode function to be used
% 
% Outputs 
%
% 1. a file 'name.m' defining a function that defines the ode equations of
%   the model for use with the ode solvers
%    (if no name supplied the model id will be used)


%<!---------------------------------------------------------------------------
% This file is part of SBMLToolbox.  Please visit http://sbml.org for more
% information about SBML, and the latest version of SBMLToolbox.
%
% Copyright (C) 2009-2012 jointly by the following organizations: 
%     1. California Institute of Technology, Pasadena, CA, USA
%     2. EMBL European Bioinformatics Institute (EBML-EBI), Hinxton, UK
%
% Copyright (C) 2006-2008 jointly by the following organizations: 
%     1. California Institute of Technology, Pasadena, CA, USA
%     2. University of Hertfordshire, Hatfield, UK
%
% Copyright (C) 2003-2005 jointly by the following organizations: 
%     1. California Institute of Technology, Pasadena, CA, USA 
%     2. Japan Science and Technology Agency, Japan
%     3. University of Hertfordshire, Hatfield, UK
%
% SBMLToolbox is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation.  A copy of the license agreement is provided
% in the file named "LICENSE.txt" included with this software distribution.
%----------------------------------------------------------------------- -->


switch (nargin)
    case 0
        error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'must have at least one argument');
    case 1
        SBMLModel = varargin{1};
        filename = 'defaultSBMLODEName';
    case 2
        SBMLModel = varargin{1};
        filename  = varargin{2};
    otherwise
        error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'does not take more than two arguments');
end;

% check input is an SBML model
if (~isValidSBML_Model(SBMLModel))
    error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'first argument must be an SBMLModel structure');
end;

% -------------------------------------------------------------
% check that we can deal with the model
% for i=1:length(SBMLModel.parameter)
%   if (SBMLModel.parameter(i).constant == 0)
%     error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with varying parameters');
%   end;
% end;
if SBMLModel.SBML_level > 2
  if ~isempty(SBMLModel.conversionFactor)
    error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with conversion factors');
  end;
  for i=1:length(SBMLModel.species)
    if ~isempty(SBMLModel.species(i).conversionFactor)
      error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with conversion factors');
    end;
  end;
end;

for i=1:length(SBMLModel.compartment)
  if (SBMLModel.compartment(i).constant == 0)
    error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with varying compartments');
  end;
end;
if (isempty(SBMLModel.species))
     error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with models with no species');
end;  
for i=1:length(SBMLModel.event)
  if (~isempty(SBMLModel.event(i).delay))
    error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with delayed events');
  end;
  if SBMLModel.SBML_level > 2
    if (~isempty(SBMLModel.event(i).priority))
      error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with event priorities');
    end;
    if (~isempty(SBMLModel.event(i).trigger) &&  ...
        (SBMLModel.event(i).trigger.initialValue == 1 || SBMLModel.event(i).trigger.persistent == 1))
      error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with persistent trigger');
    end;
  end;
end;
for i=1:length(SBMLModel.reaction)
  if (SBMLModel.reaction(i).fast == 1)
    error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with fast reactions');
  end;
end;
if (length(SBMLModel.compartment) > 1)
  error('WriteODEFunction(SBMLModel, (optional) filename)\n%s', 'Cannot deal with multiple compartments');
end;
if (SBMLModel.SBML_level > 1 && ~isempty(SBMLModel.time_symbol))
  for i=1:length(SBMLModel.rule)
    if (strcmp(SBMLModel.rule(i).typecode, 'SBML_ASSIGNMENT_RULE'))
      if (~isempty(matchName(SBMLModel.rule(i).formula, SBMLModel.time_symbol)))
        error('Cannot deal with time in an assignment rule');
      end;
    end;
  end;
end;
if (SBMLModel.SBML_level > 1 && ~isempty(SBMLModel.delay_symbol))
  for i=1:length(SBMLModel.rule)
    if (strcmp(SBMLModel.rule(i).typecode,'SBML_ASSIGNMENT_RULE'))
      if (~isempty(matchName(SBMLModel.rule(i).formula, SBMLModel.delay_symbol)))
        error('Cannot deal with delay in an assignment rule');
      end;
    end;
  end;
end;

%--------------------------------------------------------------
% get information from the model
% [stochimatrix, species] = GetStoichiometryMatrix(SBMLModel);
[ParameterNames, ParameterValues] = GetAllParameters(SBMLModel);
% [VarParams, VarInitValues] = GetVaryingParameters(SBMLModel);
% NumberParams = length(VarParams);
NumberSpecies = length(SBMLModel.species);

if NumberSpecies > 0
  % Species = AnalyseSpecies(SBMLModel);
  [Speciesnames,~] = GetSpecies(SBMLModel);
end;
% if NumberParams > 0
%   Parameters = AnalyseVaryingParameters(SBMLModel);
% end;
if ~isempty(SBMLModel.compartment)
  [CompartmentNames, CompartmentValues] = GetCompartments(SBMLModel);
else
  CompartmentNames = [];
end;

% if (NumberParams + NumberSpecies) == 0
%   error('Cannot detect any variables');
% end;

if (SBMLModel.SBML_level > 1)
    NumEvents = length(SBMLModel.event);
    NumFuncs  = length(SBMLModel.functionDefinition);

    % version 2.0.2 adds the time_symbol field to the model structure
    % need to check that it exists
    if (isfield(SBMLModel, 'time_symbol'))
        if (~isempty(SBMLModel.time_symbol))
            timeVariable = SBMLModel.time_symbol;
        else
            timeVariable = 'time';
        end;
    else
        timeVariable = 'time';
    end;
    if ((SBMLModel.SBML_level == 2 &&SBMLModel.SBML_version > 1) || ...
        (SBMLModel.SBML_level > 2))
      if (~isempty(SBMLModel.constraint))
        error('Cannot deal with constraints.');
      end;
    end;

else
    NumEvents = 0;
    timeVariable = 'time';
end;

%---------------------------------------------------------------
% get the name/id of the model

Name = '';
if (SBMLModel.SBML_level == 1)
    Name = SBMLModel.name;
else
    if (isempty(SBMLModel.id))
        Name = SBMLModel.name;
    else
        Name = SBMLModel.id;
    end;
end;

if (~isempty(filename))
    Name = filename;
elseif (length(Name) > 63)
    Name = Name(1:60);
end;

fileName = strcat(Name, '.m');
%--------------------------------------------------------------------
% open the file for writing

fileID = fopen(fileName, 'w');

% write the function declaration
% if no events and using octave
fprintf(fileID,  'function xd = %s(%s,x)\n',Name,timeVariable);

% need to add comments to output file
fprintf(fileID, '%%\t[t,x] = ode15s(@%s, [0, t_end], %s)\n', Name, Name);

% write the parameter values

for i = 1:length(ParameterNames)
    fprintf(fileID, '\t%s = %g;\n', ParameterNames{i}, ParameterValues(i));
end;

% write the initial concentration values for the species
for i = 1 : NumberSpecies
    fprintf(fileID, '\t%s = x(%u);\n', Speciesnames{i}, i);
end;

% write code to calculate concentration values
% print rate laws for each reaction 
numRxns = length(SBMLModel.reaction);
fprintf(fileID,'\tV = zeros(%d,1);\n',numRxns);

% output rate laws for each reaction
for j = 1: numRxns
   oldformula = SBMLModel.reaction(j).kineticLaw.formula;
   newformula = regexprep(oldformula,'*1(?=[^0-9])','');
   newformula = regexprep(newformula,'\1(?=[^0-9])','');
   newformula = regexprep(newformula,'/1(?=[^0-9])','');
   fprintf(fileID,'\tV(%d) = %s;\n',j,newformula);
end

% output equation for xd=SM*V
fprintf(fileID,'\txd=SM''*V;\n');
fprintf(fileID, 'end\n');
fclose(fileID);
end