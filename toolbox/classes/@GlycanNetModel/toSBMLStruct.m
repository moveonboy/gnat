function sbmlstruct = toSBMLStruct(obj,varargin)
%toSBMLStruct returns a SBML structure
%   
%  Example:   
%       glycanNetExampleFileName = 'gnat_test_wtannot.xml';
%       glycanNetModelObj        =  glycanNetSBMLread(...
%           glycanNetExampleFileName);
%       sbmlstructure            = glycanNetModelObj.toSBMLStruct;
%    
% See also GLYCANNETMODEL.

% Author: Gang Liu
% Date Lastly Updated: 4/7/14

% default SBML level and version L2V1
if(isempty(varargin))
    defaultlevel    = 2;
    defaultver      = 1;
    addGlycanstruct = false;
elseif(length(varargin)==1)
    addGlycanstruct = varargin{1};
    defaultlevel    = 2;
    defaultver      = 1;
elseif(length(varargin)==3)
    addGlycanstruct = varargin{1};
    defaultlevel    = varargin{2};
    defaultver      = varargin{3};
end

sbmlstruct  = Model_create(defaultlevel,defaultver);

% add compartment info to SBML 
if(isa(obj.getCompartment,'CellArrayList'))
    numCompt = length(obj.getCompartment);
    for i=1 : numCompt
        % set up each compartment
       ithcompt    = obj.getCompartment.get(i);
       sbmlcompt   = setSBMLCompt(ithcompt,defaultlevel,defaultver,i); 
       sbmlstruct  = Model_addCompartment(sbmlstruct, sbmlcompt);   
    end
elseif(isa(obj.getCompartment,'Compt'))
       ithcompt    = obj.getCompartment;
       sbmlcompt   = setSBMLCompt(ithcompt,defaultlevel,defaultver,1); 
       sbmlstruct  = Model_addCompartment(sbmlstruct, sbmlcompt);   
end

% add species info to SBML
numSpecies = obj.getGlycanPathway.getNSpecies;
for i = 1 : numSpecies
    % set up each species
    ithspecies     = obj.getGlycanPathway.getSpecies.get(i);
    sbmlspecies    = setSBMLSpecies(ithspecies,defaultlevel,defaultver,i,addGlycanstruct);
    sbmlstruct     = Model_addSpecies(sbmlstruct,sbmlspecies);   
end

% set up global parameter defined in the enzyme to SBML 
numEnzs = obj.getGlycanPathway.theEnzs.length;
for i = 1 : numEnzs
    ithEnz = obj.getGlycanPathway.theEnzs.get(i);
    sbmlstruct = setSBMLEnzGlobalParameter(sbmlstruct,ithEnz,...
        defaultlevel,defaultver);
end

% set up transport kinetic parameter
kt = obj.compartment.get(1).reactortype.ktransport;
qin = obj.compartment.get(1).reactortype.qflowin;
sbmlstruct = setSBMLTransportParameter(sbmlstruct,kt,qin,defaultlevel,defaultver);

% add reaction info to SBML
numRxns = obj.getGlycanPathway.getNReactions;
for i = 1 : numRxns
    % set up each species
    ithRxn          = obj.getGlycanPathway.getReactions.get(i);
    sbmlrxn         = setSBMLRxn(ithRxn,defaultlevel,defaultver,i);    
    sbmlstruct      = Model_addReaction(sbmlstruct, sbmlrxn);
end

end


% set up sbml rxn
function sbmlrxn = setSBMLRxn(ithrxn,level,ver,i)
    sbmlrxn = Reaction_create(level,ver);
    
    % set reaction name
    if(~isempty(ithrxn.name))
        sbmlrxn = Reaction_setName(sbmlrxn,ithrxn.name);
    else
        ithrxnname = sprintf('rxn%d',i);
        sbmlrxn = Reaction_setName(sbmlrxn,ithrxnname);
    end
    
    % set reaction id
    if(~isempty(ithrxn.id))
        sbmlrxn = Compartment_setId(sbmlrxn,ithrxn.id);
    else
        ithrxnid =  sprintf('rxn%d',i);
        sbmlrxn = Reaction_setId(sbmlrxn,ithrxnid);
    end
    
    % set reactant reference
    if(~isempty(ithrxn.reac))
        sbmlreactant               = SpeciesReference_create(level, ver);
        if(~isempty(ithrxn.reac.id))
           reacspeciesid           = ithrxn.reac.id;
        else
           reacspeciesid           = sprintf('species%d',i);
        end
        sbmlreactant.species  = reacspeciesid;
        sbmlrxn               = Reaction_addReactant(sbmlrxn, sbmlreactant);
    end
     
    % set product reference
    if(~isempty(ithrxn.prod))
        sbmlproduct = SpeciesReference_create(level, ver);
        if(~isempty(ithrxn.prod.id))
           prodspeciesid  = ithrxn.prod.id;
        else
           prodspeciesid  = sprintf('species%d',i);
        end
        sbmlproduct.species  = prodspeciesid;
        sbmlrxn              = Reaction_addProduct(sbmlrxn, sbmlproduct);
    end
    
    % set kinetic laws
    if ~isempty(ithrxn.rxnkinetics)
        sbmlkineticlaw = setSBMLkineticlaw(ithrxn,level,ver); 
        sbmlrxn        = Reaction_setKineticLaw(sbmlrxn,sbmlkineticlaw);    
    end
end

function sbmlstruct = setSBMLEnzGlobalParameter(sbmlstruct,ithEnz,level,ver)
    if(strcmpi(ithEnz.enzkinetics.mechanism,...
            TypeStringConst.bibiSeqOrderwSubstInhib))
        % add kinetic parameters
         parameternames  = {['Vm' ithEnz.id],['KmG' ithEnz.id],...
                          ['KmS' ithEnz.id],['Ke' ithEnz.id]};
         parametervalues = [ithEnz.enzkinetics.sumconsts.Vm,...
             ithEnz.enzkinetics.sumconsts.Kmd_G1E,...
             ithEnz.enzkinetics.sumconsts.Kmd_SNE,...
             ithEnz.enzkinetics.sumconsts.Keq];
         numparameter    = length(parameternames);
         for i = 1 : numparameter
            % create a kineticLaw sbml structure
            sbmlparameter  = Parameter_create(level, ver);
            sbmlparameter  = Parameter_setConstant(sbmlparameter, 1);
            sbmlparameter  = Parameter_setId(sbmlparameter,parameternames{i});        
            sbmlparameter  = Parameter_setName(sbmlparameter, parameternames{i});
            sbmlparameter  = Parameter_setValue(sbmlparameter,...
                   parametervalues(i));
            sbmlstruct     = Model_addParameter(sbmlstruct, sbmlparameter);
         end
    elseif(strcmpi(ithEnz.enzkinetics.mechanism,TypeStringConst.SMM) ||...
           strcmpi(ithEnz.enzkinetics.mechanism,TypeStringConst.Comp))
          % add kinetic parameters
         parameternames  = {['Vm' ithEnz.id],['Km' ithEnz.id]};
         parametervalues = [ithEnz.enzkinetics.sumconsts.Vm,...
             ithEnz.enzkinetics.sumconsts.Km];
         numparameter    = length(parameternames);
         for i = 1 : numparameter
            % create a kineticLaw sbml structure
            sbmlparameter  = Parameter_create(level, ver);
            sbmlparameter  = Parameter_setConstant(sbmlparameter, 1);
            sbmlparameter  = Parameter_setId(sbmlparameter,parameternames{i});        
            sbmlparameter  = Parameter_setName(sbmlparameter, parameternames{i});
            sbmlparameter  = Parameter_setValue(sbmlparameter,...
                   parametervalues(i));
            sbmlstruct     = Model_addParameter(sbmlstruct, sbmlparameter);
         end
    else
        error('MATLAB:GNAT:NOTSUPPORTED','TO BE WRITTEN');
    end
end


function sbmlstruct = setSBMLTransportParameter(sbmlstruct,kt,qin,level,ver)
         parameternames  = {'kt','qin'};
         parametervalues = [kt,qin];
         numparameter    = length(parameternames);
         for i = 1 : numparameter
            % create a kineticLaw sbml structure
            sbmlparameter  = Parameter_create(level, ver);
            sbmlparameter  = Parameter_setConstant(sbmlparameter, 1);
            sbmlparameter  = Parameter_setId(sbmlparameter,parameternames{i});        
            sbmlparameter  = Parameter_setName(sbmlparameter, parameternames{i});
            sbmlparameter  = Parameter_setValue(sbmlparameter,...
                   parametervalues(i));
            sbmlstruct     = Model_addParameter(sbmlstruct, sbmlparameter);
         end
end


% set up sbml kinetic laws
function sbmlkineticlaw = setSBMLkineticlaw(ithrxn,level,ver)
 sbmlkineticlaw = KineticLaw_create(level,ver);
 rxnkinsformula = ithrxn.rxnkinetics.kineticlaws.mathformula; % add formula
 sbmlkineticlaw = KineticLaw_setFormula(sbmlkineticlaw, rxnkinsformula);
end

% set up sbml species
function sbmlspecies = setSBMLSpecies(ithspecies,defaultlevel,defaultver,i,addglycanstruct)
    sbmlspecies = Species_create(defaultlevel,defaultver);
    
    % set species name
    if(~isempty(ithspecies.name))
      sbmlspecies = Species_setName(sbmlspecies,ithspecies.name);
    else
      newname =  ['species' num2str(i)];
      sbmlspecies = Species_setName(sbmlspecies,newname);
    end
    
     % set species initconc
    if(~isempty(ithspecies.initConc))
        sbmlspecies = Species_setInitialConcentration(sbmlspecies,ithspecies.initConc);
    else
        sbmlspecies = Species_setInitialConcentration(sbmlspecies,0);
    end
    
    % set species initAmount
    if(~isempty(ithspecies.initAmount))
        sbmlspecies = Species_setInitialAmount(sbmlspecies,ithspecies.initAmount);
    else
        sbmlspecies = Species_setInitialAmount(sbmlspecies,0);
    end
    
    % set species id
    if(~isempty(ithspecies.id))
        sbmlspecies = Compartment_setId(sbmlspecies,ithspecies.id);
    else
        newid =  ['s' num2str(i)];
        sbmlspecies = Species_setId(sbmlspecies,newid);
    end
    
    % set compartment name
    if(~isempty(ithspecies.compartment))
        if(~isempty(ithspecies.compartment.id))
            comptid     = ithspecies.compartment.id;
            sbmlspecies = Species_setCompartment(sbmlspecies,comptid);
        else
            comptname   = ithspecies.compartment.name;
            sbmlspecies = Species_setCompartment(sbmlspecies,comptname);
        end
    end
    
    % set structure annotation for sbml
    if(addglycanstruct)
        glycostr        = glycanStrwrite(ithspecies.glycanStruct);
        annotationstr   = setAnnotationStr(glycostr);
        sbmlspecies.annotation =  annotationstr;  
    end
end

function annotationstr=setAnnotationStr(glycostr)
    annotationstr =['<annotation>' ...
                     '<glycoct xmlns="http://www.eurocarbdb.org/recommendations/encoding">'];
    annotationstr =[annotationstr glycostr '</glycoct>' '</annotation>'];
    annotationstr =strrep(annotationstr,'<?xml version="1.0" encoding="UTF-8"?>', ' ');
end

% set up sbml compartment
function sbmlcompt = setSBMLCompt(ithcompt,defaultlevel,defaultver,i)
    sbmlcompt = Compartment_create(defaultlevel,defaultver);
    
    % set compartment name
     if(~isempty(ithcompt.name))
        sbmlcompt = Compartment_setName(sbmlcompt,ithcompt.name);
     else
        newcomptname =  ['compt' num2str(i)];
        sbmlcompt = Compartment_setName(sbmlcompt,newcomptname);
     end    
   
    % set compartment id
    if(~isempty(ithcompt.id))
        sbmlcompt = Compartment_setId(sbmlcompt,ithcompt.id);    
    elseif(~isempty(ithcompt.name))
        sbmlcompt = Compartment_setId(sbmlcompt,ithcompt.name);
    else
        newcomptid =  ['c' num2str(i)];
        sbmlcompt = Compartment_setId(sbmlcompt,newcomptid);        
    end
    
    % set compartment units
    if(~isempty(ithcompt.units))
        sbmlcompt = Compartment_setUnits(sbmlcompt,ithcompt.units);
    end   
    
    % set compartment size
    if(~isempty(ithcompt.size))
        sbmlcompt = Compartment_setSize(sbmlcompt,ithcompt.size);
    end
    
    % set compartment constant
    sbmlcompt = Compartment_setConstant(sbmlcompt,1);
     
    % set compartment spatial dimensions
    if(~isempty(ithcompt.spatialDimensions))
        sbmlcompt = Compartment_setSpatialDimensions(sbmlcompt,...
            (ithcompt.spatialDimensions));
    end 
    
end