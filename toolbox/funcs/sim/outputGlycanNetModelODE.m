function outputGlycanNetModelODE(glycanNetModelObj,fullmfilename)
%outputGlycanNetModelODE    
%
%
% See also outputGlycanNetSBML.

% Author: Gang Liu
% Date Lastly Updated: 4/8/14

narginchk(2,2);
if(isempty(glycanNetModelObj.glycanNet_sbmlmodel))
    glycanNetModelObj.glycanNet_sbmlmodel = glycanNetModelObj.toSBMLStruct;
end

if(length(glycanNetModelObj.compartment)>1)
   singleComptSBMLModel = ...
       toSingleCompt(glycanNetModelObj.glycanNet_sbmlmodel);
elseif(glycanNetModelObj.compartment.length==1)
   singleComptSBMLModel =  glycanNetModelObj.glycanNet_sbmlmodel;
else
   error('MATLAB:GNAT:NOCOMPTDEFINED','No Compartments are Defined');
end

[directorytowrite,filename,fileext]=fileparts(fullmfilename);
currpwd = pwd;
eval(['cd ' directorytowrite]);
mfilename = [filename];
WriteODEFunction(singleComptSBMLModel,mfilename);
eval(['cd ' currpwd]);
end

function sbml_structure_singlecompt = toSingleCompt(sbml_structure_multicompt)
% create a global compartment
[level, version] = GetLevelVersion(sbml_structure_multicompt);
newComptName ='Global';
nCompts =  length(sbml_structure_multicompt.compartment)       ;
for i=1:nCompts
     newComptName =[newComptName '_' Compartment_getName(sbml_structure_multicompt.compartment(1,i))];    
end

maxlengthstrvar = 40;
comptGlobal = Compartment_create(level,version);
if(length(newComptName)>=maxlengthstrvar)
    comptGlobal = Compartment_setId(comptGlobal,'Global_compts');
else
    comptGlobal = Compartment_setId(comptGlobal,newComptName);
end
comptGlobal = Compartment_setName(comptGlobal,...
    [newComptName '_' num2str(nCompts) 'compts']);

comptGlobal.constant=1;
comptGlobal.size=1;

sbml_structure_multicompt.compartment = comptGlobal;
sbml_structure_multicompt.id = sbml_structure_multicompt.name;

for i=1:length(sbml_structure_multicompt.species)
    sbml_structure_multicompt.species(1,i)= Species_setCompartment(sbml_structure_multicompt.species(1,i),...
        Compartment_getId(comptGlobal));
end

sbml_structure_singlecompt=sbml_structure_multicompt;
end