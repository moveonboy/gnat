classdef MMenRxnKinetics < RxnKinetics
    %MMENRXNKINETICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        enzkinetics;
        substratekinetics      = struct('Vm',[],'Km',[],'speciesid',[]);
        substratekineticsratio = struct('Vm_ratio',1,'Km_ratio',1,'speciesid',1);
        inhibitorkinetics %=struct('Km',[],'Km_ratio',[],'speciesid','');
        kineticlaws
        %  a structure with 4 fields: mathformula, parameternames
        %  parametervalues, speciesid
    end
    
    methods  % Constructor for MMenRxnKinetics define substrate specific rate constants
        function obj=MMenRxnKinetics(varargin)
            if(isempty(varargin))
                return;
            end
            obj.enzkinetics  = varargin{1};
            glycanspecies    = varargin{2};
            
            obj.substratekinetics.speciesid = glycanspecies.id;
                        
            % derive substrate-specific rate constant. Apply specificity rule if necessary
            if  isempty(obj.enzkinetics.substrspecdb)...
                    && (obj.enzkinetics.substrspecarray.length ==0)
                obj.substratekinetics.Vm = obj.enzkinetics.sumconsts.Vm;
                obj.substratekinetics.Km = obj.enzkinetics.sumconsts.Km;
            else
                obj.substratekinetics.Vm = obj.enzkinetics.sumconsts.Vm;
                obj.substratekinetics.Km = obj.enzkinetics.sumconsts.Km;
                if(obj.enzkinetics.specificityrule==1)
                    glycanstring = glycanspecies.glycanStruct.toString;
                    if(isKey(obj.enzkinetics.substrspecdb,glycanstring))
                        ratio = obj.enzkinetics.substrspecdb(glycanstring);
                        obj.substratekinetics.Vm_ratio =  ratio.Vm;
                        obj.substratekinetics.Km_ratio =  ratio.Km;
                    end
                elseif(obj.enzkinetics.specificityrule==2)
                    % check if the substrate is the substructure of another in the database
                    for i = 1: length(obj.enz.substrspecarray)
                        substruct = obj.enzkinetics.substrspecarray.get(i).struct;
                        ratio     = obj.enzkinetics.substrspecarray.get(i).ratio;
                        if(glycanspecies.glycanStruct.contains(substruct))
                            obj.substratekinetics.Vm_ratio = ratio.Vm;
                            obj.substratekinetics.Km_ratio = ratio.Km;
                            break;
                        end
                    end
                end
            end
        end
        
        function setRxnKineticsInhibitor(obj,theCompetingRxns) % does not check the rxns are qualified reactions
            if(isa(theCompetingRxns,'Rxn'))
                rxnkinetics    = theCompetingRxns.rxnkinetics;
                inhibkinetics  = struct('Ki',rxnkinetics.substratekinetics.Km,...
                    'Ki_ratio',  rxnkinetics.substratekineticsratio.Km_ratio,...
                    'speciesid',theCompetingRxns.reac.id);
                obj.inhibitorkinetics = inhibkinetics;
            elseif(iscell(theCompetingRxns))
                numCompRxns = length(theCompetingRxns);
                for i = 1 : numCompRxns
                    rxnkinetics    = theCompetingRxns{i}.rxnkinetics;
                    inhibkinetics  = struct('Ki',rxnkinetics.substratekinetics.Km,...
                        'Ki_ratio',  rxnkinetics.substratekineticsratio.Km_ratio,...
                        'speciesid',theCompetingRxns{i}.reac.id);
                    listofsubstinhibkinetics(i,1)=inhibkinetics;
                end
                obj.inhibitorkinetics = listofsubstinhibkinetics;
            end
        end
    end
    
    % method for mathematical formula
    methods
        function setMathFormula(obj,enzid,rxnid)
            % numinhibitors = obj.inhibitorkinetics.length;
            if(isempty(obj.inhibitorkinetics))
                obj.kineticlaws.speciesid        = {obj.substratekinetics.speciesid};
                obj.kineticlaws.parameternames   = {['Vm',enzid],['Km',enzid],...
                    'Vm_ratio','Km_ratio'};
                obj.kineticlaws.parametervalues  = [obj.substratekinetics.Vm, obj.substratekinetics.Km,...
                    obj.substratekineticsratio.Vm_ratio, obj.substratekineticsratio.Km_ratio];
                obj.kineticlaws.mathformula      = 'Vm*Vm_ratio*S/(Km*Km_ratio+S)';
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'S',obj.kineticlaws.speciesid );
                
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Vm_ratio',num2str(obj.kineticlaws.parametervalues(3)));
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Km_ratio',num2str(obj.kineticlaws.parametervalues(4)));
                
                
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Vm',obj.kineticlaws.parameternames{1});
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Km',obj.kineticlaws.parameternames{2});
                
                
                
            elseif(length(obj.inhibitorkinetics)==1)
                obj.enzkinetics.mechanism
                if(strcmpi(obj.enzkinetics.mechanism,TypeStringConst.Noncomp))
                    obj.kineticlaws.mathformula      = ...
                        'Vm*Vm_ratio/((1+Ci/Ki/Ki_ratio)*(1+Km*Km_ratio/S))';
                elseif(strcmpi(obj.enzkinetics.mechanism,TypeStringConst.Uncomp))
                    obj.kineticlaws.mathformula      = ...
                        '(Vm*Vm_ratio/(1+Ci/Ki/Ki_ratio))*S/(Km*Km_ratio/(1+Ci/Ki/Ki_ratio)+S)';
                elseif(strcmpi(obj.enzkinetics.mechanism,TypeStringConst.Comp))
                    obj.kineticlaws.mathformula      = ...
                        'Vm*Vm_ratio*S/(Km*Km_ratio*(1+Ci/Ki/Ki_ratio)+S)';
                elseif(strcmpi(obj.enzkinetics.mechanism,TypeStringConst.Subst))
                    obj.kineticlaws.mathformula      = ...
                        'Vm*Vm_ratio*S/(Km*Km_ratio/(1+Ci/Ki/Ki_ratio)+S)';
                else
                    error('MATLAB:GNAT:WRONGMechanism','Incorrect Enzyme Mechanism');
                end
                
                obj.kineticlaws.speciesid        = {obj.substratekinetics.speciesid;...
                    obj.inhibitorkinetics.speciesid};
                obj.kineticlaws.parameternames   = {['Vm',enzid],...
                    ['Km',enzid],['Km',enzid],...
                    'Vm_ratio','Km_ratio','Ki_ratio'};
                obj.kineticlaws.parametervalues  = [obj.substratekinetics.Vm,...
                    obj.substratekinetics.Km,obj.inhibitorkinetics.Ki,...
                    obj.substratekineticsratio.Vm_ratio,obj.substratekineticsratio.Km_ratio,...
                    obj.inhibitorkinetics.Ki_ratio];
                
                % replace ratio value in kinetic laws
                if(obj.kineticlaws.parametervalues(4)~=0)
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Vm_ratio',num2str(obj.kineticlaws.parametervalues(4)));
                else
                    obj.kineticlaws.mathformula = '0';
                    return;
                end                
                
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Km_ratio',num2str(obj.kineticlaws.parametervalues(5)));
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Ki_ratio',num2str(obj.kineticlaws.parametervalues(6)));
                
                % replace names
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'S', obj.kineticlaws.speciesid{1});
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Ci',obj.kineticlaws.speciesid{2});
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Vm',obj.kineticlaws.parameternames{1});
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Km',obj.kineticlaws.parameternames{2});
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'Ki',obj.kineticlaws.parameternames{3});
            elseif(length(obj.inhibitorkinetics)>1)
                if(strcmpi(obj.enzkinetics.mechanism,TypeStringConst.Comp))
                    obj.kineticlaws.speciesid        = {obj.substratekinetics.speciesid};
                    obj.kineticlaws.parameternames   = {['Vm',enzid],...
                        ['Km',enzid]};
                    obj.kineticlaws.parametervalues  = [obj.substratekinetics.Vm,...
                        obj.substratekinetics.Km];
                    
                    % construct mathematical formula
                    obj.kineticlaws.mathformula = 'Vm*Vm_ratio*S/(S+Km*Km_ratio*(1';
                    
                    % store ratio to parameter sets
                    obj.kineticlaws.parameternames{end+1}  = 'Vm_ratio';
                    obj.kineticlaws.parameternames{end+1}  = 'Km_ratio';
                    obj.kineticlaws.parametervalues(end+1) =  obj.substratekineticsratio.Vm_ratio;
                    obj.kineticlaws.parametervalues(end+1) =  obj.substratekineticsratio.Km_ratio;
                    
                    % replace ratio with values
                    if(obj.kineticlaws.parametervalues(4)~=0)
                        obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                          'Vm_ratio',num2str(obj.kineticlaws.parametervalues(3)));
                    else
                        obj.kineticlaws.mathformula ='0';
                        return;
                    end
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Km_ratio',num2str(obj.kineticlaws.parametervalues(4)));
                    
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'S',  obj.kineticlaws.speciesid{1});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Vm', obj.kineticlaws.parameternames{1});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Km', obj.kineticlaws.parameternames{2});
                    
                    % add inhibitor kinetics to formula
                    for i = 1 : length(obj.inhibitorkinetics)
                        obj.kineticlaws.speciesid{end+1}       = obj.inhibitorkinetics(i,1).speciesid;
                        obj.kineticlaws.mathformula            = [obj.kineticlaws.mathformula, ...
                            '+', obj.kineticlaws.speciesid{end}];
                        
                        obj.kineticlaws.parameternames{end+1}  = ['Km',enzid];
                        obj.kineticlaws.parametervalues(end+1) = obj.inhibitorkinetics(i,1).Ki;
                        
                        obj.kineticlaws.mathformula            = [obj.kineticlaws.mathformula ...
                            '/' obj.kineticlaws.parameternames{end}];
                        
                        obj.kineticlaws.parameternames{end+1}  = ['Km_ratio',enzid];
                        obj.kineticlaws.parametervalues(end+1) = obj.inhibitorkinetics(i,1).Ki_ratio;
                        
                        obj.kineticlaws.mathformula            = [obj.kineticlaws.mathformula ...
                            '/' num2str(obj.kineticlaws.parametervalues(end))];
                    end
                    obj.kineticlaws.mathformula=[obj.kineticlaws.mathformula,'))'];
                else
                    error('MATLAB:GNAT:MECHANISMNOTSUPPORTED','NOT SUPPORTED MECHANISM');
                end
            end
            
           obj.kineticlaws.mathformula  = regexprep(obj.kineticlaws.mathformula,'*1(?=[^0-9])','');
           obj.kineticlaws.mathformula  = regexprep(obj.kineticlaws.mathformula ,'\1(?=[^0-9])','');
           obj.kineticlaws.mathformula  = regexprep(obj.kineticlaws.mathformula ,'/1(?=[^0-9])','');
        end
    end
    
    methods
        function new = clone(obj)
            
            % Instantiate new object of the same class.
            new = feval(class(obj));
            
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                if(isa(obj.(p{i}),'handle'))
                    if(isa(obj.(p{i}),'CellArrayList'))
                        new.(p{i})=CellArrayList;
                        for j = 1 : obj.(p{i}).length
                            new.(p{i}).add(obj.(p{i}).get(j))
                        end
                    elseif(isa(obj.(p{i}),'containers.Map'))
                        allkeys    = obj.(p{i}).keys;
                        allvalues  = obj.(p{i}).values;
                        new.(p{i}) = containers.Map(allkeys,allvalues);
                    else
                        new.(p{i}) = obj.(p{i}).clone;
                    end
                else
                    new.(p{i}) = obj.(p{i});
                end
            end
        end
    end
end

