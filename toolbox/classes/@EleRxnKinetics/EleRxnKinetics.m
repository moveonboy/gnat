classdef EleRxnKinetics < RxnKinetics
    
    properties
        kineticsconst = struct('kf',[],'kr',[]);
        reac = CellArrayList; % each element has two fields: 1)speciesid;2)ordercoeffcient;
        prod = CellArrayList; % each element has two fields: 1)speciesid;2)ordercoeffcient;
        kineticlaws;
        
    end
    
    methods % Constructor for EleRxnKinetics define rate constants and oder
        function obj=EleRxnKinetics(speciesR,speciesP,rateconsts)
            obj.reac = speciesR;
            obj.prod = speciesP;
            obj.kineticsconst.kf = rateconsts.kf;
            obj.kineticsconst.kr = rateconsts.kr;
        end
    end    
    
    % method for mathematical formula
    methods
        function setMathFormula(obj,rxnid)
            obj.kineticlaws.speciesid       = {};
            obj.kineticlaws.parameternames  = {['kf','_',rxnid],['kr','_',rxnid]};
            obj.kineticlaws.parametervalues = [obj.kineticsconst.kf,obj.kineticsconst.kr];
            for i = 1 : length(obj.reac)
                coeffcientid                           = ['a' num2str(i) '_' rxnid];
                obj.kineticlaws.parameternames{end+1}  = coeffcientid;
                coeffcientvalue                        = obj.reac.get(i).coef;
                obj.kineticlaws.parametervalues(end+1) = coeffcientvalue;
                obj.kineticlaws.speciesid{end+1}       = obj.reac.get(i).speciesid;
            end
            for j = 1 : length(obj.prod)
                coefid                           = ['b' num2str(j) '_' rxnid];
                obj.kineticlaws.parameternames{end+1}  = coefid;
                coefvalue                        = obj.reac.get(j).coef;
                obj.kineticlaws.parametervalues(end+1) = coefvalue;
                obj.kineticlaws.speciesid{end+1}       = obj.prod.get(j).speciesid;
            end
            obj.kineticlaws.mathformula = 'kf';
            obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                'kf', obj.kineticlaws.parameternames{1});
            for k = 1 : length(obj.reac)
                reactant                    = obj.reac.get(k).speciesid;
                coeffcient                  = ['a' num2str(k) '_' rxnid];
                ele1                        = ['(' reactant '^' coeffcient ')'];
                obj.kineticlaws.mathformula = [obj.kineticlaws.mathformula '*' ele1];
            end
            obj.kineticlaws.mathformula = [obj.kineticlaws.mathformula '-' 'kr'];
            for l = 1 : length(obj.prod)
                product                     =  obj.prod.get(l).speciesid;
                coef                        = ['b' num2str(l) '_' rxnid];
                ele2                         = ['(' product '^' coef ')'];
                obj.kineticlaws.mathformula = [obj.kineticlaws.mathformula '*' ele2];
            end
            obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                'kr', obj.kineticlaws.parameternames{2});
        end
    end
    
end