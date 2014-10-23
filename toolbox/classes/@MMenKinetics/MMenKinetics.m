classdef MMenKinetics < EnzKinetics

    properties
        % Michaelis Menten Constants
        eleconsts   = struct('kf',[],'kr',[],'kcat',[],'enzconc',[]);
        sumconsts   = struct('Vm',[],'Km',[]);
        
        %Substrate Inhibition Constants;
        substinhiconst;
        substsumconst;
        
        %Inhibitor Rate Constants
        kieleconsts     = CellArrayList; % each element has three fields: 1)speciesname;2)kfi;3)kri;
        kisumconsts     = CellArrayList; % each element has two fields: 1)speciesname;2)Ki;
        
        form='';
        mechanism='';
        substrspecdb;%    = containers.Map;
        substrspecarray = CellArrayList;
        specificityrule;
    end
    
    methods
        function setSubstrateSpecificityK(obj,specistruct,specifiratio,specifictyrule)
            % Input: specistruct is a CellArrayList object containing a list
            % of glycan structure. specifiratio is a CellArrayList object
            % containing a list of ratio over basal value among all rate
            % constants
            %  Example: specistruct= CellArrayList;
            %           specistruct.add(m5gn);
            %           specistruct.add(m5gngn);
            %           specifiratio = CellArrayList;
            %           ratio1 = struct('vm',1,'km',1);
            %           specifiratio.add(ratio1);
            %           setSubstrateSpecificityK(enzObj,specistruct,specifiratio,1);
            %
            % See also MMenKinetics.
            if(length(specistruct)~=length(specifiratio))
                error('MATLAB:GNAT:ERRORINPUT',...
                'LENGTH OF ARGUMENT 1SHOULD BE EQUAL TO ARGUMENT 2');
            end
            
            if(~isnumeric(specifictyrule) && ~ischar(specifictyrule))
                error('MATLAB:GNAT:ERRORINPUT',...
                    'Wrong input for enzyme kinetic specificity rule');
            end
            
            for i = 1 : length(specistruct)
                glycanstruct = specistruct.get(i).toString;
                obj.substrspecdb(glycanstruct) = specifiratio.get(i);
            end
            
            obj.substrspecarray = CellArrayList;
            for i = 1 : length(specistruct)
                specificity.struct     = specistruct.get(i);
                specificity.ratio      = specifiratio.get(i);
                obj.substrspecarray.add(specificity);
            end
            
            obj.specificityrule = specifictyrule; % either 1,2 or string
        end
    end
    
    methods
        function obj=MMenKinetics(varargin)
            if(isempty(varargin))
              return;
            end
            
            rateconsts     = varargin{3};            
            obj.mechanism  = varargin{1};
            obj.form       = varargin{2};
            if(strcmpi(obj.form,'elemental'))
                try
                    obj.eleconsts.kf=rateconsts.kf;
                    obj.eleconsts.kr=rateconsts.kr;
                    obj.eleconsts.kcat=rateconsts.kcat;
                    obj.eleconsts.enzconc=rateconsts.enzconc;
                catch
                    error('MATLAB:GNAT:INCORRECTINPUT','INCORRECT INPUT')
                end
            elseif(strcmpi(obj.form,'ensemble'))
                try
                    obj.sumconsts.Vm=rateconsts.Vm;
                    obj.sumconsts.Km=rateconsts.Km;
                catch
                    error('MATLAB:GNAT:INCORRECTINPUT','INCORRECT INPUT')
                end
            else
                error('MATLAB:GNAT:NOTSUPPORTEDFORM','This form is not supported')
            end
            
            obj.substrspecdb=containers.Map;            
            if(length(varargin)==5)
                setSusbstrInhibition(obj,varargin{4},varargin{5});
            end
        end
        
        function setSusbstrInhibition(obj,form,inhibitconsts)
            if(strcmpi(obj.mechanism,TypeStringConst.Subst))
                if(strcmpi(form,'elemental'))
                    obj.substinhiconst.kfi     = inhibitconsts.kfi;
                    obj.substinhiconst.kri     = inhibitconsts.kri;
                    obj.numinhibitors          = obj.numinhibitors+1;
                elseif(strcmpi(form,'ensemble'))
                    obj.substsumconst.Ki       = inhibitconsts.Ki;
                    obj.numinhibitors         = obj.numinhibitors+1;
                else
                    error('MATLAB:GNAT:NOTSUPPORTEDFORM','This form is not supported')
                end
            end
        end
        
        function setInhibitor(obj,form,inhibitconsts)  %
            if(strcmpi(obj.mechanism,TypeStringConst.Subst))
                setSusbstrInhibition(obj,form,inhibitconsts)
                return;
            elseif(strcmpi(obj.mechanism,TypeStringConst.Noncomp)...
                   ||strcmpi(obj.mechanism,TypeStringConst.Uncomp)||...
                   strcmpi(obj.mechanism,TypeStringConst.Comp));
                if(strcmpi(form,'elemental'))
                    for i = 1 : length(inhibitconsts)
                        ithinhibitconst= inhibitconsts.get(i);
                        kieleconst.speciesname = ithinhibitconst.speciesname;
                        kieleconst.kfi = ithinhibitconst.kfi;
                        kieleconst.kri = ithinhibitconst.kri;
                        obj.kieleconsts.add(kieleconst);
                        obj.numinhibitors   =obj.numinhibitors   +1;
                    end
                elseif(strcmpi(form,'ensemble'))
                    for i = 1 : length(inhibitconsts)
                        ithinhibitconst= inhibitconsts.get(i);
                        kieleconst.speciesname = ithinhibitconst.speciesname;
                        kieleconst.Ki = ithinhibitconst.Ki;
                        obj.kieleconsts.add(kieleconst);
                        obj.numinhibitors   =obj.numinhibitors   +1;
                    end
                else
                    error('MATLAB:GNAT:NOTSUPPORTEDFORM','This form is not supported')
                end
            end
            
            if(strcmpi(form,'elemental'))
                obj.ele2ems(obj.mechanism)
            end
        end
        
        function kineticLaws = toMathExprWithSub(obj,substrateid,enzrxnpostfix)
            kineticLaws                       = obj.toMathExpr;
            kineticLaws.mathexpr              = regexprep(kineticLaws.mathexpr,'*S/',['*',substrateid,'/']);
            kineticLaws.mathexpr              = regexprep(kineticLaws.mathexpr,'+S)',['+',substrateid,')']);
            kineticLaws.speciesnames          = regexprep(kineticLaws.speciesnames,'S',substrateid);
            newvmname                         = ['Vm_' enzrxnpostfix];
            newkmname                         = ['Km_' enzrxnpostfix];
            kineticLaws.parameter.(newvmname) = kineticLaws.parameter.Vm;
            kineticLaws.parameter.(newkmname) = kineticLaws.parameter.Km;            
            kineticLaws.parameter             = rmfield(kineticLaws.parameter,'Vm');
            kineticLaws.parameter             = rmfield(kineticLaws.parameter,'Km');
        end
        
        function kineticlaws = toMathExpr(obj)
            if(obj.numinhibitors<=1)
                if(strcmpi(obj.mechanism,TypeStringConst.SMM))
                    kineticlaws.mathexpr = 'Vm*S/(Km+S)';
                    kineticlaws.speciesnames = {'S'};
                    kineticlaws.parameternames   = {'Vm','Km'};
                    kineticlaws.parametervalues  = [obj.sumconsts.Vm, obj.sumconsts.Km];
                elseif(strcmpi(obj.mechanism,TypeStringConst.Noncomp))
                    kineticlaws.sbmlmath         = 'Vm/((1+Ci/Ki)*(1+Km/S))';
                    kineticlaws.speciesnames     = {'S','Ci'};
                    kineticlaws.parameternames   = {'Vm','Km','Ki'};
                    kineticlaws.parametervalues  = [obj.sumconsts.Vm, obj.sumconsts.Km,obj.sumconsts.Ki];                  
                elseif(strcmpi(obj.mechanism,TypeStringConst.Uncomp))
                    kineticlaws.sbmlmath = '(Vm/(1+Ci/Ki))*S/(Km/(1+Ci/Ki)+S)';
                    kineticlaws.speciesnames = {'S','Ci'};
                     kineticlaws.parameternames   = {'Vm','Km','Ki'};
                    kineticlaws.parametervalues  = [obj.sumconsts.Vm, obj.sumconsts.Km,obj.sumconsts.Ki]; 
                elseif(strcmpi(obj.mechanism,TypeStringConst.Comp))
                    kineticlaws.sbmlmath = 'Vm*S/(Km*(1+Ci/Ki)+S)';
                    kineticlaws.speciesnames = {'S','Ci'};
                     kineticlaws.parameternames   = {'Vm','Km','Ki'};
                    kineticlaws.parametervalues  = [obj.sumconsts.Vm, obj.sumconsts.Km,obj.sumconsts.Ki]; 
                elseif(strcmpi(obj.mechanism,TypeStringConst.Subst))
                    kineticlaws.sbmlmath = 'Vm*S/(Km/(1+Ci/Ki)+S)';
                    kineticlaws.speciesnames = {'S','Ci'};
                    kineticlaws.parameternames   = {'Vm','Km','Ki'};
                    kineticlaws.parametervalues  = [obj.sumconsts.Vm, obj.sumconsts.Km,obj.sumconsts.Ki]; 
                else
                    error('MATLAB:GNAT:WRONGMechanism','Incorrect Enzyme Mechanism');
                end
            else
                if(strcmpi(obj.mechanism,TypeStringConst.Comp))
                    kineticlaws.speciesnames = {'S','Ci1'};
                    % kineticlaws.rtconsts = {'Vm','Km','Ki1'};
                    kineticlaws.parameter    = struct('Vm',obj.sumconsts.Vm,'Km',obj.sumconsts.Km);
                    kineticlaws.mathexpr='Vm*S/(S+Km(1+';
                    for i = 1 : length(obj.kisumconsts)
                        ele1=['+Ci_',obj.kisumconsts.get(i).speciesname];
                        ele2=['/Ki_',obj.kisumconsts.get(i).speciesname];
                        ele=[ele1, ele2];
                        kineticlaws.mathexpr = [kineticlaws.mathexpr,ele];
                        constant=['Ki',obj.kisumconsts.get(i).speciesname];
                        inhibitor=['Ci',obj.kisumconsts.get(i).speciesname];
                        kineticlaws.rtconsts{end+1}=constant;
                        kineticlaws.speciesnames{end+1}=inhibitor;
                    end
                    kineticlaws.mathexpr=[kineticlaws.mathexpr,'))'];
                end
            end
        end
        
        function ele2ems(obj)
            if(strcmpi(obj.mechanism,TypeStringConst.SMM))
                obj.sumconsts.Vm=obj.eleconsts.kcat*obj.eleconsts.enzconc;
                obj.sumconsts.Km=obj.eleconsts.kr/obj.eleconsts.kf;
            elseif(strcmpi(obj.mechanism,TypeStringConst.Subst))
                obj.sumconsts.Vm=obj.eleconsts.kcat*obj.eleconsts.enzconc;
                obj.sumconsts.Km=obj.eleconsts.kr/obj.eleconsts.kf;
                obj.substsumconst.Ki=obj.substinhiconst.kri/obj.substinhiconst.kfi;
            elseif (strcmpi(obj.mechanism,TypeStringConst.Uncomp)...
                    || strcmpi(obj.mechanism,TypeStringConst.Noncomp)...
                    || strcmpi(obj.mechanism,TypeStringConst.Subst));
                obj.sumconsts.vm=obj.eleconsts.kcat*obj.eleconsts.enzconc;
                obj.sumconsts.km=obj.eleconsts.kr/obj.eleconsts.kf;
                for i = 1 : length(obj.kieleconsts)
                    obj.kisumconsts.speciesname = obj.kieleconsts.get(i).speciesname;
                    obj.kisumconsts.Ki = obj.kieleconsts.get(i).kri/obj.kieleconsts.get(i).kfi;
                end
            end
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
                        if(~isempty(allkeys))
                          new.(p{i}) = containers.Map(allkeys,allvalues);
                        end
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

