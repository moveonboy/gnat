classdef BiBiKinetics < EnzKinetics
   properties
        %Krambeck Model
        sumconsts    = struct('Vm',[],'Kmd_G1E',[],'Kmd_SNE',[],'Keq',[]);
        %Substrate Inhibition Constants
        substsumconst
        
        %Inhibitor Rate Constants
        kisumconsts     = CellArrayList;  %each element has two fields: 1)speciesname;2)Ki;
        
        mechanism       = '';
        substrspecdb;
        substrspecarray = CellArrayList;
        specificityrule;
    end
    
    methods
        function setSubstSpecificityK(obj,specistruct,specifiratio,specifictyrule)
            % Input: specistruct is a CellArrayList object containing a list
            % of glycan structure. specifiratio is a CellArrayList object
            % containing a list of ratio over basal value among all rate
            % constants
            %  Example: specistruct= CellArrayList;
            %           specistruct.add(m5gn);
            %           specifiratio = CellArrayList;
            %           ratio1 = struct('vm',1,'Kmd_G1E',1,'Kmd_SNE',1,'Keq',1);
            %           specifiratio.add(ratio1);
            %           setSubstSpecificityK(enzObj,specistruct,specifiratio,1);
            %
            % See also BiBiKinetics.
            if(length(specistruct)~=length(specifiratio))
                error('MATLAB:GNAT:ERRORINPUT',...
                    'LENGTH OF ARGUMENT 1SHOULD BE EQUAL TO ARGUMENT 2');
            end
            
            if(~isnumeric(specifictyrule) && ~ischar(specifictyrule))
                error('MATLAB:GNAT:ERRORINPUT',...
                    'Wrong input for enzyme kinetic specificity rule');
            end
            
            myKeys  ={};
            myValue ={};
            for i = 1 : length(specistruct)
                glycanstruct = specistruct.get(i).toString;
                myKeys{end+1}=glycanstruct;
                myValue{end+1}=specifiratio.get(i);
                obj.substrspecdb=containers.Map(myKeys,myValue);
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
        function obj=BiBiKinetics(varargin)
            % enzobj,mechanism,rateconsts,varargin) % constructor
            if(isempty(varargin))
                return;
            end
            obj.mechanism  = varargin{1};
            rateconsts     = varargin{2};
            if strcmpi(obj.mechanism,TypeStringConst.bibiSeqOrderwSubstInhib)
                obj.sumconsts.Vm          = rateconsts.Vm;
                obj.sumconsts.Kmd_G1E     = rateconsts.Kmd_G1E;
                obj.sumconsts.Kmd_SNE     = rateconsts.Kmd_SNE;
                obj.sumconsts.Keq         = rateconsts.Keq;
                
                if(length(varargin)==3)
                  inhibitconsts = varargin{3};
                  setSusbstInhibition(obj,inhibitconsts)
                end
            else
                
                
            end
        end        
        
        function setSusbstInhibition(obj,inhibitconsts)
            if (strcmpi(obj.mechanism,TypeStringConst.bibiSeqOrderwSubstInhib))
                for i = 1 : length(inhibitconsts)
                    ithSubstinhibitconst = inhibitconsts.get(i);
                    obj.kisumconsts.speciesname = ithSubstinhibitconst.speciesname;
                    obj.kisumconsts.Ki = ithSubstinhibitconst.Ki;
                    obj.kisumconsts.add(obj.kisumconsts);
                 end
            end
            
        end
        
         function kineticlaws = toMathExpr(obj)
                 kineticlaws.speciesnames = {'G1','SN','G2','NT','E'};
                 kineticlaws.rtconsts     = {'Vm','kmd_G1E','kmd_SNE','Keq'};
                 kineticlaws.mathexpr = 'Vm*(SN*G1-(G2*NT)/Keq)/kmd_G1E/(kmd_SNE+SN*';
                 
                 if(obj.numcomptsubst>0)
                     kineticlaws.mathexpr = [kineticlaws.mathexpr,'(1'];
                     for i=1:length(obj.addkiconsnts)
                         inhconsnts=['Kmi',num2str(i)];
                         kineticlaws.rtconsts{end+1}=inhconsnts;
                         inhglycans=['Gi',num2str(i+2)];
                         kineticlaws.speciesnames{end+1} = inhglycans;
                         kineticlaws.mathexpr            = [kineticlaws.mathexpr,'+',B,'/',A];
                     end
                         kineticlaws.mathexpr                = [kineticlaws.mathexpr,')'];
                 end
         end
    end
end

