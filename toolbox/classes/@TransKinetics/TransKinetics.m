classdef TransKinetics < RxnKinetics
    %TRANSKINETICS derive rate law formula 
    %   
    %
    %See also RxnKinetics.
    
    properties
        reactionorder; % zero or first oder
        firstorderkinetics = struct('kt',[],'speciesid',[]);
        zeroorderkinetics  = struct('fin',[]);
        kineticlaws
    end
    
    methods  % Constructor for MMenRxnKinetics define substrate specific rate constants
        function obj=TransKinetics(reactionorder,varargin)
            obj.reactionorder = reactionorder;
            if(reactionorder==1)
               if(~ischar(varargin{1})|| ~isnumeric(varargin{2}))
                 error('MATLAB:GNAT:ERRORINPUT','INCORRECT INPUT FOR TRANSPORT KINETICS');
               end
                
               obj.firstorderkinetics.speciesid = varargin{1};
               obj.firstorderkinetics.kt        = varargin{2};
            elseif(reactionorder==0)
               if(~isnumeric(varargin{1}))
                 error('MATLAB:GNAT:ERRORINPUT','INCORRECT INPUT FOR TRANSPORT KINETICS');
               end
                
               obj.zeroorderkinetics = varargin{1};
            else
               error('MATLAB:GNAT:ERRORINPUT','INCORRECT INPUT FOR TRANSPORT KINETICS');
            end
        end
    end
    
    
    % method for mathematical formula
    methods
        function setMathFormula(obj,rxnid)
            if(obj.reactionorder==1)
                if isempty(obj.firstorderkinetics.speciesid)...
                  || isempty(obj.firstorderkinetics.kt)
                  error('MATLAB:GNAT:ERRORINPUT','EMPTY INPUT FOR TRANSPORT KINETICS');
                end
                
                obj.kineticlaws.speciesid        = {obj.firstorderkinetics.speciesid};
                obj.kineticlaws.parameternames   = {'kt_',rxnid};
                obj.kineticlaws.parametervalues  = [obj.firstorderkinetics.kt];
                obj.kineticlaws.mathformula      = 'kt*S';
                
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'kt',obj.kineticlaws.speciesid{1});
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'S',obj.kineticlaws.parameternames{1});
              
            elseif(obj.reactionorder==0)
                if isempty(obj.zeroorderkinetics.fin)
                  error('MATLAB:GNAT:ERRORINPUT','EMPTY INPUT FOR TRANSPORT KINETICS');
                end
                
                obj.kineticlaws.speciesid        = [];
                obj.kineticlaws.parameternames   = {'fin_',rxnid};
                obj.kineticlaws.parametervalues  = [obj.zeroorderkinetics.fin];
                obj.kineticlaws.mathformula      = 'fin';
                
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'fin',obj.kineticlaws.parameternames{1});
           end
        end
    end
    
end

