classdef TransportKinetics < RxnKinetics
    %TRANSPORTKINETICS describe transport Kinetics
    % 
    %See also RxnKinetics
    
    % Author: Gang Liu
    % Date Lastly Updated: 3/29/14
    
    properties(Constant)
        firstorder='First Order';
        zeroorder='Zero Order';        
    end
    
    properties
        transportorder;
    end
    
    properties
        firstorderkinetics=struct('kt',[],'speciesid',[]);
        zeroorderkinetics=struct('qin',[],'speciesid',[],'initconc',[]);
        kineticlaws
    end
    
    methods
        function obj=TransportKinetics(transorder,varargin)
            obj.transportorder = transorder;
            if(strcmpi(transorder,TransportKinetics.firstorder))
                if( length(varargin)~=2 ||...
                    ~isnumeric(varargin{1}) ||...
                    ~ischar(varargin{2}))
                    error('MATLAB:GNAT:INCORRECTINPUT',...
                    'INCORRECT TRANSPORT INPUT') 
                end
                obj.firstorderkinetics.kt        = varargin{1};
                obj.firstorderkinetics.speciesid = varargin{2};
            elseif(strcmpi(transorder,TransportKinetics.zeroorder))
                if( length(varargin)~=3 ||...
                    ~isnumeric(varargin{1})||...
                    ~isnumeric(varargin{3})||...
                    ~ischar(varargin{2}))
                    error('MATLAB:GNAT:INCORRECTINPUT',...
                    'INCORRECT TRANSPORT INPUT') 
                end
                obj.zeroorderkinetics.qin      = varargin{1};
                obj.zeroorderkinetics.speciesid = varargin{2};
                obj.zeroorderkinetics.initconc  = varargin{3};
            else
                error('MATLAB:GNAT:INCORRECTTRANSPORTKINETICS',...
                    'INCORRECT TRANSPORT KINETICS')
            end
        end    
    end
    
    methods
        function setMathFormula(obj,rxnid)
            if(strcmpi(obj.transportorder,TransportKinetics.firstorder))
                obj.kineticlaws.mathformula ='kt*S';
                obj.kineticlaws.speciesid   = obj.firstorderkinetics.speciesid;
                obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                    'S',obj.kineticlaws.speciesid);
                obj.kineticlaws.parameternames  = {'kt'};
                obj.kineticlaws.parametervalues = obj.firstorderkinetics.kt;
            elseif(strcmpi(obj.transportorder,TransportKinetics.zeroorder))
                obj.kineticlaws.mathformula ='qin*Si';
                obj.kineticlaws.parameternames  = {'qin'};
                obj.kineticlaws.parameternames{end+1}  = ['Si' obj.zeroorderkinetics.speciesid];
                obj.kineticlaws.parametervalues = [obj.zeroorderkinetics.qin;...
                       obj.zeroorderkinetics.initconc];
                if(obj.kineticlaws.parametervalues(2)~=0)
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Si',num2str(obj.kineticlaws.parametervalues(2)));
                else
                    obj.kineticlaws.mathformula = '0';
                end
            else
                error('MATLAB:GNAT:INCORRECTTRANSPORTKINETICS',...
                      'INCORRECT TRANSPORT KINETICS')            
            end
        end
    end   
end

