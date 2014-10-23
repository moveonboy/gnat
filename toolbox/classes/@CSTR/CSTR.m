classdef CSTR < Reactor
    % CSTR continuous flow stirred-tank reactor
    %
    %See also PFR,BatchReactor.
    
    % Author: Gang Liu
    % Date Lastly Updated: 03/26/14
    methods  % constructor
        function obj = CSTR(varargin)
          if(length(varargin)==3)
             obj.qflowin    = varargin{1};
             obj.qflowin    = varargin{2};
             obj.ktransport = varargin{3};
          elseif(length(varargin)==2)
             obj.qflowin    = varargin{1};
             obj.qflowout   = varargin{1};
             obj.ktransport = varargin{2};
          elseif(length(varargin)==4)
             obj.qflowin     = varargin{1};
             obj.qflowout    = varargin{2};
             obj.ktransport = varargin{3};
             obj.volume     = varargin{4};
          end
        end
    end
    
    methods 
        function setSpeciesInletDB(obj,speciesInletConc)
            obj.speciesidinletdb = containers.Map;
            for i = 1 : length(speciesInletConc)
                if(isa(speciesInletConc.get(i).species,'Species'))
                   obj.speciesidinletdb(speciesInletConc.get(i).species.id)=...
                        speciesInletConc.get(i).inletconc;
                elseif(ischar(speciesInletConc.get(i).species))
                    obj.speciesidinletdb(speciesInletConc.get(i).species)=...
                        speciesInletConc.get(i).inletconc;
                else
                    error('MATLAB:GNAT:ERRORINPUT','ERROR FIELD TYPE');
                end                
            end
        end
    end    
end