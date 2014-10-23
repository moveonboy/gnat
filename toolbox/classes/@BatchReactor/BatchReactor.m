classdef BatchReactor < Reactor
    methods
        function obj = BatchReactor(varargin)
            if(length(varargin)==1)
                obj.setspeciesid(obj,varargin{1})
            end
            obj.setqflowin;
            obj.setqflowout;
            obj.setktransport;
        end
        
        function setspeciesid(obj,speciesinitconc)
            for i = 1 : speciesinitconc.length
                species.id       = speciesinitconc.get(i).id;
                species.initconc = speciesinitconc.get(i).initconc;
                obj.speciesinitconc.add(species)
            end
        end
        
        function setqflowin(obj,varargin)
            obj.qflowin    = 0;
        end
        
        function setqflowout(obj,varargin)
            obj.qflowout   = 0;
        end
        
        function setktransport(obj,varargin)
            obj.ktransport = 0;
        end        
    end
end