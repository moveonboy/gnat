classdef EnzKinetics < handle
    properties
        enz;
    end
    
    methods (Abstract)
        kineticlaws = toMathExpr(obj,varargin);
      %  vel         = simVel(obj,varargin);
    end    
end

