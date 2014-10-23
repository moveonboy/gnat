classdef RxnKinetics < handle
% RXNKINETICS describe reaction rate velocity
%    
% See also Rxn, EnzKinetics.

% Author: Gang Liu
% Date Lastly Updated: 4/11/2014
    methods(Abstract)
       setMathFormula(obj)
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

