function glycanJava = structMat2Java(obj)
%structMat2Java convert a MATLAB GlycanStruct object to a Java Glycan object
%   
%  GLYCANJAVA = structMat2Java(GlycanStructobj) reads the GlycanStruct
%  object and returns a Java Glycan object. 
%  
%  GLYCANJAVA = GlycanStructobj.structMat2Java is equivilent to GLYCANJAVA
%  = structMat2Java(GlycanStructobj).
%
% See also GLYCANSTRUCT/GLYCANSTRUCT

% Author: Gang Liu
% Date Lastly Updated: 7/3/14

if(isempty(obj))
    glycanJava = [];
    return 
end

import org.eurocarbdb.application.glycanbuilder.*;

glycanJava = Glycan();
if(~isempty(obj.getRoot))
    root = residueMat2Java(obj.getRoot);
    glycanJava.setRoot(root);
end

if(~isempty(obj.bracket))
    bracket = obj.bracket;
    for i = 1 : length(bracket.linkageChildren)
        ithchild = bracket.linkageChildren(i,1).child;
        ithchildjava = residueMat2Java(ithchild);
        glycanJava.addAntenna(ithchildjava);
    end
end

if(~isempty(obj.getMassOptions))
    mass_opt = massOptMat2Java(obj.getMassOptions);
    glycanJava.setMassOptions(mass_opt);
end

end

