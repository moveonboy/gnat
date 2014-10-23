function obj = bracketResidue(obj,residuetobracket)
%removeResidue remove residue from the structure
%
% See also ADDRESIDUES
for j = 1 : length(residuetobracket)
    jthresiduetobracket = residuetobracket(j,1);
%     if(~isa(jthresiduetobracket,'GlycanResidue'))
%         errorReport(mfilename,'IncorrectInputType');
%     end
    % change parent linkage to the residue to be removed
    linkageParent =  jthresiduetobracket.getLinkageParent;
    parentResidue =  linkageParent.getParent;
    
    linkageChildren = parentResidue.getLinkageChildren;
    nChildren       = length(linkageChildren);
    
    i=1;
    isFindChild=false;
    terminalResiduePos =-1;
    while((i<=nChildren) &&(~isFindChild))
        ithchild = linkageChildren(i,1).getChild;
        if(ithchild==jthresiduetobracket)
            isFindChild=true;
            terminalResiduePos = i;
        end
        i=i+1;
    end
    if(terminalResiduePos==-1)
        error('MATLAB:GNAT:ERRORRESIDUE','RESIDUE INFO INCORRECT');
    end
    parentResidue.unsetLinkageChildren(terminalResiduePos);
end

obj.addBracket;
if(length(residuetobracket)>1)
    obj.addResiduesToBracket(residuetobracket);
else
    obj.addResidueToBracket(residuetobracket);
end

% update name
obj.glycanjava = obj.structMat2Java;
obj.name = char(obj.glycanjava.toStringOrdered(0));
end
        
