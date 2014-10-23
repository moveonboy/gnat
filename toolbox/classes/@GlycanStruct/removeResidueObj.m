function isResidueDeleted = removeResidueObj(obj,residuetoremove)
%removeResidueObj remove residue from the structure
%
% See also ADDRESIDUES
if(~isa(residuetoremove,'GlycanResidue'))
    errorReport(mfilename,'IncorrectInputType');
end
% change parent linkage to the residue to be removed
linkageParent =  residuetoremove.getLinkageParent;
parentResidue =  linkageParent.getParent;

linkageChildren = parentResidue.getLinkageChildren;
nChildren       = length(linkageChildren);

i=1;
isFindChild=false;
terminalResiduePos =-1;
while((i<=nChildren) &&(~isFindChild))
    ithchild = linkageChildren(i,1).getChild;
    if(ithchild==residuetoremove)
        isFindChild=true;
        terminalResiduePos = i;
    end
    i=i+1;
end
if(terminalResiduePos==-1)
    isResidueDeleted=0;
end
parentResidue.unsetLinkageChildren(terminalResiduePos);
isResidueDeleted=1;

% update name
obj.glycanjava = obj.structMat2Java;
obj.name = char(obj.glycanjava.toStringOrdered(0));
end
        
