function residueToRootChar = getresiduetoroot(obj,residue,varargin)
%getresiduetoroot return the string starting from specificied residue to the root
%    getresiduetoroot(obj,residue)
%
% See also removeNonRedEndResidue

% Author: Gang Liu
% Date Lastly Updated: 7/15/2014
if(length(varargin)==1)&&(strcmpi(varargin{1},'nolinkage'))
    nolinkage=1;
else
    nolinkage=0;
end

iterresidue = residue;
residueToRootChar='';
count = 0;
isParent = ~isempty(iterresidue.linkageParent);
while(isParent)
    residueTypeChar =  [iterresidue.stereoConfig.symbol,'-'];
    resType = iterresidue.residueType;
    if(~strcmpi(resType.name,'freeEnd'))
        residueTypeChar = [residueTypeChar,resType.name];
        residueTypeChar = [residueTypeChar,','];
        residueTypeChar = [residueTypeChar,iterresidue.ringType.rSize];
    else
        residueTypeChar=resType.name;
    end
    
    if(~isempty( iterresidue.linkageParent))
        if(~nolinkage)
            LinkageChar = iterresidue.linkageParent.bonds.posParent;
            LinkageChar = [LinkageChar,...
                iterresidue.anomer.symbol];
            LinkageChar = [LinkageChar,...
                iterresidue.anomer.carbonPos];
        else
            LinkageChar =iterresidue.anomer.symbol;
        end
    else
        LinkageChar='';
    end
    
    residueLinkageChar= [LinkageChar,residueTypeChar];
    residueToRootChar = [residueLinkageChar,residueToRootChar];
    isParent = ~isempty(iterresidue.linkageParent);
    if(isParent)
        iterresidue = iterresidue.linkageParent.parent;
        residueToRootChar = ['--',residueToRootChar];
    end;
    count =count +1;
end
end





