function branchdepth = getBranchDepth(obj,terminalresidue)
%getBranchDepth get depth of the branch
%
% See also getDepth

branchdepth = 0;
hasParent   = ~isempty(terminalresidue.linkageParent);
iterresidue = terminalresidue;
while(hasParent)&&(~strcmpi(iterresidue.residueType.name,'freeEnd'))
    iterresidue = iterresidue.linkageParent.parent;
    branchdepth = branchdepth + 1 ;
end

end

