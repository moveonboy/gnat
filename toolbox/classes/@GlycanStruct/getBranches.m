function branches = getBranches(obj)
%getBranches get information about the branch
%
% See also getDepth

branches.numbranch = 0;
counter = 0;
nonredendresidues_species1 =  obj.getNonRedEndResidue;
for i = 1 : length(nonredendresidues_species1)
    ithresidue = nonredendresidues_species1{1,i};
    if(strcmpi(ithresidue.residueType.name,'Fuc'))
       continue; 
    end
    branches.numbranch      = branches.numbranch +1;    
    ithresidue              = nonredendresidues_species1{i};    
    branch.residuestrwolink = obj.getresiduetoroot(ithresidue,'nolinkage');    
    branch.depth            = obj.getBranchDepth(ithresidue);
    branch.residuestrwlink  = obj.getresiduetoroot(ithresidue);
    branch.resiudeterminal  = ithresidue;    
    
    if(strfind(branch.residuestrwlink,'3a1D-Man,p--2b1D-GlcNAc,p'))
        branches.M3G2=branch;
        branches.nonemptyfields{counter+1}='M3G2';
    elseif(strfind(branch.residuestrwlink,'3a1D-Man,p--4b1D-GlcNAc,p'))
        branches.M3G4=branch;
         branches.nonemptyfields{counter+1}='M3G4';
    elseif(strfind(branch.residuestrwlink,'6a1D-Man,p--2b1D-GlcNAc,p'))
        branches.M6G2=branch;
         branches.nonemptyfields{counter+1}='M6G2';
    elseif(strfind(branch.residuestrwlink,'6a1D-Man,p--6b1D-GlcNAc,p'))
        branches.M6G6=branch;
        branches.nonemptyfields{counter+1}='M6G6';
    else
        error('MATLAB:GNAT:ERRORBRANCH','INCORRECT BRANCH');
    end
    
    counter = counter +1;
end

end

