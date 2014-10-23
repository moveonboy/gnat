function isBranchSym = compareBranch(branches1,branches2)
% COMPAREBRANCH compare the branches in complex N-glycan 
%    isBranchSym = compareBranch(branch1,branch2)
%
% See also getBranches.

% Author: Gang Liu
% Date: 7/14/2014

if(branches1.numbranch~=branches2.numbranch)
    isBranchSym = 0;
    return
end

if(~isequal(sort(branches1.nonemptyfields),sort(branches2.nonemptyfields)))
    isBranchSym = 0;
    return
end

branchnames = branches1.nonemptyfields;

for i = 1 : length(branchnames)
   filedname1 = branchnames{i};
   branch1 =  branches1.(filedname1);
   ismatched = 0;
   for j = 1 : length(branchnames)
       if(j==i)
         continue
       end    
       
       fieldname2 = branchnames{j};
       branch2 = branches2.(fieldname2);
       if(strcmpi(branch1.residuestrwolink,branch2.residuestrwolink))
         ismatched = 1;
       end
   end
   
   if(~ismatched)
       isBranchSym = 0;
       return
   end   
end

isBranchSym = 1;

end