function [glycancompstring,varargout]=glyCompString(glycanStructObj)
%gly2CompString 
%
% See also glycanFormula,gly1charformat.

% Date Lastly Updated: 5/22/14
glycomp = glycanStructObj.getComposition;

% Hex
monoresidues = fieldnames(glycomp);
newglycomp=struct('HexNAc',0,'Hex',0,'Fuc',0,'NeuGc',0,...
                  'NeuAc',0,'Xyl',0);
for i = 1:length(monoresidues)
    singleresidue = monoresidues{i};
    switch lower(singleresidue)
       case {'glcnac','galnac'}
          newglycomp.HexNAc   = newglycomp.HexNAc+glycomp.(singleresidue);  
       case {'glc','gal','man'}
          newglycomp.Hex      = newglycomp.Hex+glycomp.(singleresidue);
       case  {'fuc'}
          newglycomp.Fuc      = newglycomp.Fuc+glycomp.(singleresidue);
       case {'neugc'}
           newglycomp.NeuGc   = newglycomp.NeuGc+glycomp.(singleresidue); 
       case {'neuac'}
            newglycomp.NeuAc  = newglycomp.NeuAc+glycomp.(singleresidue); 
       case {'xyl'}
            newglycomp.Xyl  = newglycomp.Xyl+glycomp.(singleresidue);      
    end
end

monoresidues = fieldnames(newglycomp);
glycancompstring = [];
for i = 1:length(monoresidues)
    if(newglycomp.(monoresidues{i})~=0)
      glycancompstring = [glycancompstring, monoresidues{i},...
        int2str(newglycomp.(monoresidues{i}))];    
    end
end

if(nargout==2)
    varargout{1}=newglycomp;
end    

end

