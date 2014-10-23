function smallglypepstring = toGlycanListinSmallGlyPep(obj,removeisomer,varargin)

if(length(varargin)==1)
    outputfilename = varargin{1};
    fid = fopen(outputfilename,'w');    
end

outputglycan        = CellArrayList;
nonisomerglycanlist = CellArrayList; 
if(removeisomer)
    for i = 1 :obj.getNSpecies 
       glycancomp =  obj.theSpecies.get(i).glycanStruct.getComposition;
       if(~nonisomerglycanlist.contains(glycancomp))
           nonisomerglycanlist.add(glycancomp);
           outputglycan.add(obj.theSpecies.get(i)); 
       end
    end
else
    outputglycan = obj.theSpecies;    
end        

smallglypepstring    = '';
for i = 1 : outputglycan.length;
   glycanlinucs = outputglycan.get(i).glycanStruct.toLinucs;
   glycansmallglypep = sprintf('%s\n',linucs2SmallGlyPep(glycanlinucs,'linucs'));
   smallglypepstring = [smallglypepstring glycansmallglypep];
end

if(length(varargin)==1)
    fprintf(fid,'%s',smallglypepstring);    
end

end