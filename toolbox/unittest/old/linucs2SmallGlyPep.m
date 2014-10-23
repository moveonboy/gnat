function  SmallGlyPepnew = linucs2SmallGlyPep(pathwayinput)
theSpecies = pathwayinput.theSpecies;
for i = 1:length(theSpecies)
        intermediate{i} = theSpecies.get(i).glycanStruct.toLinucs;
end
for i = 1:length(intermediate)
SmallGlyPep = Linucs2SmallGly(intermediate{i});
capexpr='[a-z]';
 capexpr2='[{}]';
 capexpr3='{}';
 capexpr4='}[a-z]{';
 SmallGlyPeptemp=char(zeros(1,length(SmallGlyPep)));
 list=regexp(SmallGlyPep,capexpr);
 list2=regexp(SmallGlyPep,capexpr2);
 list3=regexp(SmallGlyPep,capexpr3);
 list4=regexp(SmallGlyPep,capexpr4);
 list4=[list4 length(SmallGlyPep)];
 j=0;
 for m=1:length(list3)*2
     if mod(m,2)
         tempgly=list(find(list <= list3((m+1)/2) & list > j));
         tempstr=list2(find(list2 <= list3((m+1)/2) & list2 > j));
         tempglyr=tempgly+ones(1,length(tempgly));
         tempstrr=tempstr-ones(1,length(tempstr));
         SmallGlyPeptemp(tempglyr) = SmallGlyPep(tempgly);
         SmallGlyPeptemp(tempstrr) = SmallGlyPep(tempstr);
         j=list3((m+1)/2);
     else
         SmallGlyPeptemp(list3(m/2)+1:list4(m/2))=SmallGlyPep(list3(m/2)+1:list4(m/2));
         j=list4(m/2);
     end
 end
 SmallGlyPepnew{i} = SmallGlyPeptemp;
end
SmallGlyPepnew = SmallGlyPepnew';
end