% test bracket species
storedvarpath = load('Bi_Nglycanpathway.mat');
binglycanpath = storedvarpath.nlinkedpath;
binglycanpath.resetjava;
glycanPathViewer(binglycanpath);

binglycanpath.setListOfSpecies;
newgroupspecies = bracketSpeciesList(binglycanpath.listofSpecies);

% for i = 1 : length(newgroupspecies)
%     ithgroup = newgroupspecies(i,1);
%     for j = 1 : length(ithgroup.glycanspecies)
%         glycanViewer(ithgroup.glycanspecies(j,1).glycanStruct)
%     end
%     if(~isempty(ithgroup.bracketspecies))
%         glycanViewer(ithgroup.bracketspecies.glycanStruct);
%         i
%     end
% end

for i = 1 : length(newgroupspecies)
    ithgroup = newgroupspecies(i,1);
    if(~isempty(ithgroup.bracketspecies))
        glycanViewer(ithgroup.bracketspecies.glycanStruct);
        i
    end
end
disp('end');