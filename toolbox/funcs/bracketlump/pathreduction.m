function newpathway = pathreduction(glypath)
% PATHREDUCTION simplify the pathway using bracket lumping method
%
%  newpathway = pathreduction(glypath) reduces the network 
%     to a smaller one. 
% 
%See also: bracketSpeciesList. 

% Authro: Gang Liu
% Date Lastly Updated: 7/6/2014

newpathway = Pathway;
glypath.setListOfSpecies;
newgroupspecies = bracketspecieslist(glypath.listofSpecies);

indexmap = containers.Map;
for i = 1 : length(newgroupspecies)
    speciesindex = newgroupspecies(i,1).speciesindex;
    for j = 1 : length(speciesindex)
        indexmap(int2str(speciesindex(j)))=i;
    end
end

% build reaction matrix
rxnmatrix = zeros(length(glypath.theRxns),3);
for i = 1 : length(glypath.theRxns)
    ithrxn = glypath.theRxns.get(i);
    for j = 1 : length(glypath.theSpecies)
        if(ithrxn.reac==glypath.theSpecies.get(j))
            rxnmatrix(i,1)=indexmap(int2str(j));
        end
        
        if(ithrxn.prod==glypath.theSpecies.get(j))
            rxnmatrix(i,2)=indexmap(int2str(j));
        end
    end
    
    for k = 1 : length(glypath.theEnzs)
       if(ithrxn.enz==glypath.theEnzs.get(k))
           rxnmatrix(i,3)=k;
       end
    end
end

rxnmatrix = unique(rxnmatrix,'rows');

% add group species
for i = 1 : length(newgroupspecies)
    ithgroupspecies = newgroupspecies(i,1);
    if(ithgroupspecies.singlespec)
        newpathway.addGlycan(ithgroupspecies.glycanspecies);
    else
        newpathway.addGlycan(ithgroupspecies.bracketspecies);
    end
end

% add enzymes
newpathway.theEnzs = glypath.theEnzs;

% add rxns
for i = 1 : size(rxnmatrix,1);
    reacindex = rxnmatrix(i,1);
    prodindex = rxnmatrix(i,2);
    enzindex  = rxnmatrix(i,3);
    newpathway.addRxn(newpathway.theSpecies.get(reacindex),...
        newpathway.theSpecies.get(prodindex),newpathway.theEnzs.get(enzindex));
end

% set up graph structure
for i = 1 : size(newpathway.theSpecies)
    ithspecies = newpathway.theSpecies.get(i);
    for j = 1 : length(newpathway.theRxns)
        jthrxn = newpathway.theRxns.get(j);
        if(ithspecies==jthrxn.reac)
            ithspecies.listOfReacRxns.add(jthrxn);
        end
        
        if(ithspecies==jthrxn.prod)
            ithspecies.listOfProdRxns.add(jthrxn);
        end
    end
end

% add compartment info
if(~isempty(newpathway.compartment))
    for i = 1 : length(newpathway.theSpecies)
        newpathway.theSpecies.get(i).compartment=newpathway.compartment.get(1);
    end
end

end