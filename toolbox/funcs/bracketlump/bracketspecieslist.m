function newgroupspecies = bracketspecieslist(listofSpecies)

% first step: group the glycans based on the composition
groupspecies    = classifyglycans(listofSpecies);
% second step: group the complex glycans based on the branch length
newgroupspecies = classifycomplexnglycan(groupspecies);
% third step: create the bracket based on list of glycans
newgroupspecies = buildbracketspecies(newgroupspecies);
end

function newgroupspecies = buildbracketspecies(newgroupspecies)
for i = 1 : length(newgroupspecies)
    if(~newgroupspecies(i,1).singlespec)
        glycanspecieslistinithgroup = newgroupspecies(i,1).glycanspecies;
        if(newgroupspecies(i,1).numbranch==2)
            [bracketspecies,corestruct]=createBiAntennaryBracket(glycanspecieslistinithgroup);
        elseif(newgroupspecies(i,1).numbranch==3)
            [bracketspecies,corestruct]=createTriAntennaryBracket(glycanspecieslistinithgroup);
        elseif(newgroupspecies(i,1).numbranch==4)
            [bracketspecies,corestruct]=createTetriAntennaryBracket(glycanspecieslistinithgroup);
        end
        newgroupspecies(i,1).bracketspecies = bracketspecies;
        newgroupspecies(i,1).corestruct     = corestruct;
    end
end
end

function [bracketspecies,corestruct] = createBiAntennaryBracket(glycanspecieslistinithgroup)
if(length(glycanspecieslistinithgroup)~=2)
    error('MATLAB:GNAT:ERRORLIST','THE NUMBER OF SPECIES SHOULD BE EQUALT TO 2');
end
glycanspecies1 = glycanspecieslistinithgroup(1,1);
glycanspecies2 = glycanspecieslistinithgroup(2,1);

nonredendresidues_species1 =  glycanspecies1.glycanStruct.getNonRedEndResidue;
for i = 1: length(nonredendresidues_species1)
    ithnredendres = nonredendresidues_species1{i};
    if(~strcmpi(ithnredendres.residueType.name,'Fuc'))
        branches_species1{i}   = glycanspecies1.glycanStruct.getresiduetoroot(ithnredendres);
        branchspecies1depth(i) = glycanspecies1.glycanStruct.getBranchDepth(ithnredendres);
    end
end

nonredendresidues_species2 =  glycanspecies2.glycanStruct.getNonRedEndResidue;
for i = 1: length(nonredendresidues_species2)
    ithnredendres = nonredendresidues_species2{i};
    branches_species2{i} = glycanspecies2.glycanStruct.getresiduetoroot(ithnredendres);
    branchspecies2depth(i) = glycanspecies2.glycanStruct.getBranchDepth(ithnredendres);
end

% check if two branches are the same

% take 1 species
branchdepthdiff = branchspecies1depth(1) - branchspecies1depth(2);
if(branchdepthdiff>0)
    longerbranchnum = 1;
    shorterbranchnum = 2;
else
    longerbranchnum = 2;
    shorterbranchnum = 1;
end

% create core structure
corestruct = createCoreStruct(glycanspecies1,longerbranchnum,shorterbranchnum,branchdepthdiff);

% create
bracketspecies = createBracketStruct(glycanspecies1,longerbranchnum,shorterbranchnum,branchdepthdiff);
end

function bracketspecies = createBracketStruct(glycanspecies1,longerbranchnum,shorterbranchnum,branchdepthdiff)
    bracketspecies = glycanspecies1.clone;
    nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
    longerbranchnonredendresidue = nonredendresidues_species1{longerbranchnum};
    shorterbranchnonredendresidue = nonredendresidues_species1{shorterbranchnum};

    for i = 1 : abs(branchdepthdiff)-1
        parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
        longerbranchnonredendresidue = parentresidue;
    end

    parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
    while(~strcmpi(parentresidue.residueType.name,shorterbranchnonredendresidue.residueType.name))
        parentresidue = parentresidue.getLinkageParent.getParent;
        shorterbranchnonredendresidue = shorterbranchnonredendresidue.getLinkageParent.getParent;
    end
    if(~isempty(shorterbranchnonredendresidue.getLinkageChildren))
        residuestoremove = [parentresidue.getLinkageChildren.getChild;...
            shorterbranchnonredendresidue.getLinkageChildren.getChild];
    else
        residuestoremove = parentresidue.getLinkageChildren.getChild;
    end
    
    bracketspecies.glycanStruct.bracketResidue(residuestoremove);
end

function corestruct = createCoreStruct(glycanspecies1,longerbranchnum,shorterbranchnum,branchdepthdiff)
    corestruct = glycanspecies1.clone;
    nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
    longerbranchnonredendresidue = nonredendresidues_species1{longerbranchnum};
    shorterbranchnonredendresidue = nonredendresidues_species1{shorterbranchnum};

    for i = 1 : abs(branchdepthdiff)-1
        parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
        longerbranchnonredendresidue = parentresidue;
    end

    parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
    while(~strcmpi(parentresidue.residueType.name,shorterbranchnonredendresidue.residueType.name))
        parentresidue = parentresidue.getLinkageParent.getParent;
        shorterbranchnonredendresidue = shorterbranchnonredendresidue.getLinkageParent.getParent;
    end
    if(~isempty(shorterbranchnonredendresidue.getLinkageChildren))
        residuestoremove = [parentresidue.getLinkageChildren.getChild;...
            shorterbranchnonredendresidue.getLinkageChildren.getChild];
    else
        residuestoremove = parentresidue.getLinkageChildren.getChild;
    end
    for i = 1 : length(residuestoremove)
        corestruct.glycanStruct.removeResidueObj(residuestoremove(i));
    end   
end


function [bracketspecies,corestruct] = createTriAntennaryBracket(glycanspecieslistinithgroup)

glycanspecies1 = glycanspecieslistinithgroup(1,1);
nonredendresidues_species1 =  glycanspecies1.glycanStruct.getNonRedEndResidue;
nonredendresidues_species1 = rmfucose(nonredendresidues_species1);
[species1branchindex,species1branchdepth,~] = listBranch(glycanspecies1,nonredendresidues_species1);

A = species1branchdepth(1);
B = species1branchdepth(2);
C = species1branchdepth(3);
if(A==B)&&(B==C)
    % create core structure
    corestruct = glycanspecies1.clone;
    nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
    
    for i = 1 : length(nonredendresidues_species1)
        ithresidue = nonredendresidues_species1{i};
        residuestoremove = ithresidue;
        corestruct.glycanStruct.removeResidueObj(residuestoremove);
    end
    
    % create bracketspecies
    bracketspecies = glycanspecies1.clone;
    nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
    
    for i = 1 : length(nonredendresidues_species1)
        ithresidue = nonredendresidues_species1{i};
        residuestobracket(i,1) = ithresidue;
    end
    bracketspecies.glycanStruct.bracketResidue(residuestobracket);
    
elseif(A~=B)&&(B~=C)&&(A~=C)
    [maxbranch,indexmax] = max(species1branchdepth);
    [minbranch,indexmin] = min(species1branchdepth);
    if(isequal(nonredendresidues_species1{indexmin}.residueType.name,'NeuAc'))
        for i = 1 : length(species1branchdepth)
            if(i~=indexmax)&&(i~=indexmin)
                branchdepthdiff = species1branchdepth(i)-minbranch;
                
                % create core structure
                corestruct = glycanspecies1.clone;
                nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
                longerbranchnonredendresidue = nonredendresidues_species1{i};
                for j = 1 : abs(branchdepthdiff)
                    parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                    longerbranchnonredendresidue = parentresidue;
                end
                residuestoremove = longerbranchnonredendresidue;
                corestruct.glycanStruct.removeResidueObj(residuestoremove);
                
                % create bracketspecies
                bracketspecies = glycanspecies1.clone;
                nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
                longerbranchnonredendresidue = nonredendresidues_species1{i};
                for j = 1 : abs(branchdepthdiff)
                    parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                    longerbranchnonredendresidue = parentresidue;
                end
                residuestobracket(1,1) = longerbranchnonredendresidue;
            end
        end
        
        branchdepthdiff = maxbranch-minbranch;
        
        % create core structure
        nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
        longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
        
        for i = 1 : abs(branchdepthdiff)
            parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
            longerbranchnonredendresidue = parentresidue;
        end
        
        residuestoremove = longerbranchnonredendresidue;
        corestruct.glycanStruct.removeResidueObj(residuestoremove);
        
        % create bracketspecies
        nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
        longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
        
        for i = 1 : abs(branchdepthdiff)
            parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
            longerbranchnonredendresidue = parentresidue;
        end
        
        residuestobracket(2,1) = longerbranchnonredendresidue;
        
        % create core structure
        nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
        Sianonredendresidue = nonredendresidues_species1{indexmin};
        residuestoremove = Sianonredendresidue;
        corestruct.glycanStruct.removeResidueObj(residuestoremove);
        
        % create bracketspecies
        nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
        longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
        
        for i = 1 : abs(branchdepthdiff)
            parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
            longerbranchnonredendresidue = parentresidue;
        end
        
        residuestobracket(3,1) = longerbranchnonredendresidue;
        bracketspecies.glycanStruct.bracketResidue(residuestobracket);
    else
        for i = 1 : length(species1branchdepth)
            if(i~=indexmax)&&(i~=indexmin)
                branchdepthdiff = species1branchdepth(i)-minbranch;
                
                % create core structure
                corestruct = glycanspecies1.clone;
                nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
                longerbranchnonredendresidue = nonredendresidues_species1{i};
                
                for j = 1 : abs(branchdepthdiff)-1
                    parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                    longerbranchnonredendresidue = parentresidue;
                end
                
                residuestoremove = longerbranchnonredendresidue;
                corestruct.glycanStruct.removeResidueObj(residuestoremove);
                
                % create bracketspecies
                bracketspecies = glycanspecies1.clone;
                nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
                longerbranchnonredendresidue = nonredendresidues_species1{i};
                
                for j = 1 : abs(branchdepthdiff)-1
                    parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                    longerbranchnonredendresidue = parentresidue;
                end
                
                residuestobracket(1,1) = longerbranchnonredendresidue;
            end
        end
        
        branchdepthdiff = maxbranch-minbranch;
        
        % create core structure
        nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
        longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
        
        for i = 1 : abs(branchdepthdiff)-1
            parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
            longerbranchnonredendresidue = parentresidue;
        end
        
        residuestoremove = longerbranchnonredendresidue;
        corestruct.glycanStruct.removeResidueObj(residuestoremove);
        
        % create bracketspecies
        nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
        longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
        
        for i = 1 : abs(branchdepthdiff)-1
            parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
            longerbranchnonredendresidue = parentresidue;
        end
        
        residuestobracket(2,1) = longerbranchnonredendresidue;
        bracketspecies.glycanStruct.bracketResidue(residuestobracket);
    end
    
else
    [maxbranch,indexmax] = max(species1branchdepth);
    [minbranch,indexmin] = min(species1branchdepth);
    indexmedian          = [];
    
    for i = 1 : length(species1branchdepth)
        if(i==indexmax)
            continue
        end
        if(i==indexmin)
            continue
        end
        indexmedian = i;
    end
    
    if(isequal(species1branchdepth(indexmedian),maxbranch))
        if(isequal(nonredendresidues_species1{indexmin}.residueType.name,'NeuAc'))
            branchdepthdiff = maxbranch-minbranch;
            
            % create core structure
            corestruct = glycanspecies1.clone;
            nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue1 = nonredendresidues_species1{indexmedian};
            longerbranchnonredendresidue2 = nonredendresidues_species1{indexmax};
            shortbranchnonredendresidue   = nonredendresidues_species1{indexmin};
            
            for i = 1 : abs(branchdepthdiff)
                parentresidue1 = longerbranchnonredendresidue1.getLinkageParent.getParent;
                longerbranchnonredendresidue1 = parentresidue1;
                parentresidue2 = longerbranchnonredendresidue2.getLinkageParent.getParent;
                longerbranchnonredendresidue2 = parentresidue2;
            end
            
            residuestoremove1 = longerbranchnonredendresidue1;
            residuestoremove2 = longerbranchnonredendresidue2;
            residuestoremove3 = shortbranchnonredendresidue;
            corestruct.glycanStruct.removeResidueObj(residuestoremove1);
            corestruct.glycanStruct.removeResidueObj(residuestoremove2);
            corestruct.glycanStruct.removeResidueObj(residuestoremove3);
            
            % create bracketspecies
            bracketspecies = glycanspecies1.clone;
            nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue1 = nonredendresidues_species1{indexmedian};
            longerbranchnonredendresidue2 = nonredendresidues_species1{indexmax};
            shortbranchnonredendresidue   = nonredendresidues_species1{indexmin};
            
            for i = 1 : abs(branchdepthdiff)
                parentresidue1 = longerbranchnonredendresidue1.getLinkageParent.getParent;
                longerbranchnonredendresidue1 = parentresidue1;
                parentresidue2 = longerbranchnonredendresidue2.getLinkageParent.getParent;
                longerbranchnonredendresidue2 = parentresidue2;
            end
            
            residuestobracket(1,1) = longerbranchnonredendresidue1;
            residuestobracket(2,1) = longerbranchnonredendresidue2;
            residuestobracket(3,1) = shortbranchnonredendresidue;
            bracketspecies.glycanStruct.bracketResidue(residuestobracket);
        else
            branchdepthdiff = maxbranch-minbranch;
            
            % create core structure
            corestruct = glycanspecies1.clone;
            nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue1 = nonredendresidues_species1{indexmedian};
            longerbranchnonredendresidue2 = nonredendresidues_species1{indexmax};
            
            for i = 1 : abs(branchdepthdiff)-1
                parentresidue1 = longerbranchnonredendresidue1.getLinkageParent.getParent;
                longerbranchnonredendresidue1 = parentresidue1;
                parentresidue2 = longerbranchnonredendresidue2.getLinkageParent.getParent;
                longerbranchnonredendresidue2 = parentresidue2;
            end
            
            residuestoremove1 = longerbranchnonredendresidue1;
            residuestoremove2 = longerbranchnonredendresidue2;
            corestruct.glycanStruct.removeResidueObj(residuestoremove1);
            corestruct.glycanStruct.removeResidueObj(residuestoremove2);
            
            % create bracketspecies
            bracketspecies = glycanspecies1.clone;
            nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue1 = nonredendresidues_species1{indexmedian};
            longerbranchnonredendresidue2 = nonredendresidues_species1{indexmax};
            
            for i = 1 : abs(branchdepthdiff)-1
                parentresidue1 = longerbranchnonredendresidue1.getLinkageParent.getParent;
                longerbranchnonredendresidue1 = parentresidue1;
                parentresidue2 = longerbranchnonredendresidue2.getLinkageParent.getParent;
                longerbranchnonredendresidue2 = parentresidue2;
            end
            
            residuestobracket(1,1) = longerbranchnonredendresidue1;
            residuestobracket(2,1) = longerbranchnonredendresidue2;
            bracketspecies.glycanStruct.bracketResidue(residuestobracket);
        end
    else
        branchdepthdiff = maxbranch-minbranch;
        isNeuAc = 0;
        for i = 1 : length(nonredendresidues_species1)
            if(i==indexmax)
                continue
            end
            if(isequal(nonredendresidues_species1{i}.residueType.name,'NeuAc'))
                isNeuAc = isNeuAc+1;
            end
        end
        
        if(isNeuAc)
            % create core structure
            corestruct = glycanspecies1.clone;
            nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
            shortbranchnonredendresidue1   = nonredendresidues_species1{indexmin};
            
            for i = 1 : length(species1branchdepth)
                if(i==indexmin)||(i==indexmax)
                    continue
                end
                shortbranchnonredendresidue2 = nonredendresidues_species1{i};
            end
            
            for i = 1 : abs(branchdepthdiff)
                parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                longerbranchnonredendresidue = parentresidue;
            end
            
            residuestoremove1 = longerbranchnonredendresidue;
            residuestoremove2 = shortbranchnonredendresidue1;
            residuestoremove3 = shortbranchnonredendresidue2;
            corestruct.glycanStruct.removeResidueObj(residuestoremove1);
            corestruct.glycanStruct.removeResidueObj(residuestoremove2);
            corestruct.glycanStruct.removeResidueObj(residuestoremove3);
            
            % create bracketspecies
            bracketspecies = glycanspecies1.clone;
            nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
            shortbranchnonredendresidue1   = nonredendresidues_species1{indexmin};
            
            for i = 1 : length(species1branchdepth)
                if(i==indexmin)||(i==indexmax)
                    continue
                end
                shortbranchnonredendresidue2 = nonredendresidues_species1{i};
            end
            
            for i = 1 : abs(branchdepthdiff)
                parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                longerbranchnonredendresidue = parentresidue;
            end
            
            residuestobracket(1,1) = longerbranchnonredendresidue;
            residuestobracket(2,1) = shortbranchnonredendresidue1;
            residuestobracket(3,1) = shortbranchnonredendresidue2;
            bracketspecies.glycanStruct.bracketResidue(residuestobracket);
        else
            % create core structure
            corestruct = glycanspecies1.clone;
            nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
            
            for i = 1 : abs(branchdepthdiff)-1
                parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                longerbranchnonredendresidue = parentresidue;
            end
            
            residuestoremove = longerbranchnonredendresidue;
            corestruct.glycanStruct.removeResidueObj(residuestoremove);
            
            % create bracketspecies
            bracketspecies = glycanspecies1.clone;
            nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
            longerbranchnonredendresidue = nonredendresidues_species1{indexmax};
            
            for i = 1 : abs(branchdepthdiff)-1
                parentresidue = longerbranchnonredendresidue.getLinkageParent.getParent;
                longerbranchnonredendresidue = parentresidue;
            end
            
            residuestobracket(1,1) = longerbranchnonredendresidue;
            bracketspecies.glycanStruct.bracketResidue(residuestobracket);
        end
    end
end
end

function [bracketspecies,corestruct] = createTetriAntennaryBracket(glycanspecieslistinithgroup)

glycanspecies1 = glycanspecieslistinithgroup(1,1);
corestruct = glycanspecies1.clone;
nonredendresidues_species1 = glycanspecies1.glycanStruct.getNonRedEndResidue;
nonredendresidues_species1 = rmfucose(nonredendresidues_species1);
[species1branchindex,species1branchdepth,branchstr] = listBranch(glycanspecies1,nonredendresidues_species1);

% get the minimum branch
[minbranchdepth,branchindex] = min(species1branchdepth);
minbranchstr = branchstr{branchindex};

% get other branch string
otherbranchdepth = [];
otherbranch      = cell(3,1);
counter          = 1;
for i = 1 : length(branchstr)
    if(i==branchindex)
        continue
    end
    ithbranchstr              = branchstr{i};
    ithbranchdepth            = species1branchdepth(i);
    otherbranch(counter)      = cellstr(ithbranchstr);
    otherbranchdepth(counter,1) = ithbranchdepth;
    counter = counter+1;
end

% find the common branch depth
isvalid = 0;
for i = 1 : length(otherbranch)
    ithbranch = otherbranch{i};
    if(strfind(ithbranch,minbranchstr))
        isvalid = isvalid+1;
    end
end

if(isvalid==3)
    % create core structure
    corestruct = glycanspecies1.clone;
    nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
    
    counter = 1;
    for i = 1 : length(branchstr)
        if(i==branchindex)
            continue
        end
        ithresidue                = nonredendresidues_species1{i};
        nresidues(counter,1)        = ithresidue;
        counter = counter+1;
    end
    
    for i = 1 : length(nresidues)
        ithresidue = nresidues(i,1);
        branchdepthdiff = otherbranchdepth(i,1)-minbranchdepth;
        if(branchdepthdiff==0)
            continue
        end
        for j = 1 : abs(branchdepthdiff)-1
            parentresidue = ithresidue.getLinkageParent.getParent;
            ithresidue = parentresidue;
        end
        residuestoremove = ithresidue;
        corestruct.glycanStruct.removeResidueObj(residuestoremove);
    end
    
    % create bracketspecies
    bracketspecies    = glycanspecies1.clone;
    nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
    
    counter = 1;
    for i = 1 : length(branchstr)
        if(i==branchindex)
            continue
        end
        ithresidue                = nonredendresidues_species1{i};
        nresidues(counter,1)        = ithresidue;
        counter = counter+1;
    end
    
    counter = 1;
    for i = 1 : length(nresidues)
        ithresidue = nresidues(i,1);
        branchdepthdiff = otherbranchdepth(i,1)-minbranchdepth;
        if(branchdepthdiff==0)
            continue
        end
        for j = 1 : abs(branchdepthdiff)-1
            parentresidue = ithresidue.getLinkageParent.getParent;
            ithresidue = parentresidue;
        end
        residuestobracket(counter,1) = ithresidue;
        counter                      = counter+1;
    end
    bracketspecies.glycanStruct.bracketResidue(residuestobracket);
else
    % move on to backward residue
    ismatch = 0;
    for i = 1 : length(minbranchdepth)
        minresidue     = nonredendresidues_species1{branchindex};
        backresidue    = minresidue.getLinkageParent.getParent;
        minbranchstr   = glycanSpecies1.glycanStruct.getresiduetoroot(backresidue);
        minbranchdepth = minbranchdepth-1;
        for j = 1 : length(otherbranch)
            jthbranch = otherbranch{j};
            if(strfind(jthbranch,minbranchstr))
                ismatch = ismatch+1;
            end
        end
        if(ismatch==3)
            break;
        else
            continue;
        end
    end
    
    % create core structure
    corestruct = glycanspecies1.clone;
    nonredendresidues_species1 =  corestruct.glycanStruct.getNonRedEndResidue;
    
    counter = 1;
    for i = 1 : length(branchstr)
        if(i==branchindex)
            continue
        end
        ithresidue                = nonredendresidues_species1{i};
        nresidues(counter,1)        = ithresidue;
        counter = counter+1;
    end
    
    for i = 1 : length(nresidues)
        ithresidue = nresidues(i,1);
        branchdepthdiff = otherbranchdepth(i,1)-minbranchdepth;
        if(branchdepthdiff==0)
            residuestoremove = ithresidue;
            corestruct.glycanStruct.removeResidueObj(residuestoremove);
        else
            for j = 1 : abs(branchdepthdiff)-1
                parentresidue = ithresidue.getLinkageParent.getParent;
                ithresidue = parentresidue;
            end
            residuestoremove = ithresidue;
            corestruct.glycanStruct.removeResidueObj(residuestoremove);
        end
    end
    corestruct.glycanStruct.removeResidueObj(minresidue);
    
    % create bracketspecies
    bracketspecies    = glycanspecies1.clone;
    nonredendresidues_species1 =  bracketspecies.glycanStruct.getNonRedEndResidue;
    
    counter = 1;
    for i = 1 : length(branchstr)
        if(i==branchindex)
            continue
        end
        ithresidue                = nonredendresidues_species1{i};
        nresidues(counter,1)        = ithresidue;
        counter = counter+1;
    end
    
    counter           = 1;
    for i = 1 : length(nresidues)
        ithresidue = nresidues(i,1);
        branchdepthdiff = otherbranchdepth(i,1)-minbranchdepth;
        if(branchdepthdiff==0)
            index(end+1) = i;
        else
            for j = 1 : abs(branchdepthdiff)-1
                parentresidue = ithresidue.getLinkageParent.getParent;
                ithresidue = parentresidue;
            end
            residuestobracket(counter,1) = ithresidue;
            counter                      = counter+1;
        end
    end
    residuestobracket(end+1)           = minresidue;
    for i = 1 : length(index)
        ithindex = index(i);
        residuestobracket(end+1)           = nresidues(ithindex,1);
    end
    bracketspecies.glycanStruct.bracketResidue(residuestobracket);
end

end

function newgroupspecies = classifycomplexnglycan(groupspecies)
newcounter    = 0;
newgroupspecies  = struct('composition',[],'glycanspecies',[],'spectype',[],'numbranch',[],...
    'corestruct',[],'bracketspecies',[],'singlespec',[],'glycandepth',[],'speciesindex',[]);
for i = 1 : length(groupspecies)
    if(~groupspecies(i,1).singlespec)
        subsetgroupspecies = createSusbsetGroups(groupspecies(i,1));
        for j = 1 : length(subsetgroupspecies)
            newcounter = newcounter +1;
            newgroupspecies(newcounter,1)=subsetgroupspecies(j,1);
        end
    else
        newcounter = newcounter +1;
        newgroupspecies(newcounter,1)=groupspecies(i,1);
    end
end
end

function subsetgroupspecies = createSusbsetGroups(ithgroupspecies)
subsetgroupspecies = struct('composition',[],'glycanspecies',[],'spectype',[],'numbranch',[],...
    'corestruct',[],'bracketspecies',[],'singlespec',[],'glycandepth',[],'speciesindex',[]);
counter = 0;
for i =1  : length(ithgroupspecies.glycanspecies)
    ithspecies        =  ithgroupspecies.glycanspecies(i,1);
    complexglycantype =  ithspecies.glycanStruct.iscomplex;
    if(~complexglycantype.complex)
        error('MATLAB:GNAT:ERRORNONCOMPLEX','NONCOMPLEX GLYCAN SELECTED');
    end
    isnewgroup = 1;
    numbranches = complexglycantype.branchnum;
    glycandepth = ithspecies.glycanStruct.getDepth;
    ithspeciesbranches    = ithspecies.glycanStruct.getBranches;
    
    if(numbranches == 2) % biantennary
        for ii = 1 : length(subsetgroupspecies)
            if(subsetgroupspecies(ii,1).glycandepth==glycandepth)
                subsetspeices = subsetgroupspecies(ii,1).glycanspecies(1,1);
                groupspeciesbranches = subsetspeices.glycanStruct.getBranches;
                if(compareBranch(ithspeciesbranches,groupspeciesbranches))
                    subsetgroupspecies(ii,1).glycanspecies = ...
                        [subsetgroupspecies(ii,1).glycanspecies;...
                        ithspecies];
                    subsetgroupspecies(ii,1).speciesindex = ...
                        [subsetgroupspecies(ii,1).speciesindex;...
                        ithgroupspecies.speciesindex(i)];
                    subsetgroupspecies(ii,1).singlespec = 0;
                    isnewgroup = 0;
                    break;                    
                end
            end
        end
        
        if(isnewgroup)
            counter=counter+1;
            subsetgroupspecies(counter,1).numbranch     = numbranches;
            subsetgroupspecies(counter,1).speciesindex  = ithgroupspecies.speciesindex(i);
            subsetgroupspecies(counter,1).glycandepth   = glycandepth;
            subsetgroupspecies(counter,1).glycanspecies = ithspecies;
            subsetgroupspecies(counter,1).singlespec    = 1;
            subsetgroupspecies(counter,1).composition   = ithgroupspecies.composition;
            subsetgroupspecies(counter,1).spectype      = ithgroupspecies.spectype;
        end
        
    elseif(numbranches == 3) % triantennary
        ithnreResidues = ithspecies.glycanStruct.getNonRedEndResidue;
        
        % get rid of Fuc from ithnreResidue list
        ithnreResidues = rmfucose(ithnreResidues);
        
        ithnreResiduestruct = struct('Gal',0,'GlcNAc',0,'NeuAc',0);
        for j = 1 : length(ithnreResidues)
            jthresidue = ithnreResidues{j};
            if(isequal(jthresidue.residueType.name,'Gal'))
                ithnreResiduestruct.Gal = ithnreResiduestruct.Gal+1;
            elseif(isequal(jthresidue.residueType.name,'GlcNAc'))
                ithnreResiduestruct.GlcNAc = ithnreResiduestruct.GlcNAc+1;
            elseif(isequal(jthresidue.residueType.name,'NeuAc'))
                ithnreResiduestruct.NeuAc = ithnreResiduestruct.NeuAc+1;
            end
        end
        % List each branch and get branch depth
        [ithbranchindex,ithbranchdepth,branchstr] = listBranch(ithspecies,ithnreResidues);
        
        for j = 1 : length(subsetgroupspecies)
            if(~isempty(subsetgroupspecies(j,1).glycanspecies))
                jthspecies = subsetgroupspecies(j,1).glycanspecies;
                jthnreResidues = jthspecies(1,1).glycanStruct.getNonRedEndResidue;
                % get rid of Fuc from ithnreResidue list
                jthnreResidues = rmfucose(jthnreResidues);
                % List each branch and get branch depth
                [jthbranchindex,jthbranchdepth,branchstr] = listBranch(jthspecies(1,1),jthnreResidues);
                
                isindexmatch   = 0;
                isdepthmatch   = 0;
                
                for jj = 1 : length(ithbranchindex)
                    jjthindex = ithbranchindex(jj);
                    if(strfind(num2str(jthbranchindex),num2str(jjthindex)))
                        isindexmatch = isindexmatch+1;
                    end
                end
                
                for jj = 1 : length(ithbranchdepth)
                    jjthdepth = ithbranchdepth(jj);
                    if(strfind(num2str(jthbranchdepth),num2str(jjthdepth)))
                        isdepthmatch = isdepthmatch+1;
                    end
                end
                
                jthnreResiduestruct = struct('Gal',0,'GlcNAc',0,'NeuAc',0);
                for jj = 1 : length(jthnreResidues)
                    jjthresidue = jthnreResidues{jj};
                    if(isequal(jjthresidue.residueType.name,'Gal'))
                        jthnreResiduestruct.Gal = jthnreResiduestruct.Gal+1;
                    elseif(isequal(jjthresidue.residueType.name,'GlcNAc'))
                        jthnreResiduestruct.GlcNAc = jthnreResiduestruct.GlcNAc+1;
                    elseif(isequal(jjthresidue.residueType.name,'NeuAc'))
                        jthnreResiduestruct.NeuAc = jthnreResiduestruct.NeuAc+1;
                    end
                end
                
                isresiduematch = 0;
                if(isequal(jthnreResiduestruct,ithnreResiduestruct))
                    isresiduematch = 1;
                end
                
                if(isindexmatch==3)&&(isdepthmatch==3)&&(isresiduematch)
                    subsetgroupspecies(j,1).glycanspecies = [subsetgroupspecies(j,1).glycanspecies;...
                        ithspecies];
                    subsetgroupspecies(j,1).speciesindex = [subsetgroupspecies(j,1).speciesindex;...
                        ithgroupspecies.speciesindex(i)];
                    subsetgroupspecies(j,1).singlespec = 0;
                    isnewgroup = 0;
                    break
                end
            end
        end
        
        if(isnewgroup)
            counter=counter+1;
            subsetgroupspecies(counter,1).numbranch     = numbranches;
            subsetgroupspecies(counter,1).speciesindex  = ithgroupspecies.speciesindex(i);
            subsetgroupspecies(counter,1).glycandepth   = glycandepth;
            subsetgroupspecies(counter,1).glycanspecies = ithspecies;
            subsetgroupspecies(counter,1).singlespec    = 1;
            subsetgroupspecies(counter,1).composition   = ithgroupspecies.composition;
            subsetgroupspecies(counter,1).spectype      = ithgroupspecies.spectype;
        end
        
    elseif(numbranches == 4) % tetraantennary
        ithnreResidues = ithspecies.glycanStruct.getNonRedEndResidue;
        
        % get rid of Fuc from ithnreResidue list
        ithnreResidues = rmfucose(ithnreResidues);
        
        ithnreResiduestruct = struct('Gal',0,'GlcNAc',0,'NeuAc',0);
        for j = 1 : length(ithnreResidues)
            jthresidue = ithnreResidues{j};
            if(isequal(jthresidue.residueType.name,'Gal'))
                ithnreResiduestruct.Gal = ithnreResiduestruct.Gal+1;
            elseif(isequal(jthresidue.residueType.name,'GlcNAc'))
                ithnreResiduestruct.GlcNAc = ithnreResiduestruct.GlcNAc+1;
            elseif(isequal(jthresidue.residueType.name,'NeuAc'))
                ithnreResiduestruct.NeuAc = ithnreResiduestruct.NeuAc+1;
            end
        end
        % List each branch and get branch depth
        [ithbranchindex,ithbranchdepth,ithbranchstr] = listBranch(ithspecies,ithnreResidues);
        
        for j = 1 : length(subsetgroupspecies)
            if(~isempty(subsetgroupspecies(j,1).glycanspecies))
                jthspecies = subsetgroupspecies(j,1).glycanspecies;
                jthnreResidues = jthspecies(1,1).glycanStruct.getNonRedEndResidue;
                % get rid of Fuc from ithnreResidue list
                jthnreResidues = rmfucose(jthnreResidues);
                % List each branch and get branch depth
                [jthbranchindex,jthbranchdepth,jthbranchstr] = listBranch(jthspecies(1,1),jthnreResidues);
                
                isdepthmatch   = 0;
                
                for jj = 1 : length(ithbranchdepth)
                    jjthdepth = ithbranchdepth(jj);
                    if(strfind(num2str(jthbranchdepth),num2str(jjthdepth)))
                        isdepthmatch = isdepthmatch+1;
                    end
                end
                
                jthnreResiduestruct = struct('Gal',0,'GlcNAc',0,'NeuAc',0);
                for jj = 1 : length(jthnreResidues)
                    jjthresidue = jthnreResidues{jj};
                    if(isequal(jjthresidue.residueType.name,'Gal'))
                        jthnreResiduestruct.Gal = jthnreResiduestruct.Gal+1;
                    elseif(isequal(jjthresidue.residueType.name,'GlcNAc'))
                        jthnreResiduestruct.GlcNAc = jthnreResiduestruct.GlcNAc+1;
                    elseif(isequal(jjthresidue.residueType.name,'NeuAc'))
                        jthnreResiduestruct.NeuAc = jthnreResiduestruct.NeuAc+1;
                    end
                end
                
                isresiduematch = 0;
                if(isequal(jthnreResiduestruct,ithnreResiduestruct))
                    isresiduematch = 1;
                end
                
                if(isdepthmatch==4)&&(isresiduematch)
                    subsetgroupspecies(j,1).glycanspecies = [subsetgroupspecies(j,1).glycanspecies;...
                        ithspecies];
                    subsetgroupspecies(j,1).speciesindex = [subsetgroupspecies(j,1).speciesindex;...
                        ithgroupspecies.speciesindex(i)];
                    subsetgroupspecies(j,1).singlespec = 0;
                    isnewgroup = 0;
                    break
                end
            end
        end
        if(isnewgroup)
            counter=counter+1;
            subsetgroupspecies(counter,1).numbranch     = numbranches;
            subsetgroupspecies(counter,1).speciesindex  = ithgroupspecies.speciesindex(i);
            subsetgroupspecies(counter,1).glycandepth   = glycandepth;
            subsetgroupspecies(counter,1).glycanspecies = ithspecies;
            subsetgroupspecies(counter,1).singlespec    = 1;
            subsetgroupspecies(counter,1).composition   = ithgroupspecies.composition;
            subsetgroupspecies(counter,1).spectype      = ithgroupspecies.spectype;
        end
    end
end
end

function groupspecies = classifyglycans(listofSpecies)
groupspecies  = struct('composition',[],'glycanspecies',[],'spectype',[],'numbranch',[],...
    'corestruct',[],'bracketspecies',[],'singlespec',[],'glycandepth',[],'speciesindex',[]);
counter = 0;
for i = 1 : length(listofSpecies)
    ithspecies     = listofSpecies{i,1};
    ithspeciescomp = ithspecies.glycanStruct.getComposition;
    
    % high mannose, each struct is considered a single group
    if(ithspecies.glycanStruct.ishighmannose)
        counter = counter +1;
        groupspecies(counter,1).glycanspecies = ithspecies;
        groupspecies(counter,1).composition   = ithspeciescomp;
        groupspecies(counter,1).singlespec    = 1;
        groupspecies(counter,1).spectype      = 'highman';
        groupspecies(counter,1).speciesindex  = i;
        continue;
    end
    
    % hybrid, each struct is considered as a single group
    if(ithspecies.glycanStruct.ishybrid)
        counter = counter +1;
        groupspecies(counter,1).glycanspecies  = ithspecies;
        groupspecies(counter,1).composition    = ithspeciescomp;
        groupspecies(counter,1).bracketspecies = [];
        groupspecies(counter,1).singlespec     = 1;
        groupspecies(counter,1).spectype       = 'hybrid';
        groupspecies(counter,1).speciesindex   = i;
        continue;
    end
    
    complextype = ithspecies.glycanStruct.iscomplex;
    if(complextype.complex==-1)
        counter = counter +1;
        groupspecies(counter,1).glycanspecies  = ithspecies;
        groupspecies(counter,1).composition    = ithspeciescomp;
        groupspecies(counter,1).bracketspecies = [];
        groupspecies(counter,1).singlespec     = 1;
        groupspecies(counter,1).spectype       = 'complex precursor';
        groupspecies(counter,1).numbranch      = complextype.branchnum;
        groupspecies(counter,1).speciesindex  = i;
        continue;
    end
    
    % if the composition is the same
    hasgroup = 0;
    for j = 1 : length(groupspecies)
        if(strcmpi(groupspecies(j,1).spectype,'complex precursor')|| ...
                strcmpi(groupspecies(j,1).spectype,'hybrid')|| ...
                strcmpi(groupspecies(j,1).spectype,'highman'))
            continue;
        end
        
        if(isequal(groupspecies(j,1).composition,ithspeciescomp))
            groupspecies(j,1).glycanspecies = [groupspecies(j,1).glycanspecies;ithspecies];
            groupspecies(j,1).glycandepth   = [groupspecies(j,1).glycandepth;ithspecies.glycanStruct.getDepth];
            groupspecies(j,1).speciesindex  = [groupspecies(j,1).speciesindex;i];
            groupspecies(j,1).singlespec    = 0;
            hasgroup = 1;
        end
    end
    
    if(~hasgroup)
        counter = counter +1;
        groupspecies(counter,1).glycanspecies  = ithspecies;
        groupspecies(counter,1).glycandepth    = ithspecies.glycanStruct.getDepth;
        groupspecies(counter,1).composition    = ithspeciescomp;
        groupspecies(counter,1).bracketspecies = [];
        groupspecies(counter,1).singlespec     = 1;
        groupspecies(counter,1).spectype       = 'complex';
        groupspecies(counter,1).numbranch      = complextype.branchnum;
        groupspecies(counter,1).speciesindex  =  i;
    end
end
end

function obj = rmfucose(obj)
for i = 1:length(obj)
    ithresidue = obj{i};
    if(isequal(ithresidue.residueType.name,'Fuc'))
        obj(i) = '';
    else
        continue
    end
end
end

function [branchindex,branchdepth,branchstr] = listBranch(glycanSpecies,nreResidues)
branchindex = [];
branchdepth = [];
branchstr   = cell(4,1);
counter     = 1;

for i = 1 : length(nreResidues)
    ithresidue          = nreResidues{i};
    ithbranchstr        = glycanSpecies.glycanStruct.getresiduetoroot(ithresidue);
   % ithbranchresiduestr = getonlyresiduetoroot(ithresidue);
    ithbranchresiduestr = glycanSpecies.glycanStruct.getresiduetoroot(ithresidue,'nolinkage');
    if(strfind(ithbranchstr,'3a1D-Man,p--2b1D-GlcNAc,p'))
        branchindex(counter) = 1;
    elseif(strfind(ithbranchstr,'3a1D-Man,p--4b1D-GlcNAc,p'))
        branchindex(counter) = 2;
    elseif(strfind(ithbranchstr,'6a1D-Man,p--2b1D-GlcNAc,p'))
        branchindex(counter) = 3;
    elseif(strfind(ithbranchstr,'6a1D-Man,p--6b1D-GlcNAc,p'))
        branchindex(counter) = 4;
    else
        error('MATLAB:GNAT:ERRORBRANCH','INCORRECT BRANCH');
    end
    
    ithbranchdepth         = glycanSpecies.glycanStruct.getBranchDepth(ithresidue);
    branchdepth(counter)   = ithbranchdepth;
    branchstr(counter)     = cellstr(ithbranchresiduestr);
    counter                = counter+1;
end

end

% function residueToRootChar = getonlyresiduetoroot(nreResidue)
% iterresidue = nreResidue;
% residueToRootChar='';
% count = 0;
% isParent = ~isempty(iterresidue.linkageParent);
% while(isParent)
%     residueTypeChar =  [iterresidue.stereoConfig.symbol,'-'];
%     resType = iterresidue.residueType;
%     if(~strcmpi(resType.name,'freeEnd'))
%         residueTypeChar = [residueTypeChar,resType.name];
%         residueTypeChar = [residueTypeChar,','];
%         residueTypeChar = [residueTypeChar,iterresidue.ringType.rSize];
%     else
%         residueTypeChar=resType.name;
%     end
%     
%     if(~isempty( iterresidue.linkageParent))
%         LinkageChar =iterresidue.anomer.symbol;
%     else
%         LinkageChar='';
%     end
%     residueLinkageChar= [LinkageChar,residueTypeChar];
%     residueToRootChar = [residueLinkageChar,residueToRootChar];
%     isParent = ~isempty(iterresidue.linkageParent);
%     if(isParent)
%         iterresidue = iterresidue.linkageParent.parent;
%         residueToRootChar = ['--',residueToRootChar];
%     end;
%     count =count +1;
% end
% end

