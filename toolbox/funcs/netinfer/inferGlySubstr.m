function varargout=inferGlySubstr(prodObj,enzObj)
% inferGlySubstr infer the substrate based on the enzyme and product.
%
%   [substrSpecies] = inferGlySubstr(prodObj,enz) infers
%    the substrate if the enzyme acts on the substrate to form the product
%    . If no substrate is  found, substrSpecies is empty.  
%    If more  than one substrates are found, they are stored in
%    a CellArrayList object substrSpecies.
%
%   [numSubstr,substrSpecies] = inferGlySubstr(prodObj,enz) returns the
%    number of substrates inferred. 
% 
%   [numSubstr,substrSpecies,rxns] = inferGlySubstr(prodObj,enz)
%   returns a list of reactions (a CellArrayList object rxns) if the enzyme acts
%   to form the product. If no substrate is found, rxns returns as an empty 
%   CellArrayList object.
%
%  [numSubstr,substrSpecies,rxns,pathway] = inferGlySubstr(prodObj,enz)
%   returns the pathway  if the enzyme acts to form the product prodObj. If
%   no substrate  is found, pathway is set empty.
%
%      Example 1:
%            mani  =GHEnz.loadmat('mani.mat');
%            m8species  = GlycanSpecies(glycanMLread('M8.glycoct_xml'));
%            [nsubstr, m9species] = inferGlySubstr(m8species,mani);
%            options  = displayset('showmass',true,'showLinkage',true,...
%                             'showRedEnd',true);
%            for i = 1: nsubstr
%                glycanViewer(m9species.get(i).glycanStruct,options);
%            end
%   
%      Example 2:
%             mani  =GHEnz.loadmat('mani.mat');
%            m8species  = GlycanSpecies(glycanMLread('M8.glycoct_xml'));
%            [nsubstr, m9species,m8rxns] = inferGlySubstr(m8species,mani);
%            for i = 1: nsubstr
%               glycanRxnViewer(m8rxns.get(i));
%            end 
%     
%      Example 3:
%            mani  =GHEnz.loadmat('mani.mat');
%            m8species  = GlycanSpecies(glycanMLread('M8.glycoct_xml'));
%            [nsubstr, m9species,m8rxns,m8pathway] = inferGlySubstr(m8species,mani);
%            glycanPathViewer(m8pathway);          
% 
% See also inferGlyProd,inferGlyRevrPath.

% Author: Gang Liu
% Date Last Updated: 8/2/13

if(isa(prodObj,'GlycanSpecies'))
    glycanObj =  prodObj.glycanStruct;
else
    error('MATLAB:GNAT:WrongInput','Input Wrong Type');
end

if(isa(enzObj,'GTEnz'));
    enzFuncToAdd=1;
    nreResidues =  glycanObj.getNonRedEndResidue;
    numTermRes  =  length(nreResidues);
elseif(isa(enzObj,'GHEnz'))
    enzFuncToAdd=0;
    allResidues           =  glycanObj.getAllResidues;
    numAllRes             =  length(allResidues);
else
    error('MATLAB:GNAT:NOTSUPPORTEDENZYME','NOT SUPPORTED ENZYME IN GNAT');
end

numSubstr=0;
substrSpecies=CellArrayList;
if(nargout==3)
    rxns=CellArrayList;
elseif(nargout==4)
    rxns=CellArrayList;
    path=Pathway;
end

isSubstrValid = true;

% rule 1: N-, O-, glycolipids
if(~isempty(enzObj.glycanTypeSpec))
    requirementglycantype     = enzObj.glycanTypeSpec;
    acceptortype    = glycanObj.glycanTypeSpec;
    isSubstrValid   = isequal(acceptortype,requirementglycantype);
end

isPathUpdated = false;
if(isSubstrValid)    
    if(enzFuncToAdd)      % remove the residue from terminal
        funcresrequirement    = enzObj.dispFuncResLink;
        for i = 1 : numTermRes
            nreResidue          = nreResidues{1,i};
            
            % rule 3: terminal residue must be function group defined in the
            % enzyme obj            
            terminalInfo   = nreResidue.dispResidueInfo;
            isValidRxnPos  = isequal(terminalInfo,funcresrequirement);
            
            if(~isValidRxnPos)
                continue;
            end
            
            % rule 4: if the the substrate should target the terminal residue
            if(~isempty(enzObj.isTerminalTarget))
                if(enzObj.isTerminalTarget)
                    if(length(nreResidue.getParent.linkageChildren)==1);
                        isValidRxnPos = true;
                    else
                        isValidRxnPos = false;
                    end
                else
                    isValidRxnPos = true;   
                end
            else
                isValidRxnPos = true;
            end
            
            if(~isValidRxnPos)
                continue;
            end
            
            % rule 5:resiude next to terminal residue is the same as
            % acceptorResidue
            requirement        = enzObj.dispAttachResLink;
            nexterminalResidue = nreResidue.getParent;
            if(~isempty(nexterminalResidue))
                acceptorTerminal     = nexterminalResidue.dispResidueInfo;
            else
                acceptorTerminal    = [];
            end
            
            isexactmatch = isempty(strfind(requirement,'?'));
            if(isexactmatch)
                isValidRxnPos = isequal(acceptorTerminal,requirement);
            elseif(numel(acceptorTerminal)==numel(requirement));
                quotaindex =  strfind(requirement,'?');
                requirement(quotaindex)='';
                acceptorTerminal(quotaindex)='';
                isValidRxnPos = isequal(acceptorTerminal,requirement);
            else
                isValidRxnPos = false;
            end
            
            if(~isValidRxnPos)
                continue;
            end
            
            % rule 6: the product should not be greater than the specified structure
            if(~isempty(enzObj.prodMinStruct))
                isValidRxnPos = isValidRxnPos && ...
                    glycanObj.contains(enzObj.prodMinStruct);
                if(~isValidRxnPos)
                   continue;
                end
            end
            
            % rule 7: the enzyme should not work on the children of specificed structure
            if(~isempty(enzObj.prodMaxStruct))
                isValidRxnPos = isValidRxnPos && ...
                    enzObj.prodMaxStruct.contains(glycanObj);
                if(~isValidRxnPos)
                  continue;
                end
            end
                 
            
            % rule 8: the enzyme should transfer the residue to the target branch
            if(~isempty(enzObj.targetBranch))  % if the rule exists
                chartargetbranch   = enzObj.targetBranch.name;
                charresidue        = glycanObj.getresiduetoroot(nreResidue.getParent);
                isValidRxnPos      = isValidRxnPos && ...
                    (strcmpi(chartargetbranch,charresidue));                
                if(~isValidRxnPos)
                  continue;
                end
            end
            
           
            
            % rule 9: the enzyme should transfer the residue to the target branch containing
            % specified structure
            if(~isempty(enzObj.targetbranchcontain)) % if the rule exists
                chartargetbranch = enzObj.targetbranchcontain.name;
                chartargetbranch = regexprep(chartargetbranch,'freeEnd--',' ');
                chartargetbranch = regexprep(chartargetbranch,'?',' ');
                
                charresidue= glycanObj.getresiduetoroot(nreResidue.getParent);
                charresidue= regexprep(charresidue,'freeEnd--','');
                charresidue= regexprep(charresidue,'?','');
                
                isValidRxnPos = isValidRxnPos && (~isempty(strfind(strtrim(charresidue),...
                    strtrim(chartargetbranch))));
                
                if(~isValidRxnPos)
                  continue;
                end
            end
            
            % rule 10: the enzyme should not act on the specificed branch
            if(~isempty(enzObj.targetNABranch))  % if the rule exists
                if(isa(enzObj.targetNABranch,'GlycanStruct'))
                    chartargetbranch = enzObj.targetNABranch.name;
                    charresidue           = glycanObj.getresiduetoroot(...
                        nreResidue.getParent);
                    isValidRxnPos         = isValidRxnPos && ...
                        (~strcmp(chartargetbranch,charresidue));
                elseif(isa(enzObj.targetNABranch,'CellArrayList'))
                    for ii = 1 : length(enzObj.targetNABranch)
                        chartargetbranch  = enzObj.targetNABranch.name;
                        charresidue       = glycanObj.getresiduetoroot(...
                            nreResidue.getParent);
                        isValidRxnPos     = isValidRxnPos && ...
                            (~strcmp(chartargetbranch,charresidue));
                    end
                end
                
               if(~isValidRxnPos)
                  continue;
               end
            end
            
            glycanSubstrObj     = glycanObj.clone;
            nreResidues         = glycanSubstrObj.getNonRedEndResidue;
            nreResiduetoremove  = nreResidues{1,i};
            
            isResidueDeleted    = glycanSubstrObj.removeNonRedEndResidue(...
                nreResiduetoremove);
            
            %rule 2: It does not contain the specified residue
            if(~isempty(enzObj.substNAResidue))
                NAResiduename = enzObj.substNAResidue.name;
                composition   = glycanSubstrObj.getComposition;
                isValidRxnPos  =  ~(isfield(composition,NAResiduename) ...
                    && (composition.(NAResiduename)>=1));
                
               if(~isValidRxnPos)
                  continue;
               end
            end
            
           
           % rule 11 the substrate should not contain the specificed branch
            if(~isempty(enzObj.substNABranch))  % if the rule exists
                if(isa(enzObj.substNABranch,'GlycanStruct'))
                    isValidRxnPos = isValidRxnPos &&...
                        (~glycanSubstrObj.contains(enzObj.substNABranch));
                elseif(isa(enzObj.substNABranch,'CellArrayList'))
                    for ii = 1 : length(enzObj.substNABranch)
                        isValidRxnPos = isValidRxnPos && ...
                (~glycanSubstrObj.contains(enzObj.substNABranch.get(ii)));
                    end
                end
                
                if(~isValidRxnPos)
                  continue;
                end
            end
            
            % rule 12, the enzyme should not work on the parent of specified structure
            if(~isempty(enzObj.substMaxStruct)) % if the rule exists
                isValidRxnPos = isValidRxnPos && ...
                    enzObj.substMaxStruct.contains(glycanSubstrObj);
                
               if(~isValidRxnPos)
                 continue;
               end
            end
            
            % rule 13, the enzyme should not work on the children of specified structure
            if(~isempty(enzObj.substMinStruct))  % if the rule exists
                isValidRxnPos = isValidRxnPos &&...
                    glycanSubstrObj.contains(enzObj.substMinStruct);
               if(~isValidRxnPos)
                continue;
               end
            end
            
            if(isResidueDeleted)
                %  disp('has substrate');
                numSubstr=numSubstr+1;
                substr = GlycanSpecies(glycanSubstrObj);
                substrSpecies.add(substr);
                if(nargout>=3)
                  rxns.add(Rxn(substr,prodObj,enzObj));
                  if(nargout==4)
                     isPathUpdated = isPathUpdated || isResidueDeleted;
                  end
                end
            end
        end
        
    else  % glycosidase
        requireattach2funcgroup = enzObj.dispAttachResLink; % attach to functional group
        isexactmatch            = isempty(...
               strfind(requireattach2funcgroup,'?'));

        for  i = 1  : numAllRes
            %i
            nreResidue    =  allResidues{1,i};            
            isValidRxnPos =  true;
            
            % disp(nreResidue.residueType.name);
            if(strcmpi(nreResidue.residueType.name,'freeEnd'))
                continue;
            end
            
            % check requirement 0
            if(~isempty(enzObj.resAtt2FG))  % if the rule exists
                residueTypeNextTerm       = nreResidue.dispResidueInfo;                
                if(isexactmatch)
                    isValidRxnPos         = isequal(residueTypeNextTerm,...
                        requireattach2funcgroup);
                elseif(numel(residueTypeNextTerm)==numel(requireattach2funcgroup));
                    quotaindex =  strfind(requireattach2funcgroup,'?');
                    tmprequireattach2funcgroup = requireattach2funcgroup;
                    tmprequireattach2funcgroup(quotaindex)='';
                    residueTypeNextTerm(quotaindex)='';
                    isValidRxnPos = isequal(residueTypeNextTerm,tmprequireattach2funcgroup);
                else
                    isValidRxnPos = false;
                end
            else
               error('MATLAB:GNAT:INCOMPLETEENZDEFINITION','INCOMPLETE ENZYME DEFINITION');
            end
            
            if(~isValidRxnPos)
                continue;
            end
            
            % rule 6: the product should not be greater than the specified structure
            if(~isempty(enzObj.prodMinStruct))
                isValidRxnPos = glycanObj.contains(enzObj.prodMinStruct);
                if(~isValidRxnPos)
                  continue;
                end
            end      
                        
            % rule 7: the product should not work on the children of specificed structure
            if(~isempty(enzObj.prodMaxStruct))
                isValidRxnPos = enzObj.prodMaxStruct.contains(glycanObj);
                if(~isValidRxnPos)
                  continue;
                end
            end            
          
            % rule 9: the enzyme should not transfer the residue to the target branch containing
            % specified structure
            if(~isempty(enzObj.targetbranchcontain)) % if the rule exists
                chartargetbranch=enzObj.targetbranchcontain.name;
                chartargetbranch = regexprep(chartargetbranch,'freeEnd','');
                
                charresidue= glycanObj.getresiduetoroot(nreResidue.getParent);
                charresidue= regexprep(charresidue,'freeEnd','');
                
                isValidRxnPos = (~isempty(strfind(chartargetbranch,charresidue)));
                if(~isValidRxnPos)
                   continue;
                end
            end
            
            % rule 10: the enzyme should not act on the specificed branch
            if(~isempty(enzObj.targetNABranch))  % if the rule exists
                if(isa(enzObj.targetNABranch,'GlycanStruct'))
                    chartargetbranch    = enzObj.targetNABranch.name;
                    charresidue         = glycanObj.getresiduetoroot(...
                        nreResidue.getParent);
                    isValidRxnPos       = isValidRxnPos && (~strcmp(...
                        chartargetbranch,charresidue));
                elseif(isa(enzObj.targetNABranch,'CellArrayList'))
                    for ii = 1 : length(enzObj.targetNABranch)
                        chartargetbranch = enzObj.targetNABranch.name;
                        charresidue      = glycanObj.getresiduetoroot(...
                            nreResidue.getParent);
                        isValidRxnPos    =  isValidRxnPos && ...
                            (~strcmp(chartargetbranch,charresidue));
                    end
                end
            end
            
            if(~isValidRxnPos)
                continue;
            end
            
            nBondChoice = length(enzObj.linkFG.bond);
            glycanSubstStructArray = CellArrayList;
            isOneResiduedeleted = false;
            
            if(nBondChoice>=1)
                for j = 1 : nBondChoice
                    bond=enzObj.linkFG.bond(j,1);
                    glycanSubstrObj      = glycanObj.clone;
                    nreResidues          = glycanSubstrObj.getAllResidues;
                    nreResidue           = nreResidues{1,i};
                    isResidueDeleted     = glycanSubstrObj.addResidue(nreResidue,...
                        enzObj.resfuncgroup,enzObj.linkFG.anomer,...
                        bond);
                    if(isResidueDeleted)
                        isOneResiduedeleted = true;                        
                        glycanSubstStructArray.add(glycanSubstrObj);
                    end
                end
            end
            
            if(~isOneResiduedeleted)
                %                 disp('debug 1');
                continue;
            end
            
            for k = 1: glycanSubstStructArray.length                
                glycanSubstrObj = glycanSubstStructArray.get(k);
                isValidRxnPos   = true;
                % rule 11 the substrate should not contain the specificed branch
                if(~isempty(enzObj.substNABranch))  % if the rule exists
                    if(isa(enzObj.substNABranch,'GlycanStruct'))
                        isValidRxnPos = isValidRxnPos && (~glycanSubstrObj.contains(enzObj.substNABranch));
                    elseif(isa(enzObj.substNABranch,'CellArrayList'))
                        for ii = 1 : length(enzObj.substNABranch)
                            isValidRxnPos = isValidRxnPos && (~glycanSubstrObj.contains(enzObj.substNABranch.get(ii)));
                        end
                    end
                end
                
                if(~isValidRxnPos)
                    continue;
                end
                
                % rule 12, the enzyme should not work on the parent of specified structure
                if(~isempty(enzObj.substMaxStruct))  % if the rule exists
                    isValidRxnPos = isValidRxnPos &&...
                        enzObj.substMaxStruct.contains(glycanSubstrObj);
                end
                
                if(~isValidRxnPos)
                    continue;
                end
                
                % rule 13, the enzyme should not work on the children of specified structure
                if(~isempty(enzObj.substMinStruct))  % if the rule exists
                    isValidRxnPos = isValidRxnPos && ...
                        glycanSubstrObj.contains(enzObj.substMinStruct);
                end
                
                if(~isValidRxnPos)
                    continue;
                end
                
                numSubstr=numSubstr+1;
                substr =GlycanSpecies(glycanSubstrObj);
                substrSpecies.add(substr);
                if(nargout>=3)
                   rxns.add(Rxn(substr,prodObj,enzObj));
                  if(nargout==4)
                     isPathUpdated = true;
                  end
                end
            end
        end
    end

    if(nargout==4)
      if(isPathUpdated)
        path.addGlycans(substrSpecies);
        path.addGlycan(prodObj);
        path.addRxns(rxns);
      end
    end 
end

if(nargout==1)
    varargout{1}=substrSpecies;
elseif(nargout==2)
    varargout{1}=numSubstr;
    varargout{2}=substrSpecies;
elseif(nargout==3)
    varargout{1}=numSubstr;
    varargout{2}=substrSpecies;
    varargout{3}=rxns;
elseif(nargout==4)
    varargout{1}=numSubstr;
    varargout{2}=substrSpecies;
    varargout{3}=rxns;
    varargout{4}=path;
end

end

