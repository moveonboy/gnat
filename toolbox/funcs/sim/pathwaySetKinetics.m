function glyPathway = pathwaySetKinetics(glyPathway,enzkineticsdb,option,varargin)
%pathwaySetKinetics set enzymatic kinetics for each reaction in the pathway
%
%  Input: glyPathway, kineticsdb, isSubstCompt
%      1. glyPathway, Pathway Object
%
%      2. enzkineticsdb, a database storing enzyme kinetics for each
%         comparmnet.
%
%      3. option, a structure with three fields substrCompt,a logical value
%          indicating if substrate competition is considered.
%
%  Output: glyPathway
%       glyPathway,  Pathway object returned
%

%See also GlycanNetModel.

% Author: Gang Liu
% Date Lastly Updated: 4/8/14

if(length(varargin)==1)
    sugarnucldb = varargin{1};
    enzDistr    = Containers.Map; 
elseif(length(varargin)==2)
    sugarnucldb = varargin{1};
    enzDistr    = varargin{2};
end

if(~isa(sugarnucldb,'containers.Map'))
    error('MATLAB:GNAT:ERRORINPUTTYPE','ERROR INPUT TYPE');
end

if(~isa(enzDistr,'containers.Map'))
    error('MATLAB:GNAT:ERRORINPUTTYPE','ERROR INPUT TYPE');
end

% reset listofrxns
glyPathway.setListOfRxns;
theRxns = glyPathway.listofRxns;
numRxns = length(theRxns);
theEnzs = glyPathway.theEnzs;
numEnzs = length(theEnzs);

% set enzyme kinetics for each enzyme 
for i =1: numEnzs
    theithEnz = theEnzs.get(i);
    enzForm   = theithEnz.name;  
    if(enzkineticsdb.isKey(enzForm))
        theithEnz.enzkinetics     = enzkineticsdb(enzForm);
        theithEnz.enzkinetics.enz = theithEnz;
    else
        enzForm
        error('MATLAB:GNAT:ERRORDATABASEINPUT',...
            'ENZYMET KINETICS DATABASE DOES NOT CONTAIN ENZYME KINETICS FOR THE ENZYME')
    end
    
    if(isempty(theithEnz.enzkinetics))
         error('MATLAB:GNAT:ERRORDATABASEINPUT',...
            'ENZYMET KINETICS DATABASE DOES NOT CONTAIN ENZYME KINETICS FOR THE ENZYME')
    end
end

%set up Vm, Km 
for i = 1:numRxns
    % theithrxn = glyPathway.theRxns.get(i);
    theithrxn = theRxns{i};
    if(~isempty(theithrxn.enz))
        theithrxn.setrxnkinetics;
    end
end


% if enzyme distribution of multiple compartment is provided
numEnzDistrs = enzDistr.Count;
if(numEnzs~=0)
    enznames = enzDistr.keys;
    for j = 1 : numEnzDistrs;
        jthenzname = enznames{j};
        ithEnzDist = enzDistr(jthenzname);
        % theRxns = glyPathway.theRxns;
        for i = 1 : numRxns
          % ithRxn = theRxns.get(i); 
           ithRxn = theRxns{i}; 
           if(strcmpi(ithRxn.enz.name,jthenzname))
               comptname = ithRxn.reac.compartment.name;
               if(isfield(ithEnzDist,comptname))
                  if(isa(ithRxn.rxnkinetics,'BiBiRxnKinetics'))
                    ithRxn.rxnkinetics.bibmech5_substratekineticsratio.Vm_ratio =...
                      ithEnzDist.(comptname);
                  elseif(isa(ithRxn.rxnkinetics,'MMenRxnKinetics'))
                     ithRxn.rxnkinetics.substratekineticsratio.Vm_ratio =...
                      ithEnzDist.(comptname); 
                  end
               else
                   error('MATLAB:GNAT:ERRORINPUT','INCORRECT ENZYME DISTRIBUTION');
               end
           end
        end           
    end
end    

% set substrate inhibition if necessary
comptRxnListdb = containers.Map;
for i = 1 : numEnzs
    ithEnz = theEnzs.get(i);
    comptRxnListBySameEnz     = glyPathway.findRxnsBySameEnz(ithEnz);
    comptRxnListdb(ithEnz.id) = comptRxnListBySameEnz;
end

if(option.substrCompt)
    for i = 1 : numRxns
        theithrxn = theRxns{i};
        if(~isempty(theithrxn.enz))
%             try
                theCompetingRxnsList = comptRxnListdb(theithrxn.enz.id);
                theCompetingRxnsList = theCompetingRxnsList.*(theCompetingRxnsList~=i);
                theCompetingRxnsList = theCompetingRxnsList(find(theCompetingRxnsList));
                if(~isempty(theCompetingRxnsList))
                    theCompetingRxns     = theRxns(theCompetingRxnsList);
                    theithrxn.rxnkinetics.setRxnKineticsInhibitor(theCompetingRxns);
                end
%             catch err
%                 disp('error in competing rxn');
%             end
        end
    end
end

%set up kinetic rate formula for each enzymatic reaction
for i = 1 : numRxns
    theithrxn = theRxns{i};
    if(~isempty(theithrxn.enz))
        if(isa(theithrxn.enz.enzkinetics,'MMenKinetics'))
            theithrxn.rxnkinetics.setMathFormula(theithrxn.enz.id,...
                theithrxn.id);
        elseif(isa(theithrxn.enz.enzkinetics,'BiBiKinetics'))  % if the enzyme is BiBiKinetics
            if(strcmpi(theithrxn.enz.donor,'UDP-N-acetyl-D-glucosamine'))
                sugarnucl = 'UDP_GlcNAc';
                nucl = 'UDP';
            elseif(strcmpi(theithrxn.enz.donor,'UDP-alpha-D-galactose'))
                sugarnucl = 'UDP_Gal';
                nucl = 'UDP';
            elseif(strcmpi(theithrxn.enz.donor,'CMP-N-acetylneuraminate'))% CMP-Sialic Acid
                sugarnucl = 'CMP_NeuAc';
                nucl = 'CMP';
            elseif(strcmpi(theithrxn.enz.donor,'GDP-beta-L-fucose'))% CMP-Sialic Acid
                sugarnucl = 'GDP_Fuc';
                nucl = 'GDP';
            else
                error('MATLAB:GNAT:NOTSUPPORTEDSUGAR','SUGAR NAME IS NOT SUPPORTED');
            end
            
            sugarnames  = {sugarnucl,nucl};
            try 
                sugarvalues = [sugarnucldb(sugarnucl),sugarnucldb(nucl)];
            catch err
                error('MATLAB:GNAT:NOTSUPPORTEDSUGAR','SUGAR NAME IS NOT SUPPORTED IN THE DATABASE');
            end
                            
            theithrxn.rxnkinetics.setMathFormula(theithrxn.enz.id,...
                theithrxn.id,sugarnames,sugarvalues,0);
        end
    end
end

end

