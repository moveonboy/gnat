function glyPathway = pathwaySetEnzDistribution(glyPathway,option,varargin)
%pathwaySetKinetics set enzymatic kinetics for each reaction in the pathway
%
%  Input: glyPathway, isSubstCompt
%      1. glyPathway, Pathway Object
%
%      2. option, a structure with two fields substrCompt, enzIDorName,a logical value
%          indicating if substrate competition is considered.
%
%  Output: glyPathway
%       glyPathway,  Pathway object returnned
%
%See also GlycanNetModel.

% Author: Gang Liu
% Date Lastly Updated: 4/8/14

if(length(varargin)==1)
    sugarnucldb=varargin{1};
end

% set enzyme kinetics for each enzyme
for i =1: length(glyPathway.theEnzs)
    theithEnz = glyPathway.theEnzs.get(i);
    if(option.enzIDorName==1)
       enzForm = theithEnz.id;
    elseif(option.enzIDorName==2)
       enzForm = theithEnz.name;  
    else
       error('MATLAB:GNAT:NOTSUPPORTEDFORM',...
            'NOT CURRENTLY SUPPORTED FORM')     
    end
    if(enzkineticsdb.isKey(enzForm))
        theithEnz.enzkinetics = enzkineticsdb(enzForm);
        theithEnz.enzkinetics.enz = glyPathway.theEnzs.get(i);
    else
        enzForm
        error('MATLAB:GNAT:ERRORDATABASEINPUT',...
            'ENZYMET KINETICS DATABASE DOES NOT CONTAIN ENZYME KINETICS FOR THE ENZYME')
    end
end

%set up Vm, Km
for i = 1:glyPathway.theRxns.length
    theithrxn = glyPathway.theRxns.get(i);
    if(~isempty(theithrxn.enz))
        theithrxn.setrxnkinetics;
    end
end

% set substrate inhibition if necessary
if(option.substrCompt)
    for i = 1 : glyPathway.theRxns.length
        theithrxn = glyPathway.theRxns.get(i);
        if(~isempty(theithrxn.enz))
            theCompetingRxns = glyPathway.findEnzCompetingRxns(theithrxn);
            theithrxn.rxnkinetics.setRxnKineticsInhibitor(theCompetingRxns);
        end
    end
end

%set up kinetic rate formula for each enzymatic reaction
for i = 1 : glyPathway.theRxns.length
    theithrxn = glyPathway.theRxns.get(i);
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
                error('MATLAB:GNAT:NOTSUPPORTEDSUGAR','SUGAR NAME IS NOT SUPPORTED');
            end
                            
            theithrxn.rxnkinetics.setMathFormula(theithrxn.enz.id,...
                theithrxn.id,sugarnames,sugarvalues,0);
        end
    end
end

end

