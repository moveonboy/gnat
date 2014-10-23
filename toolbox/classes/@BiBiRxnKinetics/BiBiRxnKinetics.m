classdef BiBiRxnKinetics < RxnKinetics
        
    properties
        enzkinetics;
        reacspeciesid1; 
        reacspeciesid2;
        prodspeciesid1;
        prodspeciesid2;
        bibmech5_substratekinetics = struct('Vm',[],'Kmd_G1E',[],'Kmd_SNE',[],'Keq',[])
        bibmech5_substratekineticsratio = struct('Vm_ratio',1,'Kmd_G1E_ratio',1,'Kmd_SNE_ratio',1,'Keq_ratio',1)
        substinhibkinetics; %=struct('Ki',[],'speciesid','');
        kineticlaws
    end
    
    methods % Constructor for BiBiRxnKinetics define substrate specific rate constants
        function obj = BiBiRxnKinetics(varargin)
            %enzkinetics,reac,prod)
         if(isempty(varargin))   
             return
         end
         
         obj.enzkinetics = varargin{1};
         reac = varargin{2};
         prod = varargin{3};
         
         if(isa(obj.enzkinetics.enz,'GTEnz'))
            obj.reacspeciesid1  = reac.id;  % glycan species
            obj.prodspeciesid1  = prod.id;  % glycan species
            obj.reacspeciesid2  = obj.enzkinetics.enz.donor; % donor 
            obj.prodspeciesid2  = obj.enzkinetics.enz.donorProd; % donar product
            obj.reacspeciesid2  = regexprep(obj.reacspeciesid2,...
                '-','_');
            obj.prodspeciesid2  = regexprep(obj.prodspeciesid2,...
                '-','_');            
            
          % obj.substinhibkinetics                = CellArrayList;
            
          % derive substrate-specific rate constant.
          % Apply specificity rule if necessary
            
          if(strcmpi(obj.enzkinetics.mechanism,'Sequential Order with substrate Inhibition'))
            if isempty(obj.enzkinetics.substrspecdb)...
                    && isempty(obj.enzkinetics.substrspecarray)
                obj.bibmech5_substratekinetics.Vm      = obj.enzkinetics.sumconsts.Vm;
                obj.bibmech5_substratekinetics.Kmd_G1E = obj.enzkinetics.sumconsts.Kmd_G1E;
                obj.bibmech5_substratekinetics.Kmd_SNE = obj.enzkinetics.sumconsts.Kmd_SNE;
                obj.bibmech5_substratekinetics.Keq     = obj.enzkinetics.sumconsts.Keq;
            else
                obj.bibmech5_substratekinetics.Vm      = obj.enzkinetics.sumconsts.Vm;
                obj.bibmech5_substratekinetics.Kmd_G1E = obj.enzkinetics.sumconsts.Kmd_G1E;
                obj.bibmech5_substratekinetics.Kmd_SNE = obj.enzkinetics.sumconsts.Kmd_SNE;
                obj.bibmech5_substratekinetics.Keq     = obj.enzkinetics.sumconsts.Keq;
                
                if(obj.enzkinetics.specificityrule==1)
                    glycanstring = reac.glycanStruct.toString;
                    if(isKey(obj.enzkinetics.substrspecdb,glycanstring))
                        ratio = obj.enzkinetics.substrspecdb(glycanstring);
                        obj.bibmech5_substratekineticsratio.Vm_ratio      = ratio.Vm ;
                        obj.bibmech5_substratekineticsratio.Kmd_G1E_ratio = ratio.Kmd_G1E ;
                        obj.bibmech5_substratekineticsratio.Kmd_SNE_ratio = ratio.Kmd_SNE ;
                        obj.bibmech5_substratekineticsratio.Keq_ratio     = ratio.Keq ;
                    end
                elseif(obj.enzkinetics.specificityrule==2)
                    for i = 1: length(obj.enz.substrspecarray)
                        substruct = obj.enzkinetics.substrspecarray.get(i).struct;
                        ratio     = obj.enzkinetics.substrspecarray.get(i).ratio;
                        if(glycanspecies.glycanStruct.contains(substruct))
                            obj.bibmech5_substratekineticsratio.Vm_ratio      = ratio.Vm ;
                            obj.bibmech5_substratekineticsratio.Kmd_G1E_ratio = ratio.Kmd_G1E ;
                            obj.bibmech5_substratekineticsratio.Kmd_SNE_ratio = ratio.Kmd_SNE ;
                            obj.bibmech5_substratekineticsratio.Keq_ratio     = ratio.Keq ;
                        end
                    end
                end                
            end
          else
              error('MATLAB:GNAT:NOTSUPPORTEDMECHANISM','MECHANISM NOT SUPPORTED YET');
          end
        elseif(isa(obj.enzkinetics.enz,'GHEnz'))
             % to be written
             
        end
       end
        
       function setRxnKineticsInhibitor(obj,theCompetingRxns) % does not check the rxns 
            usedspecid  = cell(length(theCompetingRxns),1);
            count       = 0;
            if(isa(theCompetingRxns,'Rxn'))
              newsubstrid = theCompetingRxns.reac.id;
              substinhibkinetic = struct('Ki_md_G1E',...
              theCompetingRxns.rxnkinetics.bibmech5_substratekinetics.Kmd_G1E,...
                 'Ki_md_G1E_ratio', theCompetingRxns.rxnkinetics.bibmech5_substratekineticsratio.Kmd_G1E_ratio,...
                 'speciesid',newsubstrid);
              obj.substinhibkinetics = substinhibkinetic;            
            elseif(iscell(theCompetingRxns));
                numCompRxns = length(theCompetingRxns);
                for i = 1 : numCompRxns
                    newsubstrid       = theCompetingRxns{i}.reac.id;
                    ithrxnkinetics    = theCompetingRxns{i}.rxnkinetics;
                    substinhibkinetic = struct('Ki_md_G1E',...
                       ithrxnkinetics.bibmech5_substratekinetics.Kmd_G1E,...
                       'Ki_md_G1E_ratio', ithrxnkinetics.bibmech5_substratekineticsratio.Kmd_G1E_ratio,...
                        'speciesid',newsubstrid);
                    if(sum(strcmpi(newsubstrid,usedspecid)~=0))
                        continue;
                    end
                    count=count+1;
                    usedspecid{count} =  newsubstrid;
                    listofsubstinhibkinetics(count,1)=substinhibkinetic;
                end
                obj.substinhibkinetics = listofsubstinhibkinetics;
            end
       end
    end
        
    %method for mathematical formula
    methods
        function setMathFormula(obj, enzid, rxnid, constspecids, constspecievalues,varargin) 
            % structure: 1)names; 2) values
            if (length(varargin)==1)
                ifnegrev = varargin{1};
            else
                ifnegrev = false;
            end
            
            if(strcmpi(obj.enzkinetics.mechanism,'Sequential Order with substrate Inhibition'))
                if(isempty(obj.substinhibkinetics))                                
                    obj.kineticlaws.parametervalues = [obj.bibmech5_substratekinetics.Vm,...
                        obj.bibmech5_substratekinetics.Kmd_G1E,obj.bibmech5_substratekinetics.Kmd_SNE,...
                        obj.bibmech5_substratekinetics.Keq];
                    obj.kineticlaws.parameternames = {['Vm',enzid],...
                        ['KmG',enzid],['KmS',enzid],...
                        ['Kq',enzid]};
                    
                    obj.kineticlaws.parameternames{end+1}=constspecids{1}; % first one is sugar donor
                    obj.kineticlaws.parameternames{end+1}=constspecids{2}; % second one is nucleotide
                    
                    obj.kineticlaws.parametervalues(end+1)=constspecievalues(1);
                    obj.kineticlaws.parametervalues(end+1)=constspecievalues(2);
                    
                    if(ifnegrev)
                        obj.kineticlaws.mathformula = ...
                        'Vm*Vm_ratio*(SUGARNUCL*Pi-(Pplus1*NUCL)/Keq/Keq_ratio)/(Kmd_G1E*Kmd_G1E_ratio*(Kmd_SNE*Kmd_SNE_ratio+SUGARNUCL)*(1+Pi/Kmd_G1E/Kmd_G1E_ratio))';
                    else
                        obj.kineticlaws.mathformula = ...
                       'Vm*Vm_ratio*SUGARNUCL*Pi/(Kmd_G1E*Kmd_G1E_ratio*(Kmd_SNE*Kmd_SNE_ratio+SUGARNUCL)*(1+Pi/Kmd_G1E/Kmd_G1E_ratio))';
                    end
                   
                    obj.kineticlaws.speciesnames = {obj.reacspeciesid1,obj.prodspeciesid1};% a cell array                    
                    
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Pi',obj.reacspeciesid1);  % glycan reactant 
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Pplus1',obj.prodspeciesid1); % glycan product
                    
                    if(obj.bibmech5_substratekineticsratio.Vm_ratio~=0)
                        obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                            'Vm_ratio',num2str(obj.bibmech5_substratekineticsratio.Vm_ratio));
                    else
                        obj.kineticlaws.mathformula  ='0';
                        return;
                    end
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_G1E_ratio',num2str(obj.bibmech5_substratekineticsratio.Kmd_G1E_ratio));
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_SNE_ratio',num2str(obj.bibmech5_substratekineticsratio.Kmd_SNE_ratio));
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Keq_ratio',num2str(obj.bibmech5_substratekineticsratio.Keq_ratio)); 
                   
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Vm',obj.kineticlaws.parameternames{1});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_G1E',obj.kineticlaws.parameternames{2});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_SNE',obj.kineticlaws.parameternames{3});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Keq',obj.kineticlaws.parameternames{4});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'SUGARNUCL',sprintf('%.3f',obj.kineticlaws.parametervalues(5)));
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'NUCL',sprintf('%.3f',obj.kineticlaws.parametervalues(6)));
                    
                else
                  %   numsubsinhib = length(obj.substinhibkinetics);
                    obj.kineticlaws.parametervalues = [obj.bibmech5_substratekinetics.Vm,...
                        obj.bibmech5_substratekinetics.Kmd_G1E,obj.bibmech5_substratekinetics.Kmd_SNE,...
                        obj.bibmech5_substratekinetics.Keq];
                    obj.kineticlaws.parameternames = {['Vm',enzid],...
                        ['KmG',enzid],['KmS',enzid],...
                        ['Ke',enzid]};
                    
                    obj.kineticlaws.parameternames{end+1}=constspecids{1}; % first one is sugar donor
                    obj.kineticlaws.parameternames{end+1}=constspecids{2}; % second one is nucleotide
                    
                    obj.kineticlaws.parametervalues(end+1)=constspecievalues(1);
                    obj.kineticlaws.parametervalues(end+1)=constspecievalues(2);
                    
                    obj.kineticlaws.speciesnames = {obj.reacspeciesid1,obj.prodspeciesid1};% a cell array                    
                    
                    if(ifnegrev)
                       obj.kineticlaws.mathformula = ...
                          'Vm*Vm_ratio*(SUGARNUCL*Pi-(Pplus1*NUCL)/Keq/Keq_ratio)/(Kmd_G1E*Kmd_G1E_ratio*(Kmd_SNE*Kmd_SNE_ratio+SUGARNUCL)*(1+Pi/Kmd_G1E/Kmd_G1E_ratio';
                    else
                      obj.kineticlaws.mathformula = ...
                         'Vm*Vm_ratio*SUGARNUCL*Pi/(Kmd_G1E*Kmd_G1E_ratio*(Kmd_SNE*Kmd_SNE_ratio+SUGARNUCL)*(1+Pi/Kmd_G1E/Kmd_G1E_ratio';
                    end
                    
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Pi',obj.reacspeciesid1);
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Pplus1',obj.prodspeciesid1);
                    
                    
                    if(obj.bibmech5_substratekineticsratio.Vm_ratio==0)
                        obj.kineticlaws.mathformula = '0';
                        return
                    else 
                        obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Vm_ratio',num2str(obj.bibmech5_substratekineticsratio.Vm_ratio));
                    end
                    
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_G1E_ratio',num2str(obj.bibmech5_substratekineticsratio.Kmd_G1E_ratio));
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_SNE_ratio',num2str(obj.bibmech5_substratekineticsratio.Kmd_SNE_ratio));
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Keq_ratio',num2str(obj.bibmech5_substratekineticsratio.Keq_ratio)); 
                    
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Vm',obj.kineticlaws.parameternames{1});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_G1E',obj.kineticlaws.parameternames{2});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Kmd_SNE',obj.kineticlaws.parameternames{3});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'Keq',obj.kineticlaws.parameternames{4});
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'SUGARNUCL',sprintf('%.3f',obj.kineticlaws.parametervalues(5)));
                    obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,...
                        'NUCL',sprintf('%.3f',obj.kineticlaws.parametervalues(6)));
                    for i = 1 : length(obj.substinhibkinetics)
                            ithsubstinhibkinetics = obj.substinhibkinetics(i,1);
                            obj.kineticlaws.parameternames{end+1}  = ['Ki_md_G1E_',enzid];
                            try
                                obj.kineticlaws.parametervalues(end+1) = ithsubstinhibkinetics.Ki_md_G1E;
                                obj.kineticlaws.parameternames{end+1}  = ['Ki_md_G1E_ratio',enzid,...
                                                   ithsubstinhibkinetics.speciesid];
                                obj.kineticlaws.parametervalues(end+1) = ithsubstinhibkinetics.Ki_md_G1E_ratio;

                                obj.kineticlaws.speciesnames{end+1}= ithsubstinhibkinetics.speciesid;

                                obj.kineticlaws.mathformula = [obj.kineticlaws.mathformula, '+',...
                                    ithsubstinhibkinetics.speciesid];
                                obj.kineticlaws.mathformula = [obj.kineticlaws.mathformula, ...
                                    '/KmG',enzid];                            
                                obj.kineticlaws.mathformula = [obj.kineticlaws.mathformula, '/',...
                                     num2str(ithsubstinhibkinetics.Ki_md_G1E_ratio)];  
                            catch
                                disp('debug in BiBi');
                            end
                            
                    end
                    obj.kineticlaws.mathformula=[obj.kineticlaws.mathformula,'))'];
%                 else
%                     error('MATLAB:GNAT:MECHANISMNOTSUPPORTED','NOT SUPPORTED MECHANISM');
                end
            else
                
            end
            
            obj.kineticlaws.mathformula = regexprep(obj.kineticlaws.mathformula,'*1(?=[^0-9])','');
            obj.kineticlaws.mathformula  = regexprep(obj.kineticlaws.mathformula ,'\1(?=[^0-9])','');
            obj.kineticlaws.mathformula  = regexprep(obj.kineticlaws.mathformula ,'/1(?=[^0-9])','');
        end
    end
    
    methods
        function new = clone(obj)
               
            % Instantiate new object of the same class.
            new = feval(class(obj));
            
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                if(isa(obj.(p{i}),'handle'))
                    if(isa(obj.(p{i}),'CellArrayList'))
                        new.(p{i})=CellArrayList;
                        for j = 1 : obj.(p{i}).length
                            new.(p{i}).add(obj.(p{i}).get(j))
                        end
                    elseif(isa(obj.(p{i}),'containers.Map'))
                        allkeys    = obj.(p{i}).keys;
                        allvalues  = obj.(p{i}).values;
                        new.(p{i}) = containers.Map(allkeys,allvalues);
                    else
                        new.(p{i}) = obj.(p{i}).clone;
                    end
                else
                    new.(p{i}) = obj.(p{i});
                end
            end
        end
    end    
end