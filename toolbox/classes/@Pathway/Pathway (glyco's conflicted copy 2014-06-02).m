classdef Pathway  < handle
    %Pathway class representing a glycosylation reaction network
    %
    % A Pathway object is a generic representation of a glycosylation
    %  reaction network. It consists of a number of species, reactions and
    %  enzymes involved.
    %
    % Pathway properties:
    %  theRxns             - a list of the reactions in the network
    %  theSpecies          - a list of the species in the network
    %  theEnzs             - a list of enzymes in the network
    %  initGlycans         - initial glycan structures
    %  finalGlycans        - final glycan structures
    %  compartment         - the compartment where the reactions take place
    %  name                - the name of the network
    %  SM                  - stoichiometry matrix
    %
    % Pathway methods:
    %  Pathway             - create a Pathway object
    %  setSM               - set stoichiometry matrix 
    %  getReactions        - retrieve the "theRxns" property
    %  getSpecies          - retrieve the "theSpecies" property
    %  getEnzymes          - retrieve the "theEnzs" property
    %  getInitGlycans      - retrieve the "initGlycans" property
    %  getFinalGlycans     - retrieve the "finalGlycans" property
    %  getCompartment      - retrieve the "compartment" property
    %  getGlycanStruct     - get the glycan structure of the ith species in the network
    %  getName             - retrieve the "name" property
    %  setReactions        - set the "theRxns" property
    %  setSpecies          - set the "theSpecies" property
    %  setEnzymes          - set the "theEnzs" property
    %  setInitGlycans      - set the "initGlycans" property
    %  setFinalGlycans     - set the "finalGlycans" property
    %  setCompartment      - set the "compartment" property
    %  setName             - set the "name" property
    %  isempty             - return true if no reactions is in the network
    %  clearSpeciesMark    - reset all the species' marks as false
    %  clearRxnsMark       - reset all the reactions' mark as false
    %  removeRxn           - delete the reaction from the network
    %  removeSpecies       - remove species from the network
    %  addGlycans          - add the glycans to the reaction network
    %  addGlycan           - add a glycan to the reaction network
    %  findSpeciesByStruct - find the species by the given structure
    %  findSpeciesByName   - find the species by the given name
    %  findComptByName     - find the compartment by the given name
    %  findEnzCompetSubstr - find the competing substrate for the enzyme
    %  setSpeciesID        - reset the ID for each species
    %  setRxnKinetics      - reset the kinetics for each reaction
    %  setSpeciesIsolated  - set the species "isolated"
    %  setComptReactorType - set compartment type
    %  setComptsQflow      - set q flow rate in all compartments
    %  setComptQflow       - set q flow rate in ith compartments
    %  setComptsKt         - set kt in all compartments
    %  setKtinReactor      - set kt in ith compartment
    %  clone               - copy a PATHWAY object
    %
    %
    % See also Rxn,GlycanSpecies,Compt.
    
    % Author: Gang Liu
    % Date Last updated: 5/31/2014
    
    methods
        function setListOfRxns(obj)  
            numRxns = obj.theRxns.length;
            lRxns = cell(numRxns,1);
            for i = 1 : numRxns
                lRxns{i} = obj.theRxns.get(i);
            end
            obj.listofRxns = lRxns;
        end
        
        function setSM(obj)
            numSpecies  = obj.theSpecies.length;
            numRxns     = obj.theRxns.length;
            sm          = zeros(numRxns,numSpecies);
            
            specieslist = cell(numSpecies,1);
            for i=1: numSpecies
              specieslist{i} = obj.theSpecies.get(i);
            end           
            
            for i = 1: numRxns;
                ithRxn = obj.theRxns.get(i);
                for j = 1 : numSpecies
                    jthSpecies = specieslist{j};
                    if(strcmp(ithRxn.reac.id,jthSpecies.id))
                        sm(i,j)=-1;
                    elseif(strcmp(ithRxn.prod.id,jthSpecies.id))
                        sm(i,j)=1;
                    else
                        sm(i,j)=0;
                    end    
                end
            end
            obj.SM = sm;
        end
        
        function setEM(obj)
            numEnzs = obj.theEnzs.length;
            numRxns = obj.theRxns.length;
            em     =  zeros(numRxns,1);
            for i = 1: numRxns;
                ithRxn = obj.theRxns.get(i);
                for j = 1 : numEnzs
                    jthEnz = obj.theEnzs.get(j);
                    if(strcmp(ithRxn.enz.id,jthEnz.id))
                        em(i)=j;
                    end    
                end
            end
            obj.EM = em;
        end
        
        function setComptsQflow(obj,qflow)
            if(~isnumeric(qflow))
                error('MATLAB:GNAT:ERRORINPUTTYPE','WRONG INPUT TYPE');
            end
            for i =1:obj.compartment.length
                setComptQflow(obj,i,qflow);
            end
        end
        
        function setComptQflow(obj,ith,qflow)
            ithreactor = obj.compartment.get(ith).reactortype;
            ithreactor.qflowin  = qflow;
            ithreactor.qflowout = qflow;
        end
        
        function setComptskt(obj,kt)
            if(~isnumeric(kt))
                error('MATLAB:GNAT:ERRORINPUTTYPE','WRONG INPUT TYPE');
            end
            for i =1:obj.compartment.length
                setComptkt(obj,i,kt);
            end
        end
        
        function setComptkt(obj,ith,kt)
            ithreactor     = obj.compartment.get(ith).reactortype;
            ithreactor.ktransport  = kt;
        end
        
        function firstcompt = findFirstCompt(obj)
            if(obj.compartment.length==1)
                firstcompt = obj.compartment.get(1);
            else
                firstcomptindex = [];
                for i = 1 : obj.compartment.length
                    theithcompt = obj.compartment.get(i);
                    if(isempty(theithcompt.priorcompt))
                        firstcomptindex=[firstcomptindex;i];
                    end
                    
                    if(length(firstcomptindex)>1)
                        error('MATLAB:GNAT:ERRORCOMPTSETUP',...
                            'COMPARTMENTS ARE NOT SET UP PROPERLY');
                    end
                end
                
                firstcompt = obj.compartment(firstcomptindex);
            end
        end
        
        function lastcompt = findLastCompt(obj)
            if(obj.compartment.length==1)
                lastcompt = obj.compartment.get(1);
            else
                lastcomptindex = [];
                for i = 1 : obj.compartment.length
                    theithcompt = obj.compartment.get(i);
                    if(isempty(theithcompt.posteriorcompt))
                        lastcomptindex=[lastcomptindex;i];
                    end
                    
                    if(length(lastcomptindex)>1)
                        error('MATLAB:GNAT:ERRORCOMPTSETUP',...
                            'COMPARTMENTS ARE NOT SET UP PROPERLY');
                    end
                end
                
                lastcompt = obj.compartment(lastcomptindex);
            end
        end
    end
    
    properties (Constant)
        VISIT_COLOR_WHITE = 1;
        VISIT_COLOR_GREY = 2;
        VISIT_COLOR_BLACK = 3;
    end
    
    properties
        %THERXNS the reactions in the network.
        %  THERXNS is an array of RXN objects.
        %
        % See also Pathway.
        theRxns
        % list of edges in graph structures
        %In graph theory, they are called list of edges in graph structure
        
        listofRxns
    end
    
    properties
        %THESPECIES the species in the network.
        %  THESPECIES  is an array of GLYCANSPECIES objects
        %
        % See also Pathway.
        theSpecies;
        % In graph theory, they are called list of verticies in graph structure
    end
    
    properties
        %INITGLYCANS the initial structures in the network.
        %  In graph theory, they are called the roots of the graph.
        %  INITGLYCANS property is an array of GLYCANSPECIES objects.
        %
        % See also Pathway.
        initGlycans
    end
    
    properties
        %FINALGLYCANS the terminal structures in the network.
        %  In graph theory, they are called the endings of the graph.
        %  FINALGLYCANS property is an array of GLYCANSPECIES objects.
        %
        % See also Pathway.
        finalGlycans % ending of the graph. biological insights what's the end product
    end
    
    properties
        %THEENZS the enzymes involved in the network
        %     THEENZS is an array of Enzs objects.
        %
        % See also Pathway.
        theEnzs
    end
    
    properties
        %COMPARTMENT the cellular location where the reactions take place
        %    COMPARTMENT is a Compt object
        %
        % See also Pathway.
        compartment;
    end
    
    properties
        % reactortype the type of the reactor for the simulation.
        %
        reactortype
    end
    
    properties
        %NAME the name of the network
        %    NAME property is a character array.
        %
        % See also Pathway.
        name;
    end
    
    properties
        %SM the stoichiometry matrix
        %    SM property is a matrix.
        %
        % See also Pathway.
        SM;
    end
    
    properties
        %EM the enzyme reaction matrix
        %    EM property is a vector.
        %
        % See also Pathway.
        EM;
    end
    
    methods (Static)
        function pathobj=loadmat(matfilename)
            pathstruct=load(matfilename);
            p = fieldnames(pathstruct);
            if(length(p)~=1)
                error('MATLAB:GNAT:WRONGINPUT','wrong variable stored');
            end
            pathobj = pathstruct.(p{1});
            if(isa(pathobj,'Pathway'))
                pathobj.resetjava;
            else
                error('MATLAB:GNAT:WRONGINPUT','wrong variable stored');
            end
        end
    end
    
    methods
        function resetjava(obj)
            for i = 1:obj.theSpecies.length
                thejthspecies = obj.theSpecies.get(i);
                thejthspecies.glycanStruct.resetjava;
            end
        end
    end
    
    methods
        function cloneNetworkinAllCompts(obj,ListOfCompts)  % assume no compts defined yet
            obj.setSpeciesIDIndex;
            obj.setRxnsIDIndex;
            obj.setEnzIDIndex;
            obj.setSM;

            if(ListOfCompts.length==1) % single compt
                for i = 1 : obj.theSpecies.length     
                    obj.theSpecies.get(i).compartment = ListOfCompts.get(1);                    
                end
                
                for  i=1:obj.theEnzs.length
                    obj.theEnzs.get(i).compartment = ListOfCompts.get(1); 
                end
                
                obj.compartment = CellArrayList;
                obj.compartment.add(ListOfCompts.get(1));
             elseif(ListOfCompts.length>1)  % multiple compts
                pathwayclonelib = CellArrayList;
                for i = 1:  ListOfCompts.length-1
                    pathwayclonelib.add(obj.copy(0))
                end
                    
                for i = 1 : obj.theSpecies.length     
                     obj.theSpecies.get(i).compartment = ListOfCompts.get(1);                     
                end
                
                for  i=1:obj.theEnzs.length
                    obj.theEnzs.get(i).compartment = ListOfCompts.get(1); 
                end
                
                obj.compartment = CellArrayList;
                obj.compartment.add(ListOfCompts.get(1));                 
                 
                 for i = 2 : ListOfCompts.length
                      theithcompt  = ListOfCompts.get(i);
                      pathwayclone = pathwayclonelib.get(i-1);  
                      
                      for j = 1 : length(pathwayclone.theSpecies)
                         pathwayclone.theSpecies.get(j).compartment ...
                              = theithcompt; 
                      end
                      
                      for j = 1 : length(pathwayclone.theEnzs)
                          pathwayclone.theEnzs.get(j).compartment ...
                               = theithcompt; 
                      end
                                            
                      obj.addUniqueGlyPath(pathwayclone);
                      obj.addCompt(theithcompt);
                  end
             else
                  return
             end             
        end
        
        function setComptReactorType(obj,reactortypestring)
            %SETCOMPTREACTORTYPE specify reactor type for the
            %  compartment in the pathway
            %
            %See also GlycanNetModel
            obj.reactortype = reactortypestring;
            for i = 1 : obj.getCompartment.length                
                if(strcmpi(obj.reactortype,TypeStringConst.BatchR))
                    obj.getCompartment.get(i).reactortype = BatchReactor;
                elseif(strcmpi(obj.reactortype,TypeStringConst.CSTR))
                    obj.getCompartment.get(i).reactortype = CSTR;
                elseif(strcmpi(obj.reactortype,TypeStringConst.PFR))
                    obj.getCompartment.get(i).reactortype = PFR;
                else
                    error('MATLAB:GNAT:ERRORREACTORTYPE','NOT SUPPORTED REACTOR TYPE');
                end
            end
            
            if(~strcmpi(obj.reactortype,TypeStringConst.BatchR))
                obj.addComptTransportRxns;
            end
        end
        
        function thespeciesinIthCompt = getSpeciesinIthCompt(obj,theithCompt)
            % thespeciesinIthCompt = CellArrayList;
            thespeciesinIthCompt = [];
            count=1;
            for i = 1 : obj.theSpecies.length
                theithspecies = obj.theSpecies.get(i);
                if(strcmpi(theithspecies.compartment.id,theithCompt.id))
                    thespeciesinIthCompt(count,1)=theithspecies;
                    count = count +1;
                end
            end
        end
        
        function thecomptroot = identifyfirstcompt(obj)
            thecompts = obj.getCompartment;
            numroots = 0;
            for i =1: length(thecompts)
                ithcompt = thecompts.get(i);
                if(isempty(ithcompt.priorcompt))
                    thecomptroot = ithcompt;
                    numroots = numroots+1;
                end
            end
            if(numroots>1) || (numroots==0)
                error('MATLAB:GNAT:ERRORCOMPTORDER','THE ORDER OF COMPARTMENTS IS NOT SET UP');
            end
        end
        
        function isComptLinkageValid = verifyLinkageBetwCompts(obj)
            thecomptroot = obj.identifyfirstcompt;
            numcompts = 1;
            theithcompt = thecomptroot;
            while(~isempty(theithcompt.posteriorcompt))
                theithcompt = thecomptroot.posteriorcompt;
                numcompts =  numcompts+1;
            end
            
            isComptLinkageValid = (numcompts==obj.getCompartment.length);
        end
        
        function addComptTransportRxns(obj)
            % first compartment
            thecomptroot = obj.identifyfirstcompt;
            theithcompt  = thecomptroot;
            numadded = obj.addComptInFluxRxns(thecomptroot);
            
            disp(['the number of influx reactions in first compartment is ',num2str(numadded)]);
            while(~isempty(theithcompt.posteriorcompt))
                numadded = obj.addInterComptTransportRxns(theithcompt,...
                    theithcompt.posteriorcompt);
                disp(['the number of inter-compartment reactions is ',num2str(numadded)]);
                theithcompt = theithcompt.posteriorcompt;
            end
            
            numadded = obj.addComptOutFluxRxns(theithcompt);
            disp(['the number of outflux reactions in last compartment is ',num2str(numadded)]);
        end
        
        function numadded = addInterComptTransportRxns(obj,theithcompt,theiplus1compt)
            numadded = 0;
            
            thespeciesinIthCompt      = obj.getSpeciesinIthCompt(theithcompt);
            thespeciesinIplus1thCompt = obj.getSpeciesinIthCompt(theiplus1compt);
            
            if(thespeciesinIthCompt.length>thespeciesinIplus1thCompt.length)
                error('MATLAB:GNAT:SPECIESINPUTERROR',...
                    'THE SPECIES IN TWO COMPARTMENTS ARE NOT DEFINED PROPERLY');
            end
            
            if(thespeciesinIplus1thCompt.length>thespeciesinIthCompt.length)
                error('MATLAB:GNAT:NOTSUPPORTED','FUNCTIONS NOT SUPPORTED YET');
            end
            
            for i = 1 : length(thespeciesinIthCompt)
                ithspeciesinIthCompt = thespeciesinIthCompt(i,1);
                for j = 1 : length(thespeciesinIplus1thCompt)
                    jthspeciesinIplus1thCompt = thespeciesinIplus1thCompt(j,1);
                    if(ithspeciesinIthCompt.glycanStruct.equalStruct(...
                            jthspeciesinIplus1thCompt.glycanStruct))
                        obj.theRxns.add(Rxn(ithspeciesinIthCompt,...
                            jthspeciesinIplus1thCompt,[]));
                        numadded = numadded +1;
                        break
                    end
                end
            end
        end
        
        function numaddedall = addComptInFluxRxns(obj,theithcompt)
            thespeciesinIthCompt = obj.getSpeciesinIthCompt(theithcompt);
            numaddedall = 0;
            for i = 1 : length(thespeciesinIthCompt)
                obj.addInFluxRxn(thespeciesinIthCompt(i,1));
                numaddedall = numaddedall +1;
            end
        end
        
        function addInFluxRxn(obj,thejthspecies)
            obj.theRxns.add(Rxn([],thejthspecies,[]));
        end
        
        function addOutFluxRxn(obj,thejthspecies)
            obj.theRxns.add(Rxn(thejthspecies,[],[]));
        end
        
        function numadded = addComptOutFluxRxns(obj,theithcompt)
            numadded = 0;
            thespeciesinIthCompt = obj.getSpeciesinIthCompt(theithcompt);
            for i = 1 : length(thespeciesinIthCompt)
                obj.addOutFluxRxn(thespeciesinIthCompt(i,1));
                numadded = numadded +1;
            end
        end
    end
    
    methods
        function isemp = isempty(obj)
            %isempty return true if the reaction network contains no reaction
            %  isempty return type is logical
            %
            %  See also Pathway.
            
            isemp = (obj.theRxns.length==0);
        end
        
        function obj = Pathway(varargin)
            %Pathway create a Pathway object.
            %
            %  GP = Pathway creates an empty Pathway object
            %
            %  GP = Pathway(SBML_MODEL,LEVEL,VERSION) converts a SBML_model structure to
            %   a Pathway object and glycan structure format uses GlycoCT format
            %
            %  GP = Pathway(SBML_MODEL,LEVEL,VERSION,GLYCANFORMAT)
            %   specify the structure format of glycans as GLYCANFORMAT
            %   string inputs.
            %
            %  See also Rxn,GlycanSpecies.
            
            % input number check
            % validate the number of input arguments
            %             if(~verLessThan('matlab','7.13'))
            %                 narginchk(0,4);
            %             else
            error(nargchk(0,4,nargin));
            if (nargin == 0)
                obj.theRxns = CellArrayList();
                obj.theSpecies = CellArrayList();
                obj.theEnzs = CellArrayList();
                obj.compartment = CellArrayList();
                obj.name = '';
                return
            elseif ((nargin ==3)||(nargin==4))  %take SBML input
                
                % check arrays of sbml_reaction
                if(~isSBML_Model(varargin{1}))
                    errorReport(mfilename,'IncorrectInputType');
                end
                
                sbmlModel = varargin{1};
                level = varargin{2};
                version = varargin{3};
                
                glycanformat = 'glycoct_xml';
                if(nargin==4)
                    glycanformat = varargin{4};
                end
                
                obj.theRxns = CellArrayList();
                obj.theSpecies = CellArrayList();
                obj.compartment = CellArrayList();
                
                % add compartments to pathway  1st step
                obj.addCompts(sbmlModel.compartment,level,version);
                
                % add array of species in sbmlspecies format 2nd step
                obj.addGlycans(sbmlModel.species,level,version,glycanformat);
                
                %fprintf(1,'number of species is : %i \n',numSpeices);
                
                % add reactions to pathway %3rd step
                obj.addRxns(sbmlModel.reaction);
                
                %fprintf(1,'number of reactions is: %i\n',numRxns);
                
                % add enzyme type
                % currently empty
            end
        end % end constructor
    end
    
    methods
        function isGlycan= find(obj,species)
            % find check if the species exist in the species list.
            %
            %  isGlycan =  find(species) check if the species exisits
            %    in the list
            %
            %  See also Rxn,GlycanSpecies
            isGlycan=false;
            for i = 1 : length(obj.theSpecies)
                if(obj.theSpecies.get(i)==species)
                    isGlycan=true;
                    return
                end
            end
        end
        
        function glycanIndex = locspecies(obj,species)
            % locspecies return location of the species  in the species list.
            %
            %  isGlycan =  locspecies(species) check if the species exisits
            %    in the list
            %
            %  See also Rxn,GlycanSpecies
            glycanIndex=[];
            for i = 1 : length(obj.theSpecies)
                if(obj.theSpecies.get(i)==species)
                    glycanIndex = [glycanIndex;i];
                end
            end
        end
        
        
        
        % set all compartments based on the sbml_compartment
        function numCompts = addCompts(obj,varargin)  % cellarraylist
            % addCompts add an array of Compt objects to the network.
            %
            %  NUMCOMPTSADDED =  addCompts(LISTOFCOMPTS) takes an
            %   CellArrayList object and adds a series of Compt
            %   objects to the network.
            %
            %  NUMCOMPTSADDED =  addCompts(SBMLCOMPTS,LEVEL,VERSION) takes an
            %   SBML Compartment structure, SBML-level and version, converts
            %   to a number of Compt objects, and add them to the network.
            %
            %  See also Rxn,GlycanSpecies.
            
            numCompts=0;
            if (isa(varargin{1},'CellArrayList'))
                comptsToAdd = varargin{1};
                for i=1:  length(comptsToAdd);
                    isAdded=obj.addCompt(comptsToAdd.get(i));
                    if(isAdded)
                        numCompts=numCompts+1;
                    end
                end
            elseif (isa(varargin{1},'struct'))  %SBML-compartment
                level = varargin{2};
                version = varargin{3};
                comptsToAdd = varargin{1};
                
                for i=1:  length(comptsToAdd);
                    isAdded =obj.addCompt(comptsToAdd(i),level,version);
                    if(isAdded)
                        numCompts=numCompts+1;
                    end
                end
            end
        end
        
        % add a compartment to pathway
        function isAdded = addCompt(obj,varargin)
            % addCompt add a Compt object to the network
            %
            %  ISCOMPTADDED =  addCompt(GCOMPT) takes an
            %   Compt object and adds it to the network
            %
            %  ISCOMPTADDED =  addCompt(SBMLCOMPARTMENT,LEVEL,VERSION) takes an
            %   SBML_COMPARTMENT structure, SBML-level and verions, and
            %   adds a series of Compt objects to the network.
            %
            %  See also Rxn,Compartment.
            
            % add exception control
            if((length(varargin)==3))
                sbmlCompt = varargin{1};
                level = varargin{2};
                version = varargin{3};
                comptToAdd = Compt(sbmlCompt,level,version); % convert from sbml-compartment to compt
            elseif(length(varargin)==1)
                comptToAdd = varargin{1};
                if(~isa(comptToAdd,'Compt'))
                    errorReport(mfilename,'IncorrectInputType');
                end
            end
            
            if(~obj.compartment.contains(comptToAdd))
                %  fprintf('size of original glycan array is: %i \n ', length(obj.theSpecies));
                obj.compartment.add(comptToAdd);
                % fprintf('size of new glycan array: %i \n ', length(obj.theSpecies));
                isAdded = true;
            else
                % fprintf('species name is %s  \n',sbml_species.name);
                % fprintf('size of glycan array remains at : %i \n ', length(obj.theSpecies));
                isAdded = false;
            end
        end
        
        % set all glycans based on the sbml_species
        function numGlycans = addGlycans(obj,varargin)  % cellarraylist
            % addGlycans add an array of GlycanSpecies objects to the network.
            %
            %  NUMGLYCANADDED =  addGlycans(LISTOFGLYCANS) takes an
            %   CellArrayList object and adds a series of GlycanSpecies
            %    objects to the network.
            %
            %  NUMGLYCANADDED =  addGlycans(SBMLSPECIES,LEVEL,VERSION) takes an
            %   SBML species object, SBML-level and version, and adds a series
            %   of GlycanSpecies objects to the network.
            %
            %  See also Rxn,GlycanSpecies.
            
            numGlycans=0;
            if (isa(varargin{1},'CellArrayList'))
                glycansToAdd = varargin{1};
                for i=1:  length(glycansToAdd);
                    isAdded=obj.addGlycan(glycansToAdd.get(i));
                    if(isAdded)
                        numGlycans=numGlycans+1;
                    end
                end
            elseif (isa(varargin{1},'struct'))  %SBML-species
                level = varargin{2};
                version = varargin{3};
                glycanformat = 'glycoct_xml';
                if(nargin==4)
                    glycanformat = varargin{4};
                end
                glycansToAdd = varargin{1};
                for i=1:  length(glycansToAdd);
                    isAdded =obj.addGlycan(glycansToAdd(i),level,version,glycanformat);
                    if(isAdded)
                        numGlycans=numGlycans+1;
                    end
                end
            end
            %numGlycans
            %  fprintf(1,'number of glycan added: %i',numGlycans);
        end
        
        %
        function locs=findsameStructGlycan(obj,glycanSpeciesObj)
            %findsameStructGlycan find a GlycanSpecies object in the species
            %list sharing the same structure
            %
            % locs =  findsameStructGlycan(glycanSpeciesObj) find the
            %   position of species sharing the same glyan structure. If not
            %   found, return 0.
            %
            %  See also Rxn,GlycanSpecies.
            locs = 0;
%             try
%                 b=glycanSpeciesObj.glycanStruct;
%             catch error
%                 disp('error');
%             end
            for i = 1 : length(obj.theSpecies)
                % a = obj.theSpecies.get(i).glycanStruct;
                if(obj.theSpecies.get(i).glycanStruct.equalStruct(...
                        glycanSpeciesObj.glycanStruct))
                    locs=i;
                    return  % assume only one same structure in the list
                end
            end
        end
        
        function locs=findsameStructRxn(obj,rxnObj)
            %findsameStructRxn find a reaction object in the reaction
            %list sharing the same reactant and product structure
            %
            % locs =  findsameStructRxn(glycanRxnObj) find the
            % position of reaction sharing the same glyan reaction. If not
            % found, return 0.
            %
            %  See also Rxn,GlycanSpecies.
            locs = 0;
            for i = 1 : length(obj.theRxns)
                %isSameRxn = false;
                theithRxn = obj.theRxns.get(i);
                react = theithRxn.getReactant;
                prod = theithRxn.getProduct;
                
                if(equalStruct(react.glycanStruct, ....
                        rxnObj.getReactant.glycanStruct)...
                        &&equalStruct(prod.glycanStruct,...
                        rxnObj.getProduct.glycanStruct))
                    locs=i;
                    return  % assume only one same structure in the list
                end
            end
        end
        
        
        % add glycan to pathway
        function isAdded = addGlycan(obj,varargin)
            %addGlycan add a GlycanSpecies object to the network
            %
            % ISGLYCANADDED =  addGlycan(GLYCANSPEC) takes an
            %  GlycanSpecies object and adds it to the network
            %
            % ISGLYCANADDED =  addGlycan(SBMLSPECIES,LEVEL,VERSION) takes an
            %  SBML Species structure, SBML-level and version, and adds a
            %  number of GlycanSpecies objects to the network.
            %
            %  See also Rxn,GlycanSpecies.
            
            % add exception control
            if((length(varargin)==3)||length(varargin)==4)
                sbml_species = varargin{1};
                level = varargin{2};
                version = varargin{3};
                if(isSBML_Species(sbml_species,level,version))
                    glycanToAdd = GlycanSpecies(sbml_species,level,version); % convert from sbml-speceis to theSpecies
                    
                    % set compartment for glycan
                    theComptName = Species_getCompartment(sbml_species);
                    theCompt = obj.findComptByName(theComptName);
                    glycanToAdd.setCompartment(theCompt);
                else
                    error('MATLAB:GNAT:WRONGSBMLVER','sbml version is not right');
                end
            elseif(length(varargin)==1)
                glycanToAdd = varargin{1};
                if(~isa(glycanToAdd,'GlycanSpecies'))
                    error('MATLAB:GNAT:IncorrectInputType','Wrong Input Type');
                end
            end
            
            if(~obj.theSpecies.contains(glycanToAdd)) % check if there is equal object
                % fprintf('size of original glycan array is: %i \n ', length(obj.theSpecies));
                locs=findsameStructGlycan(obj,glycanToAdd);
                % fprintf('find same structure pos: %i \n',locs);
                if(locs==0)
                    isAdded=true;
                else %if(isprop(obj.theSpecies.get(locs),'compartment'))
                    % fprintf('find same structure : %i \n ', locs);
                    comptToAdd  = glycanToAdd.getCompartment;
                    comptCurrent = obj.theSpecies.get(locs).getCompartment;
                    if(comptToAdd==comptCurrent)
                        isAdded =false;
                    else
                        isAdded = true;
                    end
%                 else
%                     isAdded=false;
                end
            end
            
            if(isAdded)
                % fprintf('add structure \n ');
                obj.theSpecies.add(glycanToAdd);
                if(~isempty(glycanToAdd.getCompartment)) 
                    obj.addCompt(glycanToAdd.getCompartment);
                end
            else
                % fprintf('species name is %s  \n',sbml_species.name);
                % fprintf('size of glycan array remains at : %i \n ', length(obj.theSpecies));
                %isAdded = false;
            end
        end
        
        % add a new pathway to the network
        function [numSpeciesAdded,numRxnsAdded,numEnzAdded] =  addGlyPath(obj,pathtoadd)
            % addGlyPath add a pathway to the network
            %
            %  NUMRXNSADDED =  addGlyPath(obj,pathtoadd) takes an
            %   pathway object and adds a number of Rxn and Species
            %   objects in the list to the network
            %
            %  See also Pathway.
            
            rxnsToAdd          = pathtoadd.theRxns;
            speciesToAdd       = pathtoadd.theSpecies;
            enzsToAdd          = pathtoadd.theEnzs;
            
            numSpeciesAdded    = obj.addGlycans(speciesToAdd); % check structure and object equality
            numRxnsAdded       = obj.addRxns(rxnsToAdd);
            numEnzAdded        = obj.addEnzs(enzsToAdd);            
        end
        
        % add a different pathway to the network
        function [numSpeciesAdded,numRxnsAdded,numEnzAdded] = addUniqueGlyPath(obj,pathtoadd)
            % addGlyPath add a pathway to the network
            %
            %  NUMRXNSADDED =  addUniqueGlyPath(obj,pathtoadd) takes an
            %   pathway object and adds a number of Rxn and Species
            %   objects in the list to the network
            %
            %  See also Pathway.
            
            rxnsToAdd          = pathtoadd.theRxns;
            speciesToAdd       = pathtoadd.theSpecies;
            enzsToAdd          = pathtoadd.theEnzs;
            
            for j = 1 : speciesToAdd.length
               obj.theSpecies.add(speciesToAdd.get(j));               
            end
            
            for j = 1 : rxnsToAdd.length
               obj.theRxns.add(rxnsToAdd.get(j));
            end
            
            for j = 1 : enzsToAdd.length
               obj.theEnzs.add(enzsToAdd.get(j));
            end    
            
            numSpeciesAdded =  speciesToAdd.length;
            numRxnsAdded    =  rxnsToAdd.length;
            numEnzAdded     =  enzsToAdd.length;
        end
        
        % add a new pathway to the network based on the structure
        function [newGlycans,newRxns] =  addGlyPathByStruct(obj,pathtoadd)
            % addGlyPath add a pathway to the network
            %
            %  NUMRXNSADDED =  addGlyPath(obj,pathtoadd) takes an
            %   pathway object and adds a number of Rxn and Species
            %   objects in the list to the network
            %
            %  See also Pathway.
            
            rxnsToAdd       = pathtoadd.theRxns;
            speciesToAdd    = pathtoadd.theSpecies;
            enzsToAdd       = pathtoadd.theEnzs;
            
            newGlycans   = obj.addGlycansByStruct(speciesToAdd); % check structure and object equality
            newRxns      = obj.addRxnsByStruct(rxnsToAdd);
           %  newEnzs      = obj.addEnzs(enzsToAdd);
            %obj.addRxns(rxnsToAdd);
        end
        
        % add reactions to pathway
        function newRxns =  addRxnsByStruct(obj,varargin)
            % addRxnsByStruct add an array of Rxn object(s) to the network
            %
            %  NUMRXNSADDED =  addRxnsByStruct(LISTOFRXNS) takes an
            %   CellArrayList object and adds a number of Rxn
            %   objects in the list to the network
            %
            %
            %  See also Pathway.
            
            newRxns = CellArrayList;
            if(isa(varargin{1},'CellArrayList'))
                rxnsToAdd = varargin{1};
                for i=1:rxnsToAdd.length
                    %fprintf(1,' ith reaction %f\n',i);
                    theithrxn = rxnsToAdd.get(i);
                    isAdded = obj.addRxnByStruct(theithrxn);
                    if(isAdded)
                        newRxns.add(theithrxn);
                    end
                end
            end
        end
        
        % add reaction to pathway by structure
        function isAdded = addRxnByStruct(obj,varargin)
            % addRxnByStruct add a Rxn object to the network
            %
            %   ISRXNADDED =  addRxnByStruct(THERXN) takes the RXN object as the input to create
            %     a reaction and adds it to the network
            %
            %   ISRXNADDED =  addRxnByStruct(THEREAC, THERPODUCT) takes the GlycanSpecies objects
            %    THEREACT and THERPODUCT as the input to create a reaction and adds it to
            %    the network
            %
            %   ISRXNADDED =  addRxnByStruct(THEREAC,THERPODUCT,THEENZ) takes the
            %      GlycanSpecies objects THEREACT and THERPODUCT and the
            %      Enz object THEENZ as the inputs to create a reaction and adds
            %      it to the network.
            %
            %  See also Pathway.
            
            % narginchk(1,3);
            
            if(length(varargin)==1)
                if isa(varargin{1},'Rxn')
                    rxnToAdd = varargin{1};
                    reactant = rxnToAdd.getReactant;
                    product  = rxnToAdd.getProduct;
                elseif(length(varargin)==2)
                    reactant = varargin{1};
                    product  = varargin{2};
                    rxnToAdd = Rxn(reactant,product);
                elseif(length(varargin)==3)
                    reactant   = varargin{1};
                    product    = varargin{2};
                    enzyme     =  varargin{3};
                    rxnToAdd   = Rxn(reactant,product,enzyme);
                end
                
                if(~isempty(rxnToAdd.rxnkinetics))
                    rxnkinetics = rxnToAdd.rxnkinetics;
                else
                    rxnkinetics =[];
                end
                
                reactlocs=obj.findsameStructGlycan(reactant);
                prodlocs =obj.findsameStructGlycan(product);
                
                if((~isempty(reactant))&&(reactlocs==0))
                   obj.theSpecies.add(reactant); 
                   reactlocs = length(obj.theSpecies);
                end
                
                if((~isempty(product) &&(prodlocs==0)))
                %    errorReport(mfilename,'product species is not in the network');
                   obj.theSpecies.add(product); 
                   prodlocs = length(obj.theSpecies);
                end
                
                % check if reactant or product is already included in the
                % species list
                if(reactlocs>0)  % if structure is same
                    reactant  = obj.theSpecies.get(reactlocs);
                end
                
                if(prodlocs>0)  % if structure is same check compartment
                     product  = obj.theSpecies.get(prodlocs);
                end
                
                if(isempty(rxnToAdd.enz))
                    rxnToAdd = Rxn(reactant,product);
                else
                    rxnToAdd = Rxn(reactant,product,rxnToAdd.enz);
                end
                
                rxnlocs = obj.findsameStructRxn(rxnToAdd);
                
                if(rxnlocs==0)
                    if(~isempty(rxnkinetics))
                        rxnToAdd.rxnkinetics= rxnkinetics;
                    end
                    
                    if(~isempty(rxnToAdd.enz))
                        sameEnztofind = obj.findsameECNOEnz(rxnToAdd.enz);
                        if(isempty(sameEnztofind))
                            obj.theRxns.add(rxnToAdd);
                            obj.theEnzs.add(rxnToAdd.enz);
                        else
                            rxnToAdd.enz = sameEnztofind;
                            obj.theRxns.add(rxnToAdd);
                        end
                    else
                        obj.theRxns.add(rxnToAdd);     
                    end
                    isAdded = true;
                    
                    
                    if(~isempty(reactant))
                        reactant.addRxn(rxnToAdd);
                    end
                    
                    if(~isempty(product))
                        product.addRxn(rxnToAdd);
                    end
                else
                    isAdded = false;
                end
            end
        end
        
        % set all glycans based on the structure
        function newGlycans = addGlycansByStruct(obj,varargin)
            % addGlycansByStruct add an array of GlycanSpecies objects to the existing network.
            %
            %  NUMGLYCANADDED =  addGlycansByStruct(LISTOFGLYCANS) takes an
            %   CellArrayList object and adds a series of GlycanSpecies
            %    objects to the network.
            %
            %
            %  See also Rxn,GlycanSpecies.
            
            newGlycans = CellArrayList;
            if (isa(varargin{1},'CellArrayList'))
                glycansToAdd = varargin{1};
                for i=1:  length(glycansToAdd);
                    theithglycan = glycansToAdd.get(i);
                    isAdded=obj.addGlycanByStruct(glycansToAdd.get(i));
                    if(isAdded)
                        newGlycans.add(theithglycan);
                    end
                end
            end
        end
        
        % add glycan to pathway by structure
        function isAdded = addGlycanByStruct(obj,varargin)
            %addGlycanByStruct add a GlycanSpecies object to the network
            % by structure
            %
            % ISGLYCANADDED =  addGlycanByStruct(GLYCANSPEC) takes an
            %  GlycanSpecies object and adds it to the network
            %
            %
            %  See also Rxn,GlycanSpecies.
            
            % add exception control
            if(length(varargin)==1)
                glycanToAdd = varargin{1};
                %                 if(~isa(glycanToAdd,'GlycanSpecies'))
                %                     error('MATLAB:GNAT:IncorrectInputType','Incorrect Input Type');
                %                 end
            end
            
            locs=findsameStructGlycan(obj,glycanToAdd);
            if(locs==0)
                isAdded=true;
            else
                isAdded=false;
            end
            
            if(isAdded)
                obj.theSpecies.add(glycanToAdd);
                if(isprop(glycanToAdd,'Compartment') ...
                        &&(~isempty(glycanToAdd.getCompartment)) )
                    obj.addCompt(glycanToAdd.getCompartment);
                end
            else
                % fprintf('species name is %s  \n',sbml_species.name);
                % fprintf('size of glycan array remains at : %i \n ', length(obj.theSpecies));
                %isAdded = false;
            end
        end
        
        % add an array of Enzymes to pathway
        function newEnzs =  addEnzs(obj,varargin)
            % addRxns add an array of Enz object(s) to the network
            %
            %  numEnzs =  addEnzs(LISTOFENZS) takes an
            %   CellArrayList object and adds a number of Enz
            %   objects in the list to the network
            %
            %  See also Pathway.
            
            newEnzs = CellArrayList;
            if(isa(varargin{1},'CellArrayList'))
                enzsToAdd = varargin{1};
                for i=1:enzsToAdd.length
                    %fprintf(1,' ith reaction %f\n',i);
                    theIthEnz = enzsToAdd.get(i);
                    isAdded = obj.addEnz(theIthEnz);
                    if(isAdded)
                        newEnzs.add(theIthEnz);
                    end
                end
            end
        end
        
        % add an  Enzymes to the pathway
        function numEnzs =  addEnz(obj,varargin)
            % addEnz add an Enz object to the network
            %
            %  numEnzs =  addEnz(enz) takes an
            %   Enz object to the network
            %
            %  See also Pathway.
            
            numEnzs = 0;
            enztoadd = varargin{1};
            if(~obj.theEnzs.contains(enztoadd))
                obj.theEnzs.add(enztoadd);
            end
        end
        
        % add reactions to pathway
        function numRxns =  addRxns(obj,varargin)
            % addRxns add an array of Rxn object(s) to the network
            %
            %  NUMRXNSADDED =  addRxns(LISTOFRXNS) takes an
            %   CellArrayList object and adds a number of Rxn
            %   objects in the list to the network
            %
            %  NUMRXNSADDED =  addRxns(SBML_REACTION,LEVEL,VERSION) takes an
            %   SBML reaction structure, SBML-level and version, and adds a
            %   series of Rxn objects to the network.
            %
            %  See also Pathway.
            
            numRxns = 0;
            if(isa(varargin{1},'CellArrayList'))
                rxnsToAdd = varargin{1};
                for i=1:rxnsToAdd.length
                    %fprintf(1,' ith reaction %f\n',i);
                    isAdded = obj.addRxn(rxnsToAdd.get(i));
                    if(isAdded)
                        numRxns = numRxns+1;
                    end
                end
            elseif(isa(varargin{1},'struct'))  %ccheck valid of sbml-reaction
                rxnsToAdd = varargin{1};
                for i=1:length(rxnsToAdd)
                    % fprintf(1,' ith reaction %f\n',i);
                    isAdded = obj.addRxn(rxnsToAdd(1,i));
                    if(isAdded)
                        numRxns = numRxns+1;
                    else
                        disp('reaction is already included in the pathway');
                    end
                end
            end
        end
        
        % add reaction to pathway
        function isAdded = addRxn(obj,varargin)
            % addRxn add a Rxn object to the network
            %
            %   ISRXNADDED =  addRxn(SBMLRXN) takes a SBML reaction structure SBMLRXN
            %    as an input and adds it to the network
            %
            %   ISRXNADDED =  addRxn(THEREAC, THERPODUCT) takes the GlycanSpecies objects
            %    THEREACT and THERPODUCT as an inputs to create a reaction and adds it to
            %    the network
            %
            %   ISRXNADDED =  addRxn(THEREAC,THERPODUCT,THEENZ) takes the
            %      GlycanSpecies objects THEREACT and THERPODUCT and the
            %      Enz object THEENZ as the inputs to create a reaction and adds
            %      it to the network
            %
            %  See also Pathway.
            
            % narginchk(1,3);
            
            if(length(varargin)==1)
                if isa(varargin{1},'Rxn')
                    rxnToAdd = varargin{1};
                    reactant = rxnToAdd.getReactant;
                    product  = rxnToAdd.getProduct;
                elseif ( ...% isValid(varargin{1}) &&...
                        isSBML_Reaction(varargin{1},GetLevel(varargin{1}))...
                        )
                    sbmlRxn =varargin{1};
                    %[level,version]=GetLevelVersion(sbmlRxn);
                    
                    if(~isempty(sbmlRxn.reactant))
                        reactant = obj.findSpeciesByID(sbmlRxn.reactant(1).species);
                    else
                        reactant = [];
                    end
                    
                    if(~isempty(sbmlRxn.product))
                        product  = obj.findSpeciesByID(sbmlRxn.product(1).species);
                    else
                        product = [];
                    end
                    
                    if((isempty(reactant))&&(isempty(product)))
                        errorReport(mfilename,'wrong reaction assignment')
                    end
                    
                    if(isempty(sbmlRxn.name))
                        if~isempty(reactant)
                            reactantname =reactant.name;
                        else
                            reactantname = 'null';
                        end
                        
                        if~isempty(product)
                            productname =product.name;
                        else
                            productname = 'null';
                        end
                        
                        rxnName = [reactantname '_to_' productname];
                    else
                        rxnName = sbmlRxn.name;
                    end
                    
                    if(isempty(rxnName))
                        rxnName ='unknown';
                    end
                    
                    rxnToAdd = Rxn(reactant,product,rxnName);
                end
            elseif(length(varargin)==2)
                reactant = varargin{1};
                product  = varargin{2};
                rxnToAdd = Rxn(reactant,product);
            elseif(length(varargin)==3)
                reactant     = varargin{1};
                product     = varargin{2};
                enzyme     =  varargin{3};
                rxnToAdd = Rxn(reactant,product,enzyme);
            end
            
            if(~isempty(reactant))
                reactlocs=obj.findsameStructGlycan(reactant);
            end
            
            if(~isempty(product))
                prodlocs =obj.findsameStructGlycan(product);
            end
            
            if((~isempty(reactant)) &&(~obj.theSpecies.contains(reactant))...
                    && (reactlocs==0))
                errorReport(mfilename,'reactant species is not in the network');
                return
            end
            
            if((~isempty(product) &&~obj.theSpecies.contains(product))...
                    && (prodlocs==0))
                errorReport(mfilename,'product species is not in the network');
            end
            
            % check if reactant or product is already included in the
            % species list
            if(~isempty(reactant))
                isReactsInList = obj.theSpecies.contains(reactant);
                if(isReactsInList)
                    ithpos =  locationsOf(obj.theSpecies,reactant);
                    reactant =obj.theSpecies.get(ithpos);
                elseif(reactlocs>0)  % if structure is same check compartment
                    reactant =obj.theSpecies.get(reactlocs);                   
                end
            end
            
            if(~isempty(product))
                isProdInlist = obj.theSpecies.contains(product);
                if(isProdInlist)
                    ithpos =  locationsOf(obj.theSpecies,product);
                    product =obj.theSpecies.get(ithpos);
                elseif(prodlocs>0)  % if structure is same check compartment
                    product =obj.theSpecies.get(prodlocs);
                end
            end
            
            if(~isempty(rxnToAdd.rxnkinetics))
                rxnkineticscopy = rxnToAdd.rxnkinetics;
            else
                rxnkineticscopy = [];
            end
            
            if(isempty(rxnToAdd.enz))
                rxnToAdd = Rxn(reactant,product);
            else
                rxnToAdd = Rxn(reactant,product,rxnToAdd.enz);
            end
            
            if(~isempty(rxnkineticscopy))
                rxnToAdd.rxnkinetics = rxnkineticscopy;
            end
            
            if(~obj.theRxns.contains(rxnToAdd))
                obj.theRxns.add(rxnToAdd);
                isAdded = true;
                if(~isempty(reactant))
                    reactant.addRxn(rxnToAdd);
                end
                if(~isempty(product))
                    product.addRxn(rxnToAdd);
                end
            else
                isAdded = false;
            end
        end
        
        % remove glycan from the pathway
        function isRemoved = removeSpecies(obj,speciesToRemove)
            % removeSpecies remove the species from the network
            %
            %  isRemoved = removeSpecies(THEPATHWAY,SPECIESTOREMOVE) removes
            %   the GlycanSpecies speciesToRemove from the Pathway Object
            %   THEPATHWAY.
            %
            % See also Pathway
            
            if( ~obj.theSpecies.contains(speciesToRemove))
                isRemoved = false;
                return
            else
                speciesLoc =  obj.theSpecies.locationsOf(speciesToRemove);
            end
            
            obj.theSpecies.remove(speciesLoc);
            
            % check the rxn list which contains the species to remove and
            % remove it
            for i=1:obj.theRxns.length
                theithRxn = obj.theRxns.get(i);
                theReactant = theithRxn.getReactant;
                if(isequal(speciesToRemove,theReactant))
                    obj.theRxns.remove(i);
                    continue;
                end
                
                theProduct = theithRxn.getProduct;
                if(isequal(speciesToRemove,theProduct))
                    obj.theRxns.remove(i);
                end
            end
            
            %             %remove the reaction associated with species from the network
            %             %being a reactant
            %             for i=1:speciesToRemove.getNumReacRxns
            %                 reacRxn = speciesToRemove.getReacRxn(i);
            %                 speciesToRemove.removeRxn(reacRxn);
            %                 prod = reacRxn.getProduct;
            %                 prod.removeRxn(reacRxn);
            %                 reacRxnLoc = obj.theRxns.locationsOf(reacRxn);
            %                 obj.theRxns.remove(reacRxnLoc);
            %             end
            %
            %             %remove the reaction associated with species from the network
            %             %being a product
            %             for i=1:speciesToRemove.getNumProdRxns
            %                 prodRxn = speciesToRemove.getProdRxn(i);
            %                 speciesToRemove.remove(prodRxn);
            %                 reac = prodRxn.getReac;
            %                 reac.remove(prodRxn);
            %                 obj.theRxns.remove(prodRxn);
            %             end
            
            isRemoved = true;
        end
        
        % remove rxn from the pathway
        function isRemoved = removeRxn(obj,reac,prod)
            % removeRxn remove a reaction from the network
            %   STATUS =  REMOVERXN(REACTANT,PRODUCT) remove a reaction
            %    where the reactant is REACTANT and the product is
            %    PRODUCT.
            %
            % See also Pathway
            rxn = reac.findRxn(prod);
            if(isempty(rxn))
                isRemoved = false;
                return;
            else
                reac.remove(rxn);
                prod.remove(rxn);
                obj.theRxns.remove(rxn);
                isRemoved=true;
            end
        end
        
        function isequal = equalPath(obj,obj2)
            if(obj.theSpecies.length~=obj2.theSpecies.length)
                isequal = false;
                return
            end
            
            for i = 1 : obj.theSpecies.length
                ithspecies         =  obj.theSpecies.get(i);
                isspeciesfound = obj2.findsameStructGlycan(ithspecies);
                if(~isspeciesfound)
                    isequal=false;
                    return
                end
            end
            
            isequal = true;
        end
        
        function clearSpeciesMark(obj)
            % clearSpeciesMark reset all marks as false
            %
            % See also Pathway
            for i=1:obj.theSpecies.length
                obj.theSpecies.get(i).clearMark;
                obj.theSpecies.get(i).clearMarkState;
            end
        end
        
        function clearRxnsMark(obj)
            % clearRxnsMark remove all reactions'a marks from the network
            %
            % See also Pathway
            for i=1:obj.theRxns.length
                obj.theRxns.get(i).clearMark;
            end
        end
    end
    
    methods
    end
    
    methods
        function speciesToFindByName = findSpeciesByName(obj,specName)
            % findSpeciesByName search species in the network by its name
            %
            %   SP = FINDSPECIESBYNAME(PATHWAY,SPNAME) takes SPNAME as an
            %    input and returns the reference to the species which has the same
            %    name. If no same name is found, SP is returns as an empty
            %    value
            %
            % See also Pathway.
            speciesToFindByName =[];
            for i=1: obj.theSpecies.length
                theithSpecies = obj.theSpecies.get(i);
                if(strcmpi(theithSpecies.getName,specName))
                    speciesToFindByName = theithSpecies;
                    return;
                end
            end
        end
        
        function sameECNOEnz = findEnzByECNO(obj,ecno)           
            sameECNOEnz =[];
            for i=1: obj.theEnzs.length
                theithEnz = obj.theEnzs.get(i);
                if(isequal(theithEnz.ecno,ecno))
                    sameECNOEnz = theithEnz;
                    return;
                end
            end
        end
        
        function thesameEnz  = findsameECNOEnz(obj,theEnztofind)
            thesameEnz = obj.findEnzByECNO(theEnztofind.ecno);
        end
        
        function comptToFindByName = findComptByName(obj,comptName)
            % findComptByName search compartment in the network by its name
            %
            %   COMPT = FINDCOMPTBYNAME(PATHWAY,COMPTNAME) takes
            %    COMPTNAME as an input and returns the reference to the
            %    compartment which has the same name. If no same name is found,
            %     COMPT is returns as an empty value.
            %
            % See also Pathway.
            comptToFindByName =[];
            for i=1: obj.compartment.length
                theithcompartment = obj.compartment.get(i);
                if (strcmpi(theithcompartment.getName,comptName))
                    comptToFindByName = theithcompartment;
                    return;
                end
            end
        end
        
        % search the species by glycan name
        function speciesToFindByID = findSpeciesByID(obj,id)
            % findSpeciesByID search species in the network by its id
            %
            %   SP = FINDSPECIESBYID(PATHWAY,SPID) takes SPID as an
            %    input and returns the reference to the species which has the same
            %    ID. If no same name is found, SP is returns as an empty
            %    value
            %
            % See also Pathway.
            speciesToFindByID =[];
            for i=1: obj.theSpecies.length
                theithSpecies = obj.theSpecies.get(i);
                if(strcmpi(theithSpecies.getID,id))
                    speciesToFindByID = theithSpecies;
                    return;
                end
            end
        end
        
        function isStructinPath = isStructinPath(obj,spec)
            isStructinPath = ~isempty(obj.findSpeciesByStruct(spec));
        end
        
        function savepathway(obj,filename)
            pathtosave = obj.clone;
            for i=1: pathtosave.getNSpecies
                ithspecies=pathtosave.theSpecies.get(i);
                ithspecies.glycanStruct.glycanjava=[];
            end
            save(filename,'pathtosave');
        end
        
        % search the species by
        function [speciesToFindByStruct,varargout] = findSpeciesByStruct(obj,spec,varargin)
            % findSpeciesByStruct search species in the network by its structure
            %
            %     SP = FINDSPECIESBYSTRUCT(PATHWAY,SPSTRUCT) takes SPSTRUCT as an
            %    input and returns the reference to the species which has the same
            %    structure. If no same structure is found, SP is returns as empty
            %    value
            %
            % See also Pathway.
            
            if(length(varargin)==1)
                comptname=varargin{1};
            else
                comptname=[];
            end   
            
            speciesToFindByStruct =[];
            for i = 1: obj.theSpecies.length
                ithspecies = obj.theSpecies.get(i);
                if(spec.glycanStruct.equalStruct(ithspecies.glycanStruct))
                    if(isempty(comptname)) 
                        speciesToFindByStruct = ithspecies;
                        break;
                    elseif(strcmpi(ithspecies.compartment.name,comptname))
                        speciesToFindByStruct = ithspecies;
                        break;
                    end
                end
            end
            
            if(length(nargout)==1)
               varargout{1}=i;
            end
        end
        
        %set ID for each species
        function setSpeciesIDIndex(obj)
            for i=1: obj.theSpecies.length
                obj.theSpecies.get(i).id = ['G' int2str(i)];
            end
        end
        
        %set ID for each species
        function setRxnsIDIndex(obj)
            for i=1: obj.theRxns.length
                obj.theRxns.get(i).id = ['GR' int2str(i)];
            end
        end
        
        %set ID for each enzyme
        function setEnzIDIndex(obj)
            for i = 1 : obj.theEnzs.length
                obj.theEnzs.get(i).id = ['E' int2str(i)];              
            end
        end
        
        %set initial concentration zero for all species
         function setSpeciesInitConcZero(obj)
            for i=1: obj.theSpecies.length
                obj.theSpecies.get(i).initConc = 0;
            end
        end
        
        
        function setAllRxnsKinetics(obj)
            obj.setEnzRxnsKinetics
            obj.setTransportRxnKinetics;
        end
        
        
        %set reaction kinetic for each enzymatic reaction
        function setEnzRxnsKinetics(obj,varargin)
            %set up Vm, Km
            for i = 1:obj.theRxns.length
                theithRxn = obj.theRxns.get(i);
                if(~isempty(theithRxn.enz))
                    theithRxn.setrxnkinetics;
                end
            end
            
            %set up inhibition parameter if exist
            for i = 1 : obj.theRxns.length
                theithRxn = obj.theRxns.get(i);
                if(~isempty(theithRxn.enz))
                    theCompetingRxns = findEnzCompetingRxns(obj,theithRxn) ;
                    theithRxn.rxnkinetics.setRxnKineticsInhibitor(theCompetingRxns)
                end
            end
            
            %set up formula
            for i = 1 : obj.theRxns.length
                theithRxn = obj.theRxns.get(i);
                if(~isempty(theithRxn.enz))
                    if(isa(theithRxn.enz.enzkinetics,'MMenKinetics'))
                        theithRxn.rxnkinetics.setMathFormula(theithRxn.enz.id,...
                            theithRxn.id);
                    elseif(isa(theithRxn.enz.enzkinetics,'BiBiKinetics'))  % if the enzyme is BiBiKinetics
                        if(length(varargin)==1)
                            sugarnucldb = varargin{1};                            
                            
                            if(strcmpi(theithRxn.enz.donor,'UDP-N-acetyl-D-glucosamine'))
                                sugarnucl = 'UDP_GlcNAc';
                                nucl = 'UDP';
                            elseif(strcmpi(theithRxn.enz.donor,'UDP-alpha-D-galactose'))
                                sugarnucl = 'UDP_GaL';
                                nucl = 'UDP';
                            elseif(strcmpi(theithRxn.enz.donor,'CMP-N-acetylneuraminate'))% CMP-Sialic Acid
                                sugarnucl = 'CMP-NeuAc';
                                nucl = 'CMP';
                            elseif(strcmpi(theithRxn.enz.donor,'GDP-beta-L-fucose'))% CMP-Sialic Acid
                                sugarnucl = 'GDP-Fuc';
                                nucl = 'GDP';
                            else
                                error('MATLAB:GNAT:NOTSUPPORTEDSUGAR','SUGAR NAME IS NOT SUPPORTED');
                            end
                            
                            sugarnames  = {sugarnucl,nucl};
                            sugarvalues = [sugarnucldb(sugarnucl),sugarnucldb(nucl)];
                            
                            theithRxn.rxnkinetics.setMathFormula(theithRxn.enz.id,...
                                theithRxn.id,sugarnames,sugarvalues);
                        elseif(length(varargin)==2)
                            sugarnucldb = varargin{1};
                            ifnegrev    = varargin{2};
                            
                            if(strcmpi(theithRxn.enz.donor,'UDP-N-acetyl-D-glucosamine'))
                                sugarnucl = 'UDP_GlcNAc';
                                nucl = 'UDP';
                            elseif(strcmpi(theithRxn.enz.donor,'UDP-alpha-D-galactose'))
                                sugarnucl = 'UDP_GaL';
                                nucl = 'UDP';
                            elseif(strcmpi(theithRxn.enz.donor,'CMP-N-acetylneuraminate'))% CMP-Sialic Acid
                                sugarnucl = 'CMP_NeuAc';
                                nucl = 'CMP';
                            elseif(strcmpi(theithRxn.enz.donor,'GDP-beta-L-fucose'))% CMP-Sialic Acid
                                sugarnucl = 'GDP_Fuc';
                                nucl = 'GDP';
                            else
                                error('MATLAB:GNAT:NOTSUPPORTEDSUGAR','SUGAR NAME IS NOT SUPPORTED');
                            end
                            
                            sugarnames  = {sugarnucl,nucl};
                            sugarvalues = [sugarnucldb(sugarnucl),sugarnucldb(nucl)];
                            
                            theithRxn.rxnkinetics.setMathFormula(theithRxn.enz.id,...
                                theithRxn.id,sugarnames,sugarvalues,ifnegrev);
                        end
                        
                    end
                end
            end
        end
        
        %set transport reaction kinetics for each transport reaction
        function setTransportRxnKinetics(obj)
            % set up transport kinetics
            % set up transport reaction id
            
            for i = 1:obj.theRxns.length
                theithRxn = obj.theRxns.get(i);
                if(isempty(theithRxn.enz))
                    theithRxn.setrxnkinetics;
                    theithRxn.id=['glycantransport',num2str(i)];
                end
            end
            
            %set up formula
            for i = 1 : obj.theRxns.length
                theithRxn = obj.theRxns.get(i);
                if(isempty(theithRxn.enz))
                    try 
                        theithRxn.rxnkinetics.setMathFormula(theithRxn.id);
                    catch error
                      disp('debug');
                    end
                end
            end
        end
        
        function rxnlistbysameEnz = findRxnsBySameEnz(obj,ithEnz)
            rxnlistbysameEnz = [];
            for i = 1 : obj.theRxns.length
                 ithRxnToCompare = obj.listofRxns{i};
                 if(strcmp(ithRxnToCompare.enz.id,ithEnz.id))
                    rxnlistbysameEnz = [rxnlistbysameEnz;i]; 
                 end    
            end
        end        
        
        function theCompetingRxns = findEnzCompetingRxns(obj,ithRxn)  % to be tested
            theCompetingRxns = CellArrayList;
            for i = 1 : obj.theRxns.length
                ithRxnToCompare = obj.theRxns.get(i);
                if(ithRxnToCompare==ithRxn)
                    continue;
                end
                
                if(isempty(ithRxnToCompare.enz))
                    continue;
                end
                
                if(ithRxnToCompare.reac.compartment...
                      ~=ithRxn.reac.compartment)
                   continue;
                end
               
                if(ithRxnToCompare.enz.ecno==ithRxn.enz.ecno)
                    theCompetingRxns.add(ithRxnToCompare);
                end
            end
        end
    end
    
    methods  %getter methods
        function res = getInitGlycans(obj)
            % getInitGlycans get the initial glycan species
            %
            % See also Pathway
            res = obj.initGlycans;
        end
        
        function res = getFinalGlycans(obj)
            % getFinalGlycans get the final glycan species
            %
            % See also Pathway
            res = obj.initGlycans;
        end
        
        function setInitGlycans(obj,initGlycans)
            % setInitGlycans set the initial glycan species
            %
            % See also Pathway
            obj.initGlycans= initGlycans;
            
            if (~obj.theSpecies.contains(initGlycans))
                obj.theSpecies.add(initGlycans)
            end
        end
        
        function setFinalGlycans(obj,finalGlycans)
            % setFinalGlycans set the initial glycan species
            %
            % See also Pathway
            obj.finalGlycans= finalGlycans;
            if (~obj.theSpecies.contains(finalGlycans))
                obj.theSpecies.add(finalGlycans)
            end
        end
        
        function setCompartment(obj,compartment)
            % setCompartment set the compartment in the network
            %
            % See also Pathway
            obj.compartment = compartment;
        end
        
        function res = getCompartment(obj)
            % getCompartment get the compartment in the network
            %
            % See also Pathway
            res = obj.compartment;
        end % end
        
        function res = getReactions(obj)
            % getReactions get the reactions in the network
            %
            % See also Pathway
            res = obj.theRxns;
        end % end
        
        function res = getNReactions(obj)
            % getNReactions get the number of the reactions in the network
            %
            % See also Pathway
            res = obj.theRxns.length;
        end % end
        
        function res = getEnzymes(obj)
            % getEnzymes get the enzymes in the network
            %
            % See also Pathway
            res = obj.enzs;
        end % end
        
        function res = getName(obj)
            % getName get the name of  the network
            %
            % See also Pathway
            res = obj.name;
        end % end
        
        function res = getNSpecies(obj)
            % getNSpecies get the number of the species in  the network
            %
            % See also Pathway
            res = obj.theSpecies.length;
        end % end
        
        function res = getSpecies(obj)
            % getSpecies get the species in  the network
            %
            % See also Pathway
            res = obj.theSpecies;
        end % end
        
        function res = getGlycanStruct(obj,i)
            % getGlycanStruct get the glycan structure of ith species in  the network
            %
            % See also Pathway
            res = obj.theSpecies.get(i).getGlycanStruct;
        end % end
    end
    
    methods  % setter methods
        function setReactions(obj,reactions)
            % SETREACTIONS set the reactions in the network
            %
            %See also PATHWAY
            obj.theRxns = reactions;
        end % end
        
        function setSpecies(obj,species)
            % SETSPECIES set the species in the network
            %
            %See also PATHWAY
            obj.theSpecies = species;
        end % end
        
        function setName(obj,name)
            % SETNAME set the name of the network
            %
            % See also PATHWAY
            obj.name = name;
        end % end
        
        function setEnzymes(obj,enzymes)
            % SETENZYMES set the enzymes in the network
            %
            % See also PATHWAY
            obj.enzs = enzymes;
        end % end
        
        function setSpeciesIsolated(obj,glycanstructobj)
            % setSpeciesIsolated set the species in the network isolated
            %
            % See also PATHWAY
            for i = 1 : obj.getNSpecies
                ithspecies = obj.theSpecies.get(i);
                
                if(ithspecies.glycanStruct.equalStruct(glycanstructobj))
                    ithspecies.listOfReacRxns=CellArrayList();
                    ithspecies.listOfProdRxns=CellArrayList();
                    continue;
                end
                
                removereactindex = [];
                for ii = 1: ithspecies.listOfReacRxns.length
                    theprod = ithspecies.listOfReacRxns.get(ii).getProduct;
                    if(theprod.glycanStruct.equalStruct(glycanstructobj))
                        removereactindex =[removereactindex;ii];
                    end
                end
                if(~isempty(removereactindex))
                    ithspecies.listOfReacRxns.remove(removereactindex);
                end
                
                removeprodtindex = [];
                for ii = 1: ithspecies.listOfProdRxns.length
                    thereact = ithspecies.listOfProdRxns.get(ii).getReactant;
                    if(thereact.glycanStruct.equalStruct(glycanstructobj))
                        removeprodtindex =[removeprodtindex;ii];
                    end
                end
                
                if(~isempty(removeprodtindex))
                    ithspecies.listOfProdRxns.remove(removeprodtindex);
                end
            end
            
            removeindex = [];
            for i = 1 : obj.getNReactions
                ithrxn = obj.theRxns.get(i);
                theReact = ithrxn.getReactant.glycanStruct;
                theProd  = ithrxn.getProduct.glycanStruct;
                if(theReact.equalStruct(glycanstructobj))
                    removeindex = [removeindex;i];
                    continue;
                end
                
                if(theProd.equalStruct(glycanstructobj))
                    removeindex = [removeindex;i];
                    continue;
                end
            end
            obj.theRxns.remove(removeindex);
        end
        
    end
    
    methods
        glycanNetJava = pathwayMat2Java(obj )
        listofglycansinSmallGlyPep = toGlycanListinSmallGlyPep(obj,varargin);
    end
    
    methods
        function indexComptSubstr = findEnzSubstr(obj,enztofind)  % to be tested
            indexComptSubstr = [];
            for i = 1 : obj.theRxns
                theithRxn = obj.theRxns.get(i);
                theithenz = theithRxn.enz;
                if(isequal(theithenz,enztofind))
                    indexComptSubstr = [indexComptSubstr;...
                        obj.theGlycans.findElt(theithRxn.react)];
                end
            end
        end
        
        function indexComptSubstr = findEnzComptSubstr(obj,enztofind,thesubstr)  % to be tested
            indexComptSubstr = findEnzSubstr(obj,enztofind);
            thesubstrindex   = obj.theGlycans.findElt(thesubstr);
            indexComptSubstr(~(indexComptSubstr-thesubstrindex))=[];
        end
        
        function setAllRxnKinetics(obj)
            for i = 1 : length(obj.theRxns)
                obj.theRxns.get(i).setRxnKinetic;
            end
        end
        
        function obj=addInhibitors(obj)  % to be tested
            % apply inhibition rules
            indexComptSubstr = findEnzComptSubstr(thePathway,obj.enz,obj.react);
            theComptSubstr   = thePathway.theGlycans.get(indexComptSubstr);
            inhibitconsts    = CellArrayList;
            for i = 1 : length(theComptSubstr)
                inhibit.speciesnames = theComptSubstr{i}.id;
                inhibitconsts.add(inhibit);
            end
        end
    end
    
    methods
        function new = clone(obj,varargin)
            %CLONE deep copy a Pathway object
            %  clone(RATHWAYobj) creates a new PATHWAY object
            %
            %  See also Pathway
            
            % Instantiate new object of the same class.
            meta = metaclass(obj);
            new = feval(class(obj));
            for i = 1:length(meta.Properties)
                prop = meta.Properties{i};
                if ~(prop.Dependent || prop.Constant) ...
                        && ~(isempty(obj.(prop.Name)) ...
                        && isempty(new.(prop.Name)))
                    if ~isa(obj.(prop.Name),'handle')
                        new.(prop.Name) = obj.(prop.Name);
                    end
                end
            end
            
            new.theRxns     = CellArrayList();
            new.theSpecies  = CellArrayList();
            new.compartment = CellArrayList();
            new.theEnzs     = CellArrayList();
            
            % add enzs to the pathway
            for i=1:length(obj.theEnzs)
                theEnz = obj.theEnzs.get(i).clone;
                new.theEnzs.add(theEnz);
            end
            
            % add species to the pathway
           if(length(varargin)==1)
              isstructclone = varargin{1};
           else
              isstructclone = 1;  
           end
           
           for i=1:length(obj.theSpecies)
               theGlycan = obj.theSpecies.get(i).clone(isstructclone);
               new.theSpecies.add(theGlycan);
           end
            
            % add reactions to the pathway
            for i=1:length(obj.theRxns)
                theRxnToCopy = obj.theRxns.get(i);
                
                % find reactant and add it to a new reaction
                theReactant = theRxnToCopy.getReactant;
                if(~isempty(theReactant))
                    reacLoc = obj.theSpecies.locationsOf(theReactant);
                    reac = new.theSpecies.get(reacLoc);
                else
                    reac = [];
                end
                
                % find product and add it to a new reaction
                theProduct = theRxnToCopy.getProduct;
                if(~isempty(theProduct))
                    prodLoc = obj.theSpecies.locationsOf(theProduct);
                    prod    = new.theSpecies.get(prodLoc);
                else
                    prod = [];
                end
                
                % find eznyme and add it to a new reaction
                theEnz = theRxnToCopy.getEnzyme;
                if(~isempty(theEnz))
                    enzLoc = obj.theEnzs.locationsOf(theEnz);
                    enz    = new.theEnzs.get(enzLoc);
                else
                    enz = [];
                end
                if(~isempty(enz))
                    theRxnClone = Rxn(reac,prod,enz);
                else
                    theRxnClone = Rxn(reac,prod);
                end
                new.addRxn(theRxnClone);
                %new.theRxns.add(theRxnClone);
            end
            
            % add compartments to the pathway
            for i=1:length(obj.compartment)
                new.addCompts(obj.compartment.get(i).clone);
            end
        end
      
    function new = copy(obj,varargin)
            %COPY copy a Pathway object
            %  copy(RATHWAYobj) creates a new PATHWAY object
            %
            %  See also Pathway
            
            % set id for each species
            obj.setSpeciesIDIndex;
            obj.setRxnsIDIndex;
            obj.setEnzIDIndex;
            
            % set Stoichiometry matrix
            if(isempty(obj.SM))
               obj.setSM;            
            end
            
            if(isempty(obj.EM))
               obj.setEM;            
            end
            
            % Instantiate new object of the same class.
            meta = metaclass(obj);
            new  = feval(class(obj));
            for i = 1:length(meta.Properties)
                prop = meta.Properties{i};
                if ~(prop.Dependent || prop.Constant) ...
                        && ~(isempty(obj.(prop.Name)) ...
                        && isempty(new.(prop.Name)))
                    if ~isa(obj.(prop.Name),'handle')
                        new.(prop.Name) = obj.(prop.Name);
                    end
                end
            end
            
            new.theRxns     = CellArrayList();
            new.theSpecies  = CellArrayList();
            new.compartment = CellArrayList();
            new.theEnzs     = CellArrayList();
            
            % add enzs to the pathway
            for i=1:length(obj.theEnzs)
                theEnz = obj.theEnzs.get(i).clone;
                new.theEnzs.add(theEnz);
            end
            
            % add species to the pathway
           if(length(varargin)==1)
              isstructclone = varargin{1};
           else
              isstructclone = 1;  
           end
           
           for i=1:length(obj.theSpecies)
               theGlycan = obj.theSpecies.get(i).clone(isstructclone);
               new.theSpecies.add(theGlycan);
           end
            
            % set up reactions
            thePathwaySM = obj.SM;
            theEnzymeEM  = obj.EM;
            numRxns = length(obj.theRxns);
            for i=1 : numRxns
                ithrow = thePathwaySM(i,:);
                reacLoc  = find(ithrow==-1,1);
                % find reactant and add it to a new reaction
                if(~isempty(reacLoc))
                    reac    = new.theSpecies.get(reacLoc);
                else
                    reac    = [];
                end
                
                % find product and add it to a new reaction
                prodLoc  = find(ithrow==1,1);
                if(~isempty(prodLoc))
                    prod = new.theSpecies.get(prodLoc);
                else
                    prod = [];
                end
                
                % find eznyme and add it to a new reaction
                enzLoc = theEnzymeEM(i);
                if(enzLoc>0)
                    enz = new.theEnzs.get(enzLoc);
                else
                    enz = [];
                end
                if(~isempty(enz))
                    try 
                        theRxnClone = Rxn(reac,prod,enz);
                    catch error
                        disp('error');
                    end
                else
                    theRxnClone = Rxn(reac,prod);
                end
                new.theRxns.add(theRxnClone);
            end
            
            % add compartments to the pathway
            for i=1:length(obj.compartment)
                new.addCompts(obj.compartment.get(i).clone);
            end
        end
    end
end