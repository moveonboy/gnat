classdef ResidueType < hgsetget
    %ResidueType class representing properties for the monosaccharide type
    %
    % ResidueType properties:
    %  name                - name of the residue type
    %  desc                - detailed description of the residue
    %  synonyms            - the synonyms of the residue type
    %  superclass          - residue type superclass
    %  IUPACNAME           - IUPAC name
    %  anomericCarbon      - anomeric carbon position
    %  chirality           - the chirality
    %  ringSize            - the size of the monosaccharide ring structure
    %  isSaccharide        - the logical value if the residue is saccharide
    %  nAcetyls            - number of acetyl groups
    %  nMethyls            - number of methyl groups
    %  nlinkages           - number of linkages
    %  linkagePos          - the positions of linkages
    %  chargesPos          - charge positions
    %  compositionClass    - the composition class
    %  resMassMain         - mass used currently for the residue
    %  resMassAvg          - average mass used currently for the residue
    %  dropMthylated       - if the iron drop methylated
    %  dropAcetylated      - if the iron drop acetylated
    %  isCleavable         - if the structure cleavable
    %  isLabile            - if the structure alterable
    %  canRedend           - if the residue can be reducing end
    %  canParent           - if the residue can be parent
    %
    %
    % ResidueType  methods:
    %  ResidueType         -  create a ResidueType  object
    %  resTypetoString     -  export residue type properties to a string
    %  resTypeMat2Java     -  convert a MATLAB ResidueType object to a Java ResidueType object
    %  clone               -  copy a ResidueType object
    %
    % See also GlycanResidue,GlycanStruct.
    
    % Author: Gang Liu
    % CopyRight 2012 Neelamegham Lab
    
    properties
        %NAME the name of the residue type
        %   NAME is a character array
        %
        % See also ResidueType
        name                             %1
    end
    
    properties
        % DESC description of the residue type
        %   DESC is a character array
        %
        % See also ResidueType
        desc                               %2
    end
    
    properties
        % SYNONYMS residue name synonyms
        %     SYNONYMS is a character array
        %
        % See also ResidueType
        synonyms                     %3
    end
    
    properties
        % SUPERCLASS residue type superclass
        %     superclass is a character array
        %
        % See also ResidueType
        superclass                   %4
    end
    
    properties
        % IUPACNAME name in IUPAC format
        %     IUPACNAME is a character array
        %
        % See also ResidueType
        iupacName                 %5
    end
    
    properties
        % ANOMERICCARBON the position of anomeric carbon
        %    ANOMERICCARBON is a character array
        %
        % See also ResidueType
        anomericCarbon        %6\
    end
    
    properties
        % CHIRALITY the stereisomer property
        %   CHIRALITY  is a char
        %
        % See also ResidueType
        chirality                          %7
    end
    
    properties
        %  RINGSIZE  the size of the ring structure
        %     RINGSIZE is a numeric type (double)
        %
        % See also ResidueType
        ringSize                         %8
    end
    
    properties
        %  ISSACCHARIDE describes if the residue is a saccharide or not
        %    ISSACCHARIDE is a logic value
        %
        % See also ResidueType
        isSaccharide               %9
    end
    
    properties
        %  NACETYLS the number of the acetyles
        %     NACETYLS is a numeric type (double)
        %
        %
        % See also ResidueType
        nAcetyls                        %10
    end
    
    properties
        %  NMETHYLS the number of methyles
        %     NMETHYLS is a numeric type (double)
        %
        %
        % See also ResidueType
        nMethyls                       %11
    end
    
    properties
        %  NLINKAGES the number of linkages
        %      NLINKAGES is a numeric type
        %
        % See also ResidueType
        nlinkages                     %12
    end
    
    properties
        %  LINKAGEPOS the position of the linkage
        %     LINKAGEPOS is a character array
        %
        % See also ResidueType
        linkagePos                  %13
    end
    
    properties
        % CHARGEPOS the position of the charge in the structure
        %     CHARGEPOS is a numeric scalar (double)
        %
        % See also ResidueType
        chargesPos                %14
    end
    
    properties
        % COMPOSITIONCLASS the class of the composition
        %     COMPOSITIONCLASS  is a character array
        %
        % See also ResidueType
        compositionClass;   %15
    end
    
    properties
        %  RESMASSMAIN monoisotopic mass of the residue
        %      RESMASSMAIN is a numeric scalar (double)
        %
        % See also ResidueType
        resMassMain           %16
    end
    
    properties
        %  RESMASSAVG average  mass of the residue
        %      RESMASSAVG is a numeric scalar (double)
        %
        % See also ResidueType
        resMassAvg             %17
    end
    
    properties
        %   DROPMTHYLATED
        %       DROPMTHYLATED is a logic scalar
        %
        % See also ResidueType
        dropMthylated          %18
    end
    properties
        %   DROPACETYLATED
        %       DROPACETYLATED is a logic scalar
        %
        % See also ResidueType
        dropAcetylated         %19
    end
    
    properties
        %  ISCLEAVABLE indicate if the residue is cleavable
        %      ISCLEAVABLE is a logic scalar
        %
        % See also ResidueType
        isCleavable             %20
    end
    
    properties
        %  ISLABILE indicate if the residue is easily altered
        %      ISLABILE is a logic scalar
        %
        % See also ResidueType
        isLabile                    %21
    end
    
    properties
        % CANREDEND indicate if the residue can be reducing end
        %     CANREDEND is a logic scalar
        %
        % See also ResidueType
        canRedend             %22
    end
    
    properties
        % CANPARENT indicate if the residue can be parent residue
        %     CANPARENT is a logic scalar
        %
        %
        % See also ResidueType
        canParent               %23
    end
    
    methods  % constructor
        function obj= ResidueType(varargin)
            %ResidueType create a ResidueType object
            %
            %      RT = ResidueType() creates a default ResidueType.
            %
            %      RT = ResidueType(RTJAVA) converts a Java ResidueType object
            %      to a MATLLAB ResidueType object.
            %
            %  See also GlycanResidue,GlycanStruct.
            
            % validate the number of input arguments
            error(nargchk(0,1,nargin));
                       
            if(nargin==0)
                return
%                 obj.name = '#empty';   %1
%                 obj.desc = 'Empty'; %2
%                 obj.synonyms ='Empty'; %3
%                 obj.superclass = 'special'; %4
%                 obj.iupacName = '-'; %5
%                 obj.anomericCarbon = '?'; %6
%                 obj.chirality = '?'; %7
%                 obj.ringSize = '?'; %8
%                 
%                 obj.isSaccharide = 0; %9
%                 obj.nAcetyls = 0; %10
%                 obj.nMethyls = 0; %11
%                 obj.nlinkages = 0;%12
%                 
%                 obj.linkagePos ='?'; %13
%                 obj.chargesPos ='?'; %14
%                 obj.compositionClass = '?'; %15
%                 
%                 obj.resMassMain = 0.0;%17
%                 obj.resMassAvg =0.0 ;%18
%                 obj.dropMthylated = 0;%19
%                 obj.dropAcetylated = 0;%20
%                 
%                 obj.isCleavable = 0;%21
%                 obj.isLabile = 0; %22
%                 obj.canRedend = 0;%23
%                 obj.canParent = 0;%24
            elseif(nargin==1)
                firstArg = varargin{1};
                import org.eurocarbdb.application.glycanbuilder.ResidueType;
                if(isa(firstArg,'org.eurocarbdb.application.glycanbuilder.ResidueType'))
                    obj.name = char(firstArg.getName); %1
                    obj.superclass =  char(firstArg.getSuperclass); %4
                    obj.compositionClass = char(firstArg.getCompositionClass);  %15
                    obj.synonyms = char(firstArg.getSynonyms); %3
                    obj.iupacName = char(firstArg.getIupacName); %5
                    obj.anomericCarbon = char(firstArg.getAnomericCarbon);%6
                    obj.chirality = char(firstArg.getChirality); %7
                    obj.ringSize =char(firstArg.getRingSize); %8
                    obj.isSaccharide = firstArg.isSaccharide; %9
                    obj.isCleavable = firstArg.isCleavable; %21
                    obj.isLabile = firstArg.isLabile; %22
                    obj.resMassMain =   firstArg.getResidueMassMain; %17
                    obj.resMassAvg   =  firstArg.getResidueMassAvg; %18
                    obj.nMethyls = firstArg.getNoMethyls; %11
                    obj.dropMthylated = firstArg.isDroppedWithMethylation;%19
                    obj.nAcetyls = firstArg.getNoAcetyls; %10
                    obj.dropAcetylated = firstArg.isDroppedWithAcetylation; %20
                    obj.nlinkages       = firstArg.getMaxLinkages;      %12
                    obj.linkagePos    = firstArg.getLinkagePositions; %13
                    obj.chargesPos   = firstArg.getChargePositions; %14
                    obj.canRedend = firstArg.canBeReducingEnd;%23
                    obj.canParent = firstArg.canHaveParent;   %24
                    obj.desc = char(firstArg.getDescription); %2
                end
            end
        end
    end
    
    methods
        resTypeJava = resTypeMat2Java(obj)
    end
    
    methods
        function toString = resTypetoString(obj)
            % resTypetoString export residue type properties to a string
            %
            %     RTSTRING = resTypeToString(obj) output residue type
            %     properties to a string.
            %
            %
            % See also RESIDUETYPE
            
        toString = sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%i\t%i\t%i\t%i\t%s\t%i\t%s\t%i\t%s\t%s\t%s\t%s\t%s',...
              obj.name,isEmptyString(obj.superclass),obj.compositionClass,...
              isEmptyString(obj.synonyms),obj.iupacName,obj.anomericCarbon,...];%6
              obj.chirality,obj.ringSize,parseBool2Str(obj.isSaccharide), ...];%9
              parseBool2Str(obj.isCleavable),parseBool2Str(obj.isLabile), 0,...
              obj.resMassMain,obj.resMassAvg,obj.nMethyls,...
             parseBool2Str(obj.dropMthylated),obj.nAcetyls,...
             parseBool2Str(obj.dropMthylated),obj.nlinkages,...
             isEmptyString(obj.linkagePos),isEmptyString(obj.chargesPos),...
             parseBool2Str(obj.canRedend),parseBool2Str(obj.canParent),...
             obj.desc);
                    
            function boolStr = parseBool2Str(varBool)
                if(varBool)
                    boolStr='true';
                else
                    boolStr='false';
                end
            end
            
            function results = isEmptyString(str)
                if(isempty(str))
                    results='?';
                else
                    if(size(str,1)==1)
                        results=str(1);
                    else
                        results=str';
                        results=results(1);
                    end
                end
            end
        end
    end
    
    methods
        function new = clone(obj)
            %CLONE copy a ResidueType object
            %  CLONE(RESIDUETYPEobj) creates a new ResidueType object
            %
            %  See also ResidueType
            
            new = ResidueType;
            new.name     = obj.name;   %1
            new.desc     = obj.desc; %2
            new.synonyms =    obj.synonyms; %3
            new.superclass =  obj.superclass; %4
            new.iupacName  =  obj.iupacName ; %5
            new.anomericCarbon =  obj.anomericCarbon; %6
            new.chirality = obj.chirality; %7
            new.ringSize =     obj.ringSize ; %8
                
             new.isSaccharide =    obj.isSaccharide; %9
             new.nAcetyls =   obj.nAcetyls; %10
             new.nMethyls =    obj.nMethyls; %11
             new.nlinkages =   obj.nlinkages;%12
                
             new.linkagePos =   obj.linkagePos; %13
             new.chargesPos =   obj.chargesPos; %14
             new.compositionClass =    obj.compositionClass; %15
                
             new.resMassMain =    obj.resMassMain;%17
             new.resMassAvg =   obj.resMassAvg;%18
             new.dropMthylated =   obj.dropMthylated;%19
             new.dropAcetylated =   obj.dropAcetylated; %20
                
             new.isCleavable =   obj.isCleavable;%21
             new.isLabile =  obj.isLabile; %22
             new.canRedend =  obj.canRedend;%23
             new.canParent =  obj.canParent;%24
            
            
    %  name                - name of the residue type
    %  desc                - detailed description of the residue
    %  synonyms            - the synonyms of the residue type
    %  superclass          - residue type superclass
    %  IUPACNAME           - IUPAC name
    %  anomericCarbon      - anomeric carbon position
    %  chirality           - the chirality
    %  ringSize            - the size of the monosaccharide ring structure
    %  isSaccharide        - the logical value if the residue is saccharide
    %  nAcetyls            - number of acetyl groups
    %  nMethyls            - number of methyl groups
    %  nlinkages           - number of linkages
    %  linkagePos          - the positions of linkages
    %  chargesPos          - charge positions
    %  compositionClass    - the composition class
    %  resMassMain         - mass used currently for the residue
    %  resMassAvg          - average mass used currently for the residue
    %  dropMthylated       - if the iron drop methylated
    %  dropAcetylated      - if the iron drop acetylated
    %  isCleavable         - if the structure cleavable
    %  isLabile            - if the structure alterable
    %  canRedend           - if the residue can be reducing end
    %  canParent           - if the residue can be parent
            
            % Instantiate new object of the same class.
%             new = feval(class(obj));
%             
%             % Copy all non-hidden properties.
%             p = properties(obj);
%             for i = 1:length(p)
%                 if(isa(obj.(p{i}),'handle'))
%                     new.(p{i}) = obj.(p{i}).clone;
%                 else
%                     new.(p{i}) = obj.(p{i});
%                 end
%             end
        end
    end
    
    
end

