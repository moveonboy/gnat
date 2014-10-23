function mspeak = synMS(glycanReactorModel,glycanConcDB,varargin)
%SYNMS generate a theoretical MS based on a  
%  glycosylation reaction network
% 
%  mspeak = synMS(glycanReactorModel,glycanConcDB,msoption)
% 
%See also msprocess,readMS.

% Author: Gang Liu
% Date Lastly Updated: 5/21/14

narginchk(2,3);

if(nargin==2)
    peakdifftol = 0.2;
else
    msoption    = varargin{1};
    peakdifftol = msoption.peakdifftol;
end    

if (~isa(glycanReactorModel,'GlycanNetModel'))||...
   (~isa(glycanConcDB,'containers.Map'))     
    error('MATLAB:GNAT:ERRORINPUTTYPE','INCORRECT INPUT TYPE');
end

allglycans = glycanReactorModel.glycanpathway.theSpecies;
nglycans   = allglycans.length;

mspeak = [];
for i = 1 : nglycans
    ithglycan              = allglycans.get(i);
    glycanformula          = ithglycan.glycanStruct.computeFormula;
    glycandistr            = isotopicdist(glycanformula);
    glycanfraction         = glycanConcDB(ithglycan.id);
    glycandistr(:,2)       = glycandistr(:,2)*glycanfraction*100;
    mspeak                 = joinMSpeaks(mspeak,glycandistr,...
                              peakdifftol);
end

% sort data based on peak order
mspeak = sort(mspeak,1);

end

function mspeak = joinMSpeaks(mspeak,relativeglycandistr,peakdifftol)
  if(isempty(mspeak))
      mspeak =  [mspeak;relativeglycandistr];
  else
      for i = 1 : size(relativeglycandistr,1)
        ms      =  relativeglycandistr(i,1);
        mzmatch =  find(abs(mspeak-ms)<peakdifftol);
        if(isempty(mzmatch))
           mspeak = [mspeak;relativeglycandistr(i,:)];
        else 
           if(length(mzmatch)==1)
               
               mspeak(mzmatch,2)= mspeak(mzmatch,2) + relativeglycandistr(i,2); 
           else
               error('MATLAB:GNAT:ERRORGLYCANPEAKINPUT','INCORRECT GLYCAN INPUTS');
           end    
        end
      end 
  end
end

