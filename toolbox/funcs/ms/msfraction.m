function [glycanlist]= msfraction(peaklist,pfwhh,expecGlycan,glycanmwarray,varargin)
%msfraction
% 
% [glycanlist]= msfraction(peaklist,pfwhh,expecGlycan,glycanmwarray) uses 
%   default options which are 0.3 for oversegmentationfilter, 5 for 
%   heightfilter, true for showplot. The peak list p is a matrix with two
%   columns of peak locations, and peak intensities. The matrix pfwhh has
%   two columns indicating the left and right location of the full width 
%   height.
%
% Example:
%     
%See also readMS,readMsprocess.

%Date Lastly Updated: 02/17/14
relativeabundance  = zeros(length(expecGlycan),1);
isotopicnumarray   = zeros(length(expecGlycan),1);
allpeaktotalarea   =0;
monoisomw=[];
abundantmw=[];


for i=1:length(expecGlycan)
    mddouble=glycanmwarray{i,1};
    [peakarea, isotopicenvelope] = findpeakarea(mddouble,peaklist,pfwhh);
    relativeabundance(i,1)       = peakarea;
    allpeaktotalarea             = allpeaktotalarea+peakarea;
    isotopicnumarray(i,1)        = isotopicenvelope;
    monoisomw{i,1}=glycanmwarray{i,1}(1,1);
    maxfraction=max(glycanmwarray{i,1}(:,2));
    numrow=find(~(glycanmwarray{i,1}(:,2)-maxfraction));
    abundantmw{i,1}=glycanmwarray{i,1}(numrow,1);
    glycanresiduestring=expecGlycan{i,1};
    glycan1letstring=gly1charformat(glycanresiduestring);
    options.mono=false;
    avaragemw{i,1}=glycanMolWt(glycan1letstring,options);
end    

for i = 1 : length(expecGlycan)
    relativeabundance(i,1)       = relativeabundance(i,1)/allpeaktotalarea;
end

for i=1:length(expecGlycan)
    glycanlist(i).relativeabundance=relativeabundance(i,1);
    glycanlist(i).glycancompos     = expecGlycan(i);
end

if(length(varargin)==1)
    outputfilename=varargin{1};
    A1=cellstr('Composition');
    B1=cellstr('Fraction');
    C1=cellstr('Monoisotopic mass');
    D1=cellstr('Abundant mass');
    E1=cellstr('Avarage mass');
    xlswrite(outputfilename,relativeabundance,1,'B2');
    xlswrite(outputfilename,expecGlycan,1,'A2');
    xlswrite(outputfilename, monoisomw,1,'C2');
    xlswrite(outputfilename, abundantmw,1,'D2');
    xlswrite(outputfilename, avaragemw,1,'E2');
    xlswrite(outputfilename,A1,1,'A1');
    xlswrite(outputfilename,B1,1,'B1');
    xlswrite(outputfilename,C1,1,'C1');
    xlswrite(outputfilename,D1,1,'D1');
    xlswrite(outputfilename,E1,1,'E1');
end




    
    



