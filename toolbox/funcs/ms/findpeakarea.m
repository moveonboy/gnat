function  [peaktotalarea,isotopicenvelope]= findpeakarea(md,peaklist,pfwhh)
mz   = md(:,1);
dist = md(:,2);

isotopicmz = mz(1,1);
isotopicdist = dist(1,1);

isotopicindex = findclosetpeak(isotopicmz,peaklist);

expisotopicpeak = peaklist(isotopicindex,1);
%disp(expisotopicpeak);

if(isempty(isotopicindex))  % if the isotopic mass is not found, use most aboundant mass
    mostabtmz   = mz(2,1);  % default: the second peak as default
    isotopicindex = findclosetpeak(mostabtmz,peaklist);
    startithpos = 3;
else
    startithpos = 2;
end

if(isempty(isotopicindex))
   peaktotalarea=0;
   isotopicenvelope=0;
   return 
end


isotopicindexarea = abs(pfwhh(isotopicindex,1)-pfwhh(isotopicindex,2))...
    *peaklist(isotopicindex,2);

if(length(isotopicindex)~=1)
    error('ms data not right');
end

peaktotalarea =isotopicindexarea;
isotopicenvelope = 0;

selectedpeak = [isotopicindex];


for i = startithpos : size(mz,1)
    otherpeakmz        = mz(i,1);
    otherdist          = dist(i,1);
    
    checkpeakmz        = abs(peaklist(:,1)-otherpeakmz)<0.5;
    peakindex          = find(checkpeakmz);
    
    if(peakindex==isotopicindex)
        continue;
    end
    
    if(length(find(~(peakindex-selectedpeak)))>1)
        continue;
    end
    
    
    if(isempty(peakindex))
        continue;
    end
    
    numpeaks = length(peakindex);
    peakindexarea=0;
    
    for j=1:numpeaks
        singlepeak = peakindex(j,1);
        peakindexarea =peakindexarea+ abs(pfwhh(singlepeak,1)-pfwhh(singlepeak,2))...
            *peaklist(singlepeak,2);
    end
    
    theorratio=otherdist/isotopicdist;
    peakratio = peakindexarea/isotopicindexarea;
    
%     if((abs(peakratio-theorratio)/theorratio*100)<30)
        peaktotalarea =  peaktotalarea + peakindexarea;
        isotopicenvelope=isotopicenvelope+1;
%     end   
end

end

function isotopicindex = findclosetpeak(isotopicmz,peaklist)
checkisotopicmz = abs(peaklist(:,1)-isotopicmz)<0.99;
isotopicindex   = find(checkisotopicmz);

 % if multiple peaks are found close to isotopic mass, select the closet peak as the match
if(length(isotopicindex)>1) 
        isotopicmzsdif         = abs(peaklist(checkisotopicmz)-isotopicmz);
        isotopicmultipleindex  = find(checkisotopicmz);
        isotopicnewindex       = find(~abs(isotopicmzsdif-min(isotopicmzsdif)));
        isotopicindex          = isotopicmultipleindex(isotopicnewindex);
end

end