function[sumbiconc,sumtriconc] = getBiTriConcFromBaileyModel(data,biindex, triindex)
    sumspeciesintgn = sum(data(end,100:132));  

    x             = data(:,cell2mat(biindex));
    conc_bi.conc  = x(size(x,1),:);
    sumbiconc     = sum(conc_bi.conc)/sumspeciesintgn;

    x              = data(:,cell2mat(triindex));  
    conc_tri.conc  = x(size(x,1),:);
    sumtriconc     = sum(conc_tri.conc)/sumspeciesintgn;
end
