enzdbmatfilename   = 'glyenzDB.mat';
 enzdb = enzdbmatLoad(enzdbmatfilename);
 mgat2 = enzdb('mgat2');
 mgat3 = enzdb('mgat3');
m3gnspecies    = GlycanSpecies(glycanMLread('m3gn.glycoct_xml')) ;
fprintf(1,'display initial substrate structure\n');
glycanViewer(m3gnspecies.glycanStruct);
substrateArray = CellArrayList;
enzArray = CellArrayList;
substrateArray.add(m3gnspecies);
enzArray.add(mgat2);
enzArray.add(mgat3);
[isPath,nglycanpath]=inferGlyForwPath(substrateArray, enzArray);
if(isPath)
    glycanPathViewer(nglycanpath);
end
     
             
             
          