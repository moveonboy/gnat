function outputGlycanNetODEModel(GlycannNetSBMLfilename,outputODEFileName)
%outputGlycanNetODEModel write an ODE M File for glycosylation reaction
% network
%
%  
%
%See also glycanDSim, glycanSSim.

% Author: Gang Liu
% Date Lastly Updated: 5/18/14
glycanNetModel = glycanNetSBMLread(GlycannNetSBMLfilename);

if(glycanNetModel.glycanpathway.theCompts.length>1)
   glycanNetModel.glycanNet_sbmlmodel = convertToSingleCompt(glycanNetModel.glycanNet_sbmlmodel);
end

WriteODEFunction(CHOModel.glycanNet_sbmlmodel,outputODEFileName);

end
