function glycanPathViewer(path)
%glycanPathViewer read a Pathway object and return a graphical representation  
% of the glycosylation network.
%
% glycanPathViewer(pathwayObj) visualizes the glycosylation pathway.
%
% Example 1:
%  oglypath=Pathway.loadmat('glypathexample.mat');
%  glycanPathViewer(oglypath);  
%
%
% See also glycanViewer,glycanFileViewer,glycanNetFileViewer,displayset,displayget.

if(isempty(path.compartment))
  golgi  = Compt('golgi');
  name = 'Glycosylation network';
  % set same compartment for each species
  for i = 1 : length(path.theSpecies) 
    path.theSpecies.get(i).compartment = golgi;      
  end
  theCompts = CellArrayList;
  theCompts.add(golgi);
  path.compartment = theCompts; 
  gtestModel = GlycanNetModel(theCompts,path,name);
  glycanNetViewer(gtestModel);
  
  for i = 1 : length(path.theSpecies) 
    path.theSpecies.get(i).compartment = [];      
  end
  path.compartment = CellArrayList;
  
else
  theCompts = path.compartment;
  name = 'Glycosylation network';
  gtestModel = GlycanNetModel(theCompts,path,name);
  glycanNetViewer(gtestModel);
end

end

