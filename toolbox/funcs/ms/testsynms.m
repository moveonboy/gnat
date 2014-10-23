

glymodel = GlycanNetModel.loadmat('glycanNetReactorModel.mat');
glydb = load('specconcdb.mat');
msoption.peakdifftol = 0.2;
mspeak = synMS(glycanReactorModel,glycanConcDB,msoption);
