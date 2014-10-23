%% Simulate Bailey Model

%% Read Bailey Model from SBML 
% We start this demo by reading the relevant SBML file using glycanNetSBMLread.
% UB1997Model is the resulting GlycanNetModel object. This command will
% take a few moments. Please wait.
clc;
baileyModelSBMLfileName = 'ub1997v11.xml';
UB1997Model = glycanNetSBMLread(baileyModelSBMLfileName);

%% Simulation 
% Below, we show how to simulate the glycosylation reaction network using glycanDSim  
% Note that the current version of glycanDSim uses ode15s to simulate
% the reaction network. This function does not support multicompartment model yet.
% Thus, we write a custom function to convert the model to a single
% compartment including the transport events. 
UB1997Model.glycanNet_sbmlmodel = convertToSingleCompt(UB1997Model.glycanNet_sbmlmodel);
% ubODEFileName = regexprep(UB1997Model.glycanNet_sbmlmodel.name,'.m','');  
ubODEFileName = 'ub1997modelv11';
if(exist(ubODEFileName,'file')~=2)
   WriteODEFunction(UB1997Model.glycanNet_sbmlmodel,ubODEFileName);
end

% set up inputs
endTime = 100;
tsteps  = 100;
options = odeset('RelTol',1e-12,'AbsTol',1e-6);

%run a simulation test based on the parameter value defined in SBML 
rehash toolboxcache;
[tspan, speciesConc, speciesnames] = glycanDSim(UB1997Model,ubODEFileName,endTime, tsteps,options,1);

binames={'M3Gn2species_tgn';'M3Gn2Gspecies_tgn';'M3Gn2Gnbspecies_tgn';'M3Gn2GnbGspecies_tgn'}; %,'M3GnGnb_TGN','M3GnGnbG_TGN'};
biindex   =  findbynameByBaileyModel(speciesnames,binames);
    
triprimenames={'M3Gn3bspecies_tgn';'M3Gn3bGspecies_tgn';'M3Gn3bGnbspecies_tgn';'M3Gn3bGnbGspecies_tgn'}; 
triindex  =  findbynameByBaileyModel(speciesnames,triprimenames);

[sumbiconc,sumtriconc] = getBiTriConcFromBaileyModel(speciesConc,biindex,triindex);
fprintf(1,' the mole fraction of biantennary glycan structures is %f\n',sumbiconc);
fprintf(1,' the mole fraction of tri-prime-anntenary glycan structures is %f\n',sumtriconc);
 
 %% Display Results From Published Paper
 %  Lastly, we run a series of dynamic simulations with different values of
 %  parameter q (productivity rate). The simulation results are the same 
 %  as Fig. 4 in Umana et al 1997. These examples show how simple MATLAB
 %  scripts can be used to augment the functionality available in GNAT!
 qpname = 'qin';
 qpArray=[0:100:2000]/24;
 qpArray(1,1)=1/24;
 plotBi = zeros(1,length(qpArray));
 plotTri =zeros(1,length(qpArray)); 

for i=1:length(qpArray);
     qp = qpArray(i);
     % adjust m file
     adjustBaileyMFile(ubODEFileName, qpname, qp); 
     options = odeset('RelTol',1e-12,'AbsTol',1e-6);
     [tspan, speciesConc, speciesnames] = glycanDSim(UB1997Model,ubODEFileName,endTime, tsteps,options);
     [sumbiconc,sumtriconc]             = getBiTriConcFromBaileyModel(speciesConc,biindex, triindex);
     plotBi(i)                          = sumbiconc;
     plotTri(i)                         = sumtriconc;
     rehash toolboxcache;
end

adjustBaileyMFile(ubODEFileName, qpname, 1/24); 

figure
newqparray = qpArray*24;
plot(newqparray,plotBi,'r',newqparray,plotTri,'b','Linewidth',2);
axis([0 2000 0 1]);
set(gca,'YTick',0:0.25:1);
set(gca,'YTickLabel',{'0','0.25','0.5','0.75','1'});
xlabel('Glycoprotein productivity','fontweight','b','fontsize',12);
ylabel('Mole fraction','fontweight','b','fontsize',12);
legend('Biantennary','Triantennary'); 

%% Reference
%  1. Umana, P. and J.E. Bailey, A mathematical model of N-linked 
%  glycoform biosynthesis. Biotechnol Bioeng, 1997. 55(6): p. 890-908.