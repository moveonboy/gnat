% test bracket species
clc;
clear;
storedvarpath = load('Tri_Nglycanpathway.mat');
binglycanpath = storedvarpath.nlinkedpath;
binglycanpath.resetjava;
% glycanPathViewer(binglycanpath);
fprintf(1,'network size: %d Species *%d Reactions\n',...
    binglycanpath.getNSpecies,binglycanpath.getNReactions);
newpathway = pathreduction(binglycanpath);
glycanPathViewer(newpathway);
fprintf(1,'network size: %d Species *%d Reactions\n',...
    newpathway.getNSpecies,newpathway.getNReactions);