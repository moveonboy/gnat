clc; clear;
glycanwobracket = glycanMLread('glycanwobracket.glycoct_xml');
glycanwobracket.addBracket;
residueMap = load('residueTypes.mat');
galresidue = residueMap.allresidues('Gal');
glycanwobracket.addResidueToBracket(galresidue);
glycanViewer(glycanwobracket);
