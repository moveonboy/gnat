clc; clear;
brackettest = glycanMLread('brackettestglycan2.glycoct_xml');
glycanjava  = brackettest.structMat2Java;
brackettest.glycanjava = glycanjava;
glycanMLwrite(brackettest,'brackettestglycan.xml');
glycanViewer(brackettest);
