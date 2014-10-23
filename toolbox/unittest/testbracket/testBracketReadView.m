brackettest = glycanMLread('brackettest.glycoct_xml');
glycanjava  = brackettest.structMat2Java;
brackettest.glycanjava = glycanjava;
glycanMLwrite(brackettest,'brackettestglycan.xml');
glycanViewer(brackettest);
