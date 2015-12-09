function [dataNI, Ffilter] = processData_BSR()

imgNamesFile='./BSR_Test_Index.txt';
imgDirectory='./BSRLumFiles/';
imgDirectory_filter = './NaturalImages/';
imgNamesFile_filter='./Natural_Images_Large_NoSky_Index.txt';
%imgNamesFile_filter='./Natural_Images_TestFilter_Index.txt';
analysesDirectory='./BSRAnalyses/'; %directory where statistics file will be stored
analysesFileName='Segmented_Images'; %name of output statistics file

% which image to be considered most in-focus
focusImg = 25;

%collect statistics (only need to do this once)
generateStatFiles

%construct data structure from all analyses & run 2-Gaussian decomposition
loadImageAnalyses
runFocusScript
selectFocusGaussian

end
