% selectFocusGaussian Select the Gaussian from the Gaussian mixture that is
% most likely to correspond to in-focus images.
%   This script uses the hint given by the variable 'focusImg' (which
%   should be predefined when the script is called). This variable
%   identifies the index of an image that is judged to be mostly in-focus.
%   The Gaussian mixture corresponding to in-focus image is chosen as the
%   one to which most of the focus image's patch correspond.
%
%   A member 'focus' is added to each of the entries in the dataNI.indA
%   structure array. It is a structure with members
%    chist:
%       A vector describing, for each patch of the focus image, which of
%       the Gaussian mixture components it belongs to.
%    component:
%       The Gaussian mixture component that is judged to be in-focus. This
%       is the one that has the majority of entries in the chist vector.
%    ev:
%       Subset of the 'ev' matrix corresponding to the in-focus patches.
%       (see analyzeImageSetModNoPC)
%    covM:
%       Covariance matrix for the ev matrix restricted to the in-focus
%       patches.

% image 'cd03A/DSC_0049.JPG' is an in-focus image that we can use to tag the
% "in-focus" gaussian group; this is image 10 in 'Natural_Images_Test_Index.txt'
testImg = focusImg; 

% select images w/in "in-focus" Gaussian group

for i=1:length(dataNI.indA)
    cs = dataNI.indA(i).cx(dataNI.indA(i).ic.image == testImg);
    dataNI.indA(i).focus.chist = cs;
    c = round(dataNI.indA(i).mn(testImg));
    
    dataNI.indA(i).focus.component = c;
    dataNI.indA(i).focus.ev = dataNI.indA(i).ev(dataNI.indA(i).cx==c,:);
    dataNI.indA(i).focus.covM = squeeze(dataNI.indA(i).obj.Sigma(:,:,c));
end

clear c testImg i cs