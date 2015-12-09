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
