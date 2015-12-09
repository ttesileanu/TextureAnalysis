function [img_out] = binarize(img,gamma)
% Binarize img so that the fraction of 'on' pixels is (1+gamma)/2

if (gamma < -1) || (gamma > 1)
    display('gamma must lie between -1 and 1');
    return;
end

p = (1+gamma)/2;
threshold = quantile(img(:),p);


img_out = double(img >= threshold);

% [N,M] = size(img);
% 
% vect = reshape(img,N*M,1);
% vect = sort(vect);
% 
% if p==1
%     index = N*M;
% else 
%     index = 1 + floor(p*N*M);
% end
% 
% threshold = vect(index);
% clear vect;

