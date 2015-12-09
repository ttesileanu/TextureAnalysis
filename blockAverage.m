function [out] = blockAverage(img,n,type)

if nargin<3
    type = 'avg';
end

if (strcmpi(type,'sub'))
    out = img(1:2^(n-1):end,1:2^(n-1):end);
    return;
end

C=ones(n)/(n^2);
outN=floor(size(img)/n);
R=conv2(img,C,'valid');
out=R( 1 + n*(0:(outN(1)-1)), 1 + n*(0:(outN(2)-1)) );

% 
% img = double(img);
% 
% %N = size(img);
% outN = floor(size(img)/n);
% out = zeros(outN);
% img = img(1:n*outN(1),1:n*outN(2));
% [ixX,ixY] = meshgrid(1:n:n*outN(1),1:n:n*outN(2));
% for i=0:(n-1)
%     for j=0:(n-1)
%         ind = sub2ind(size(img),ixY+j,ixX+i);
%         out = out + img(ind);
%     end
% end
% out = out / n^2;

if strcmp(type,'wta')
    out = out >= 0.5;
end

% 
% block = zeros(n);
% %count = zeros(outN^2,1);
% index = 1;
% for i = 1:outN
%     for j = 1:outN
%         x = (i-1)*n + 1;
%         y = (j-1)*n + 1;
%         block = img(x:x+n-1,y:y+n-1);
%         out(i,j) = mean(mean(block));            % Mean
%   %      count(index) = sum(sum(block));
%      %   out(i,j) = sum(sum(block)) >= floor(n^2 / 2); % Winner takes all
%         index = index+1;
%     end
% end
