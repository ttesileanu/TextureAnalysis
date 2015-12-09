function [ P ev entropy ] = getStats( imgs )
%getStats Take in an ensemble of images, returns the statistics associated
%with them
%   Detailed explanation goes here
P=zeros(size(imgs,3),16);
ev=zeros(size(imgs,3),10);
entropy=zeros(size(imgs,3),1);
for i=1:size(imgs,3)
    [P(i,:) ev(i,:) entropy(i)]=processBlock(imgs(:,:,i));
end

end

