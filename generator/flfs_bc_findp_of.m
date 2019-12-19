function [bc_mix,b,c,eivec]=flfs_bc_findp_of(p,flfs,configs,factors)
% [bc_mix,b,c,eivec]=flfs_bc_findp_of(p,flfs,configs,factors):  objective function for flfs_bc_findp
%
% this finds the mixture of falling sticks that  yields a given ratio of horizontal
% and vertical correlations
% flfs_bc_findp_of
%
% p: in [0,1].  p is weight assigned to horizontal sticks.  1-p assigned to vertical sticks
% flfs: setup for the falling-leaves/falling-sticks model, must be 2x2 standard at present
%    computed if not supplied (faster if supplied)
% configs: result of flfs_enumerate(flfs); computed if not supplied (faster if supplied)
% factors: result of flfs_bc_factors(configs); computed if not supplied (faster if supplied)
%
% bc_mix:  b/(b+c), where b and c are the weights of the horiz and vertical pairwise correlations
%   that remain in the falling-sticks texture generated from horizontal sticks and vertical
%   sticks of unit correlation
% b, c: see above
% eivec: eigenvector of the flfs model
%
%  See also: FLFS_DEFINE, FLFS_ENUMERATE, FLFS_FINDEIVEC, FLFS_BC_FACTORS, FLFS_BC_FINDP.
%
if (nargin<=1)
    flfs=[];
end
if (nargin<=2)
    configs=[];
end
if (nargin<=3)
    factors=[];
end
if isempty(flfs)
    flfs=flfs_define([],mtc_define(2));
end
if isempty(configs)   
    configs=flfs_enumerate(flfs);
end
if isempty(factors)
    factors=flfs_bc_factors(configs);
end
ncats=length(flfs.cats);
nsubcats=length(flfs.subcats);
p_mix=zeros(1,nsubcats);
p_mix(flfs.cats{1}.subcats)=p;
p_mix(flfs.cats{2}.subcats)=1-p;
eivec=flfs_findeivec(configs,p_mix,flfs);
b=factors.b'*eivec;
c=factors.c'*eivec;
bc_mix=b/(b+c);
return
    