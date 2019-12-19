function [p,bc,fzero_outs,eivec]=flfs_bc_findp(bc_mix,flfs,configs,factors)
% [p,bc,fzero_outs,eivec]=flfs_bc_findp(bc_mix,flfs,configs,factors) finds the
% weighting for a falling-sticks model to achieve a desired mixture of horizontal and
% vertical pairwise correlations.
%
% bc_mix:  b/(b+c), where b and c are the weights of the horiz and vertical pairwise correlations
%   that remain in the falling-sticks texture generated from horizontal sticks and vertical
%   sticks of unit correlation
%   If length(bc_mix) is 2, then it is interpreted as b/(b+c)
% flfs: setup for the falling-leaves/falling-sticks model, must be 2x2 standard at present
%    computed if not supplied (faster if supplied)
% configs: result of flfs_enumerate(flfs); computed if not supplied (faster if supplied)
% factors: result of flfs_bc_factors(configs); computed if not supplied (faster if supplied)
%
% p: in [0,1].  p is weight assigned to horizontal sticks.  1-p assigned to vertical sticks
% bc: [b,c] where b is horizontal correlation in the falling-sticks construction and
%    c is the vertical correlation
% fzero_outs.[fval,exitflag,output]:  corresponding outputs of fzero
% eivec: eigenvector of the flfs model for this mixture
%
%  See also: FLFS_DEFINE, FLFS_ENUMERATE, FLFS_FINDEIVEC, FLFS_BC_FACTORS, FLFS_BC_FINDP_OF.
%
exitflag=[];
fval=[];
output=[];
if length(bc_mix)==2
    bc_mix=bc_mix(1)/sum(bc_mix);
end
if (bc_mix==0) | (bc_mix==1)
    p=bc_mix;
    bc=[p 1-p];
else
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
    %
    f=@(p)flfs_bc_findp_of(p,flfs,configs,factors)-bc_mix;
    [p,fval,exitflag,output]=fzero(f,0.5);
end
[bc_mix_check,b,c,eivec]=flfs_bc_findp_of(p,flfs,configs,factors);
bc=[b c];
fzero_outs.fval=fval;
fzero_outs.exitflag=exitflag;
fzero_outs.output=output;
return
    