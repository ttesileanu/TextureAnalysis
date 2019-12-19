function [cgs,mtcs_used]=mtc_probs2cgs(probs,mtcs)
% [cgs,mtcs_used]=mtc_probs2cgs(probs,mtcs) maps a set of  general block probabilities into 
% values for each coordinate group
%
% probs: an array of dimension equal to the number of checks, length on each dimension
%    is number of gray levels
% mtcs:  a structure, typically created by MTC_DEFINE.  If absent, this is
%    created with a number of gray levels equal to size(probs,1) and 2x2 region.  But creation is
%    time-consuming so best to supply
%  mtcs.ng is number of gray levels
%  mtcs.nchecks is number of checks
%
% cgs: [mtcs.coord_group_count,mtcs.ng]: the total probability for 
%    each value of each linear combination, based on the name for the coordinate group
%    every entry of sum(cgs(:,2)) should be equal to sum(probs(:)), which
%    should in turn be 1
% mtcs_used: mtcs, or, if mtcs not supplied,
%    the structure created by mtc_define for 2x2 block probabilities and the number of gray levels in probs
%
% ***correspondence with getcorrs_p2x2 is as follows**
% prand=rand(2,2,2,2); prand=prand/sum(prand(:)); [cgvr,mtcs_used]=mtc_probs2cgs(prand);
% (cgvr(:,1)-cgvr(:,2))'
%    0.0597    0.1362   -0.1951    0.1097   -0.0944    0.0907    0.0902    0.0827    0.0494   -0.1128
% mtcs_used.coord_names
%   'A_1'    'AB_1_1'    'AC_1_1'    'BC_1_1'    'AD_1_1'    'ABC_1_1_1'    'ABD_1_1_1'    'ACD_1_1_1'    'BCD_1_1_1'    'ABCD_1_1_1_1'
% corrs=getcorrs_p2x2(prand); corrs 
%           alpha: -0.1128
%            beta: [0.1362 -0.1951 -0.0944 0.1097]
%           theta: [-0.0907 -0.0902 -0.0827 -0.0494]
%           gamma: -0.0597
% btc_corrs2vec(corrs)
%   -0.0597    0.1362   -0.1951   -0.0944    0.1097   -0.0494   -0.0827   -0.0907   -0.0902   -0.1128
% btc_vec2letcode(btc_corrs2vec(corrs))
%     g: -0.0597 b: 0.1362 c: -0.1951 d: -0.0944 e: 0.1097 t: -0.0494 u: -0.0827 v: -0.0907 w: -0.0902 a: -0.1128
% note sign flips of odd-order coordinates and that the ordering differs
%
% See also:  MTC_DEFINE, MTC_CGS2PROBS, MTC_MAKE_COORD_TABLE, TRICO_CVT, GTC_DEFINE, GETCORRS_P2X2, MTC_CGS2BGC_XFORM.
%
if (nargin<=1)
    mtcs=[];
end
if isempty(mtcs)
    mtcs=mtc_define(size(probs,1));
end
mtcs_used=mtcs;
ng=mtcs.ng;
nchecks=mtcs.nchecks;
nlc=ng^nchecks;
coord_groups=mtcs.coord_groups;
cgs=zeros(mtcs.coord_groups_count,ng);
probs_reshaped=probs(:);
%
coord_table=mtc_make_coord_table(ng,nchecks);
%
lincomb_count=0;
for icg=1:mtcs.coord_groups_count
    %icg
    %coord_groups{icg}.name
    coord_num_primary=coord_groups{icg}.coord_num(1);
    %
    %look across all linear combinations that lead to this coordinate
    %(these differ only in the translated position of the involved checks)
    %
    token_weights=zeros(1,ng);
    lincomb_ptrs=find(mtcs.lincomb_to_coord==coord_num_primary(1));
    for ilincomb=1:length(lincomb_ptrs)
        lincomb=mtcs.lincomb_list(lincomb_ptrs(ilincomb),:);
        %compute the weight matrix for this
        tok_array=coord_table{1}*lincomb(1);
        for idim=2:nchecks
            tok_array=tok_array+coord_table{idim}*lincomb(idim);
        end
        tok_array=mod(tok_array,ng);
        tok_array=tok_array(:);
        %add up all the probs for each token
        for ig=1:ng
            token_weights(ig)=token_weights(ig)+sum(probs_reshaped(find(tok_array==(ig-1))));
        end
    end
    %now average across all linear combs
    cgs(icg,:)=token_weights/length(lincomb_ptrs);
    lincomb_count=lincomb_count+length(lincomb_ptrs)*length(coord_groups{icg}.coord_num);
end
if ~(lincomb_count==nlc-1)
    warning(sprintf('total number of linear combinations is %1.0f, should be 1.0f',lincomb_count,nlc-1));
end
%example to show complexity of indexing, ng composite, valid multipliers 1 2 4 5 
%ng=6;probs=ones(ng,ng,ng,ng)/ng^4;
%[cgs,mtcs_used]=mtc_probs2cgs(probs,mtc_define(6,gtc_define([0 0;0 1;0 2])));
%mtcs_used.coord_groups{70}
% ans = 
%        name_firstcoord: 'ABE_442'
%                   name: 'ABE_224'
%     coord_num_allmults: [84 132 84 132]
%              coord_num: [132 84]
% mtcs_used.coord_groups_multlist(70,:)
% ans =
%    132    84     0   132    84
% the translations then go through the coordinate groups and average 
% coord=mtcs_used.coord_groups{1}.coord_num(1)
return
