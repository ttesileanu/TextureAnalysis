function cgs_reduced=mtc_cgst_reduce(cgs,mtcs)
% cgs_reduced=mtc_cgst_reduce(cgs,mtcs) reduces a set of coordinate group values
% to a set suitable for Fourier inversion.  This uses the convolution by
% the generalized totient function in graylevs_mtc.doc.
%
% This is a nontrivial operation if the number of
% gray levels is composite, since a set of coordinate group values will
% contain the values from the higher harmonics
%
% The mean is not subtracted
%
% cgs: [mtcs.coord_group_count,mtcs.ng]: the total probability for 
%    each value of each linear combination, based on the name for the coordinate group
%    every entry of sum(cgs(:,2)) be the same, typically 1
% mtcs:  a structure, typically created by MTC_DEFINE. 
%  mtcs.ng is number of gray levels
%
% cgs_reduced: same sized as cgs; the higher harmonics have been removed
%
% See also:  MTC_DEFINE, MTC_CGST2PROBS, TOTIENT, ,MTC_CGSS_REDUCE.
%
if (nargin<=1)
    mtcs=[];
end
ng=mtcs.ng;
factors=mtcs.ng_factors;
cgs_reduced=cgs;
for icgs=1:mtcs.coord_groups_count
    ilincomb=mtcs.coord_indices(mtcs.coord_groups{icgs}.coord_num(1));
    period=mtcs.lincomb_period(ilincomb,:);
    period_index=find(mtcs.period_list==period);
    totient_modified=((1+mtcs.totients_spaced(:,period_index))/ng)';
    %mtcs.coord_groups{icgs}.name
    %totient_modified
    convolution=conv([cgs(icgs,:) cgs(icgs,:)],totient_modified);
    cgs_reduced(icgs,:)=convolution(ng+[1:ng]);
end %icgs
return
