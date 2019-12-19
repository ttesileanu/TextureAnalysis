function [cgs_vals,cgs_list]=mtc_cgsstruct2cgs(cgs_struct,mtcs)
% [cgs_vals,cgs_list]=mtc_cgsstruct2cgs(cgs_struct,mtcs) is a utility to convert
% a structure of coordinate groups to an array
%
% coordinate groups that are missing are replaced by unbiased values
%
% mtcs: structure returned by mtc_define
% cgs_struct:  a structure with fields from mtcs.coord_groups{*}.name, each of which 
%   is a vector of size [1 mtcs.ng]
% cgs_vals: an array of size [length(mtcs.coord_groups),ng]
% cgs_list:  a list of the entries that are present in cgs_struct (numeric;
%    the indices into coord_groups)
%
% This does not check for missing or extra fields in cgs_struct
%
%   See also:  MTC_DEFINE, MTC_AUGCOORDS, MTC_PROBS2CGS, MTC_CGSS2PROBS, MTC_CGST2PROBS, MTC_CGS2CGSTRUCT2,
%      MTC_CGSSTRUCT_FILL, MTC_CGSSTRUCT_CLEAN, MTC_UNBIASED.
%
ng=mtcs.ng; %number of gray levels
unbiased=mtc_unbiased(ng);
cgs_list=[];
for icgs=1:length(mtcs.coord_groups)
    cgsname=mtcs.coord_groups{icgs}.name;
    if isfield(cgs_struct,cgsname)
        cgs_vals(icgs,:)=cgs_struct.(cgsname);
        cgs_list=[cgs_list,icgs];
    else
        baseharm=mtcs.coord_groups{icgs}.baseharm;
        cgs_vals(icgs,:)=unbiased(baseharm,:);
    end
end
return

 