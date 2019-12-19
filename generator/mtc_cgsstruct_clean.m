function [cleaned_struct,opts_used]=mtc_cgsstruct_clean(cgs_struct,mtcs,opts)
% [cleaned_struct,opts_used]=mtc_cgsstruct_clean(cgs_struct,mtcs,opts) cleans (removes)
% all unbiased coordinate groups from a structure
%
% cgs_struct: a coordinate structure
% mtcs: result of mtc_defopts
% opts: a set of mtc options, including tolerances, taken from mtc_defopts if not supplied
% 
% cleaned_struct:  a coordinate structure with all coords of cgs_struct removed if they differ from the unbiased value by < opts.tol_match
% opts_used: options used
%
% For compatibility with other modules, the 'structure with empty fields' is converted to []
% This should be the inverse of MTC_CGSSTRUCT_FILL.
% This should have no effect on a conversion to a cgs array by MTC_CGSSTRUCT2CGS.
%
%    No checking for normalization or length of fields or fields that are not legitimate coordinate group names.
%
%    See also:  MTC_DEFOPTS, MTC_CGS2CGSSTRUCT, MTC_CGSSTRUCT_FILL, MTC_UNBIASED.
%
if (nargin<=2)
    opts=[];
end
opts=mtc_defopts(opts);
opts_used=opts;
cleaned_struct=cgs_struct;
if (isempty(cgs_struct))
    return
end
unbiased_list=mtc_unbiased(mtcs.ng);
cnames=fieldnames(cgs_struct);
for iname=1:length(cnames)
    cname=cnames{iname};
    cnum=mtcs.coord_group_ptrs.(cname);
    baseharm=mtcs.coord_groups{cnum}.baseharm;
    unbiased=unbiased_list(baseharm,:);
    if max(abs(cgs_struct.(cname)-unbiased))<opts.tol_match
        cleaned_struct=rmfield(cleaned_struct,cname);
    end
end
if isempty(fieldnames(cleaned_struct))
    cleaned_struct=struct();
end
return

