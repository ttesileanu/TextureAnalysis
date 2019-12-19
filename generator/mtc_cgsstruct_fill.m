function [filled_struct,cgs_list]=mtc_cgsstruct_fill(cgs_struct,mtcs)
% [filled_struct,cgs_list]=mtc_cgsstruct_fill(cgs_struct,mtcs) fills in all missing fields in a coordinate structure
%
% cgs_struct: a coordinate structure
% mtcs: an mtc definition, from mtc_define
% 
% filled_struct:  a coordinate structure with all coords of cgs_struct filled with 1/ng if they had been blank
% cgs_list:  a list of the entries that are present in cgs_struct (numeric;
%    the indices into coord_groups)
%
% For compatibility with other modules, the 'structure with empty fields' is converted to []
% This should be the inverse of MTC_CGSSTRUCT_CLEAN.
% This should have no effect on a conversion to a cgs array by MTC_CGSSTRUCT2CGS.
%
%    No checking for normalization or length of fields or fields that are not legitimate coordinate group names.
%
%    See also:  MTC_DEFINE, MTC_CGS2CGSSTRUCT, MTC_CGSSTRUCT_CLEAN, MTC_UNBIASED.
%
ng=mtcs.ng;
unbiased_list=mtc_unbiased(ng);
filled_struct=[];
cgs_list=[];
for icgs=1:length(mtcs.coord_groups)
    cname=mtcs.coord_groups{icgs}.name;
    baseharm=mtcs.coord_groups{icgs}.baseharm;
    if isfield(cgs_struct,cname)
        filled_struct.(cname)=cgs_struct.(cname);
        cgs_list=[cgs_list,icgs];
    else
        filled_struct.(cname)=unbiased_list(baseharm,:);
    end
end
return
