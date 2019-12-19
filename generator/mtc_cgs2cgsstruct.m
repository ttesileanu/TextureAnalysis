function cgs_struct=mtc_cgs2cgsstruct(cgs_vals,mtcs)
% cgs_struct=mtc_cgs2cgsstruct(cgs_vals,mtcs) is a utility to convert
% an array of coordinae groups to a structure of coordinate groups
%
% mtcs: structure returned by mtc_define
% cgs_vals: an array of size [length(mtcs.coord_groups),ng]
% cgs_struct:  a structure with fields from mtcs.coord_groups{*}.name, each of which 
%   is a vector of size [1 mtcs.ng]
%
% This does not check for missing or extra fields in cgs_struct
%
%   See also:  MTC_DEFINE, MTC_AUGCOORDS, MTC_PROBS2CGS, MTC_CGSS2PROBS, MTC_CGST2PROBS, MTC_CGSSTRUCT2CGS.
%
ng=mtcs.ng; %number of gray levels
cgs_struct=struct();
for icgs=1:length(mtcs.coord_groups)
    cgsname=mtcs.coord_groups{icgs}.name;
    cgs_struct.(cgsname)=cgs_vals(icgs,:);
end
return

 