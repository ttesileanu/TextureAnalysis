function [removed_struct,rm_list,rm_indices]=mtc_cgsstruct_rmclass(cgs_struct,checks,mtcs)
% [removed_struct,rm_list,rm_indices]=mtc_cgsstruct_remove(cgs_struct,checks,mtcs,counts)
% removes a coordinate class from a coordinate group structure
%
% cgs_struct: a coordinate structure
% checks: specification of class to remove, as in 'AC', which will remove all AC_1_1, AC_1_2, ...
%    checks can also be a cell array, in which case all of the listed classes are removed
% mtcs: an mtc definition, from mtc_define.
%
% removed_struct:  cgs_struct with requested coordinate groups removed
% rm_list: list of coord groups removed
% rm_indices: pointers to coord groups removed
%
%    See also:  MTC_DEFINE, MTC_CGSSTRUCT_REMOVE, MTC_COORD_RENAME, MTC_CGS_RMCLASS, MTC_AUGCOORDS.
%
removed_struct=cgs_struct;
rm_list=[];
rm_indices=[];
if iscell(checks)
    whichsubset=[];
    for ip=1:length(checks)
        whichsubset(ip)=strmatch(checks{ip},mtcs.checks,'exact');
    end
else
    whichsubset=strmatch(checks,mtcs.checks,'exact');
end
if isempty(whichsubset)
    return
end
for icg=1:length(mtcs.coord_groups)
    subset=mtcs.coord_subset(mtcs.coord_groups{icg}.coord_num(1));
    if (ismember(subset,whichsubset))
        removed_struct=rmfield(removed_struct,mtcs.coord_groups{icg}.name);
        rm_list{end+1}=mtcs.coord_groups{icg}.name;
        rm_indices(end+1)=icg;
    end
end

    
 