function [removed_struct,counts_used]=mtc_cgsstruct_remove(cgs_struct,orders,mtcs,counts)
% [removed_struct,counts_used]=mtc_cgsstruct_remove(cgs_struct,orders,mtcs,counts)
% removes one or more orders from a coordinate group structure
%
% cgs_struct: a coordinate structure
% orders: orders of coordinate groups to remove (a subset of 1 to mtcs.nchecks
% mtcs: an mtc definition, from mtc_define.
% counts:  result of counts=mtc_countparams(mtcs); computed if not supplied
%
% removed_struct:  cgs_struct with requested orders removed
% counts_used: result of counts=mtc_countparams(mtcs)
%
%    See also:  MTC_DEFINE, MTC_COUNTPARAMS, MTC_CGS_REMOVE, MTC_CGSSTRUCT2CGS_ME.
%
if (nargin<=3)
    counts=mtc_countparams(mtcs);
end
counts_used=counts;
removed_struct=cgs_struct;
orders=intersect(orders,1:mtcs.nchecks);
for iorder_ptr=1:length(orders)
    iorder=orders(iorder_ptr);
    cg_list=counts.coord_group_list{iorder};
    for icg=cg_list
        cname=mtcs.coord_groups{icg}.name;
        if isfield(removed_struct,cname)
            removed_struct=rmfield(removed_struct,cname);
        end
    end
end
return
