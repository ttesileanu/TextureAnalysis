function coord_table=mtc_make_coord_table(ng,nchecks)
% coord_table=mtc_make_coord_table(ng,nchecks) is a utility to
% make a cell array, whose idim-th entry of size ng^nchecks,
% and each of these is a multidimensional array whose values are the idim-th coordinate
%
%  See also:  MTC_PROBS2CGS, MTC_CGS2PROBS.
%
coord_table=cell(nchecks,1);
if nchecks==1
    coord_table{1}=[0:ng-1]';
else
    for idim=1:nchecks
        coord_v=[0:ng-1];
        reshape_v=ones(1,nchecks);
        reshape_v(idim)=ng;
        coord_v=reshape([0:ng-1],reshape_v);
        repmat_v=repmat(ng,1,nchecks);
        repmat_v(idim)=1;
        coord_table{idim}=repmat(coord_v,repmat_v);
    end
end
return

