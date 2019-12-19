function factors=flfs_bc_factors(configs)
%
% determine the factors by which the b-type (horizontal 1x2) and c-type
% (vertical 1x2) correlations appear in a falling-sticks model in a 2x2
% block, with stable probabilities determined by eivecs
%
% configs: output of flfs_enumerate for the 2x2 configuration, rows and columns flfs setup
%
% each of the following is a length(configs)x1 vector
% factors.b:  fraction of 2x2 region covered by a complete row in stable probability distribution
% factors.c:  fraction of 2x2 region covered by a complete col in stable probability distribution
% factors.bb: fraction of 2x2 region covered by a pair of rows
% factors.cc: fraction of 2x2 region covered by a pair of rows
%
% eventually this could be recoded in a much more general way
%
%  See also: FLFS_DEFINE, FLFS_ENUMERATE, FLFS_FINDEIVEC, FLFS_BC_FACTORS, FLFS_BC_FINDP_DEMO.
%
nconfigs=length(configs);
factors.b=zeros(nconfigs,1);
factors.c=zeros(nconfigs,1);
factors.bb=zeros(nconfigs,1);
factors.cc=zeros(nconfigs,1);
for iconfig=1:length(configs)
    s=configs(iconfig).covered_cells;
    na=double(sum(s=='a')==2);
    nb=double(sum(s=='b')==2);
    nc=double(sum(s=='c')==2);
    nd=double(sum(s=='d')==2);
    factors.b(iconfig)=(na+nb)/2;
    factors.c(iconfig)=(nc+nd)/2;
    factors.bb(iconfig)=na*nb;
    factors.cc(iconfig)=nc*nd;
end
return
