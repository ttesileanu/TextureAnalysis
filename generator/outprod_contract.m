function [u,v_sum,w_sum]=outprod_contract(v,w,vwcd)
% [u,v_sum,w_sum]=outprod_contract(v,w,vwcd) computes the "contracted" outer product of v and w.
%
% v, w: matrices, two or more dimensions.  Avoid trailing length-1 dimensions
% vwcd: vcd=vwcd(1,:), wcd=vwcd(2,:), where vcd, wcd: lists from [1:dim(v)] and [1:dim(w)],
%   and indicate the dimensions to contract over
% u: an array of dimension dim(v)+dim(w)-size(vwcd,2), size = [non-contracted v, contracted dimensions, non-contracted w]
% v_sum: sum of v over its contracted dimensions
% w_sum: sum of w over its contraxted dimensions
%
% Dimensions of length 1 are not ignored.
% The dimensions vcd of v and wcd of w must match in length.
% The first ndims(v)-length(vcd) dimensions of u match the corresponding dimensions of w
% The next length(vcd)=length(wcd) dimensions of u are the paired dimensions
% The final ndims(w)-length(wcd) dimensions of u match the corresponding dimensions of w.
%
% For u, on the contracted dimensions, entries are multiplied
%
% size(u)=ndims(v)+ndims(w)-size(vwcd,2)
%
% If vcd and wcd are empty, then u=outprod(v,w)
% If v and w are matrices and vwcd=[2;1], then squeeze(sum(u,2))= ordinary matrix product of v and w
%
% Note: this can be used for Markov calculations, e.g,
% p(ABC)=p(AB)p(BC)/P(B), or p(ABE)=p(ABCD)p(CDE)/p(CD) etc.  If the denonimator
%   can be assumed constant, this is proportional to u.  If not, the denominator
%   is given by v_sum or u_sum (should be consistent)
%
%   See also:  OUTPROD, MTC_AUGCOORDS.
%
if (nargin<=2)
    vwcd=[];
end
if (isempty(vwcd))
    vcd=[];
    wcd=[];
else
    vcd=vwcd(1,:);
    wcd=vwcd(2,:);
end
%
dv=length(size(v));
v_unpaired=setdiff([1:dv],vcd);
v_perm=permute(v,[v_unpaired vcd]);
v_perm_size=size(v_perm);
v_unpaired_d=v_perm_size(1:length(v_unpaired));
v_unpaired_tot=prod(v_unpaired_d);
v_paired_d=v_perm_size(length(v_unpaired)+[1:length(vcd)]);
v_paired_tot=prod(v_paired_d);
v_reshaped=reshape(v_perm,[v_unpaired_tot v_paired_tot]);
v_reshaped_size=size(v_reshaped);
%
dw=length(size(w));
w_unpaired=setdiff([1:dw],wcd);
w_perm=permute(w,[wcd w_unpaired]);
w_perm_size=size(w_perm);
w_paired_d=w_perm_size(1:length(wcd));
w_paired_tot=prod(w_paired_d);
w_unpaired_d=w_perm_size(length(wcd)+[1:length(w_unpaired)]);
w_unpaired_tot=prod(w_unpaired_d);
w_reshaped=reshape(w_perm,[w_paired_tot w_unpaired_tot]);
w_reshaped_size=size(w_reshaped);
%
if length(v_paired_d)~=length(w_paired_d) %this should never happen
    error(sprintf('number of paired dimensions, %1.0f and %1.0f, does not match.',length(v_paired_d),length(w_paired_d)));
    return
end
if ~all(v_paired_d==w_paired_d)
    disp(v_paired_d)
    disp(w_paired_d)
    error('lengths of paired dimensions does not match.');
    return
end
%
%v_paired_tot should equal w_unparied_tot
u_raw=zeros([v_unpaired_tot v_paired_tot w_unpaired_tot]);
for ipair=1:v_paired_tot
    u_raw(:,ipair,:)=reshape(v_reshaped(:,ipair)*w_reshaped(ipair,:),[v_unpaired_tot 1 w_unpaired_tot]);
end
u=reshape(u_raw,[v_unpaired_d,w_paired_d,w_unpaired_d]);
vw_raw=sum(sum(u_raw,3),1);
if length(w_paired_d)<=1
    vw=vw_raw(:);
else
    vw=reshape(vw_raw,w_paired_d);
end
if (nargout>1)
    v_sum_raw=sum(v_reshaped,1);
    if length(v_paired_d)<=1
        v_sum=v_sum(:);
    else
        v_sum=reshape(v_sum_raw,v_paired_d);
    end
    w_sum_raw=sum(w_reshaped,2);
    if length(w_paired_d)<=1
        w_sum=w_sum(:);
    else
        w_sum=reshape(w_sum_raw,w_paired_d);
    end
end
return


