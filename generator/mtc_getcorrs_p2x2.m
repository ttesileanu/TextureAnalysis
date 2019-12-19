function [corrs,opts_used]=mtc_getcorrs_p2x2(p2x2,opts)
% [corrs,opts_used]=mtc_getcorrs_p2x2_mtc(p2x2,opts) gets entropies and checks normalization in
%  2x2 block probabilities.  Name is kept for compatibility with
%  getcorrs_p2x2, but this does not extract the correlations
% p2x2: an array of size [ng ng ng ng], indicating the probability of each 2x2 block
%   p2x2(ia+1,ib+1,ic+1,id+1) is the probability of [ia ib; ic id]
%
% opts: tolerances, from mtc_defopts.  If an integer (0 or 1), it is nowarn, determining
%   whether there is a warning if probs are negative or > 1 (for compatiblity with getcorrs_p2x2)
%   defaults to 0.  Otherwise this is controlled by opts.nowarn_getcorrs
%
% corrs.norm: normalization (should be 1); sum of values in p2x2
% corrs.pmin: minimum of values in p2x2
% corrs.pmax: maximum of values in p2x2
% corrs.ok:   1 if all values in p2x2 are within range [0, 1]
% corrs.entropy:   entropy of the resulting MRF, assuming existence of MRF consistent with p2x2
% corrs.entropy_area: entropy of the 2x2 probability distributions (not entropy per unit area)
%
% corrs.cig_conds: ng*(ng-1)^2 by 4 arraycorrs.cig_conds(u_ia_ib,k)=extent to which, for a cell in position k
%    that the neighbors of k are conditionally independent
%    there are ng combinations for the cell in question, and (ng-1)*(ng-1)
%    joint values for its neighbors that must be checked  (others are linearly dependent).
%    Pickard conditions: all 0's in column 1 and column 4, or all 0's in columns 2 and 3.
%    These are the cig_conds of getcorrs_p2x2 generalized to ng>2
%
% the following are empty and included just for compatibilitiy
% corrs.alpha: the fourth-order statistic, 1=even, -1=odd
% corrs.gamma: same as luminance bias, 1=all white, -1=all black
% corrs.beta(1): horizontal second-order correlation
% corrs.beta(2): vertical second-order correlation
% corrs.beta(3): diagonal (upper left to lower right) third-order correlation
% corrs.beta(4): diagonal (upper right to lower left) fourth-order correlation
% corrs.theta(1): third-order correlation of checks A and its flankers, B,C
% corrs.theta(2): third-order correlation of checks B and its flankers, A,D
% corrs.theta(3): third-order correlation of checks C and its flankers, A,D
% corrs.theta(4): third-order correlation of checks D and its flankers, B,C
%
% opts_used: options used (see mtc_defopts)
%
%   Example:
%   corrs=getcorrs_p2x2(getp2x2_pabcde(getpabcde_ag(.1,.4)))
%
%   See also:  GETCORRS_P2X2, MTC_DEFINE.
%
corrs.ok=0;
corrs.alpha=[];
corrs.beta=[];
corrs.theta=[];
corrs.gamma=[];
corrs.cig_conds=[];
%
if (nargin<=1) opts=[]; end
if ~isempty(opts) & isnumeric(opts) %nowarn passed as a numeric second argument
    nowarn=opts;
    opts=[];
    opts.nowarn_getcorrs=nowarn;
end
opts=mtc_defopts(opts);
opts_used=opts;
%
if ~(length(size(p2x2))==4)
   if (opts.nowarn_getcorrs==0)
        warning(' p2x2 does not have 4 dimensions.');
   end
   return
end
ng=min(size(p2x2));
if ~(max(size(p2x2))==ng)
   if (opts.nowarn_getcorrs==0)
       warning(' p2x2 is not of size [ng ng ng ng]');
   end
   return
end
p2x2r=reshape(p2x2,1,ng^4);
corrs.norm=sum(p2x2r);
corrs.pmin=min(p2x2r);
corrs.pmax=max(p2x2r);
corrs.ok=1;
if (corrs.pmin<-opts.tol_nonneg)
   if (opts.nowarn_getcorrs==0)
       warning(' p2x2 has values < 0');
   end
   corrs.ok=0;
end
if (corrs.pmax>1+opts.tol_nonneg)
   if (opts.nowarn_getcorrs==0)
       warning(' p2x2 has values > 1');
   end
   corrs.ok=0;
end
%
%calculate overall entropy
pv=reshape(p2x2,1,ng^4)/corrs.norm;
pvnz=pv(find(pv>0));
corrs.entropy_area=-sum(pvnz.*log(pvnz))/log(2);
%turn this into entropy per unit area
p2x2norm=p2x2/corrs.norm;
for id=1:2
   if (id==1)
       pv2=squeeze(sum(sum(p2x2norm,2),4));
   end
   if (id==2)
       pv2=squeeze(sum(sum(p2x2norm,3),4));
   end
   pvnz=pv2(find(pv2>0));
   ent2x1(id)=-sum(pvnz.*log(pvnz))/log(2);
end
pv1=squeeze(sum(sum(sum(p2x2norm,2),3),4));
pvnz=pv1(find(pv1>0));
ent1x1=-sum(pvnz.*log(pvnz))/log(2);
%
corrs.entropy=corrs.entropy_area-sum(ent2x1)+ent1x1;
%
% calculate something equivalent to conditional independence of neighbors of k
for iu=1:ng
   for k=1:4
      if (k==1) p_ku=squeeze(sum(p2x2(iu,:,:,:),4)); end
      if (k==2) p_ku=squeeze(sum(p2x2(:,iu,:,:),3)); end
      if (k==3) p_ku=squeeze(sum(p2x2(:,:,iu,:),2)); end
      if (k==4) p_ku=squeeze(sum(p2x2(:,:,:,iu),1)); end
      for ia=1:ng-1
          for ib=1:ng-1
              corrs.cig_conds((iu-1)*(ng-1)^2+(ia-1)*(ng-1)+ib,k)=p_ku(ia,ib)*sum(sum(p_ku))-sum(p_ku(ia,:))*sum(p_ku(:,ib)); 
          end
      end
    end
end
return
