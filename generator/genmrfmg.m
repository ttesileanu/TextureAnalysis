function [img,stats,optsused]=genmrfmg(opts,area)
% [img,stats,optsused]=genmrfmg(opts,area) creates a noise image
%       via a general Markov random field, not necessarily binary
%
%      CAUTION
% Generating the block probabilities does not guarantee that the Markov texture
% has the right statistics -- see eobmrf.
%
% opts.pblocks: an array of size [ng ng ng ng], indicating the probability of each 2x2 block
%   opts.pblocks(ia+1,ib+1,ic+1,id+1) is the probability of [ia ib; ic id]
%   sum(sum(sum(sum(opts.pblocks))))=1; homogeneity conditions assumed to hold
%
% if opts IS a structure:
% opts.show: set to 1 to draw the image
% opts.statsub:  number of statistical subdivisions on each axis (defaults to 2, must be >0)
% opts.firstcol: holds the first column, filled randomly if not present, length(opts.firstcol)=area(1)
% opts.firstrow: holds the first row, filled randomly if not present, length(opts.firstrow)=area(2)
%   must have opts.firstcol(1)=opts.firstrow(1)
% luminance bias and vertical and horizontal pairwise correlations are calculated from pblocks
% opts.noftps: 1 to suppress calculation of Fourier transform parameters (much faster for ng >= 5)
%
% if opts IS NOT a structure, then it is the array pblocks referred to above, and opts.show is set to 1
%
%  Block probabilities, and (unless Pickard or CIG conditions are met) may be incompatible
%  See  Champagnat,F., Idier,J., and Goussard, Y. IEEE Trans. Inf. Theory 44, 1998, 2901-2916.
%
% following combinations of arguments may be omitted:
%
% area: size of img generated, defaults to [256 256]
%
% uses recursive nonvectorized algorithm, and is therefore slow
%
% img: the image; size(img) = area; values are [0, 1, ..., ng-1]
% stats: some image statistics
% optsused: the options used -- the following are calculated based on block probabilities:
% optsused.pixelprobs(ig,1): probabilities of the pixel value ig
% optsused.dyadprobs(ig1,ig2): probabilities of the joint pixel values ig1, ig2 (supplanted by dyadprobs_[row|col]
% optsused.dyadprobs_row(ig1,ig2): probabilities of the joint pixel values ig1, ig2 in a row (22-Sep-12)
% optsused.dyadprobs_col(ig1,ig2): probabilities of the joint pixel values ig1, ig2 in a column (22-Sep-12)
% optsused.ftps: Fourier transform parameters, from GETFTPS_P2X2
%
% modified 22 Sep 12 to take into account that statistics along first row and first column
% may differ.  optsused.dyadprobs_row and optsused.dyadprobs_col
% introduced.  optsused.dyadprobs retained for backward compatibility but not used here.
%
% modified 25 Jul 14 to add noftps=1

%   See also GETFTPS_P2X2, GENMRFMG_1D, GENMRFM, MTC_PLANE_DEMO.
%
if (~isstruct(opts))
   optsm.pblocks=opts;
   optsm.show=1;
else
   optsm=opts;
end
optsm=filldefault(optsm,'noftps',0);
if (optsm.noftps==0)
    optsm.ftps=getftps_p2x2(optsm.pblocks);
else
    optsm.ftps='not requested';
end
ng=size(optsm.pblocks,1);
%
optsm.pixelprobs=reshape(sum(sum(sum(optsm.pblocks,2),3),4),[1 ng]);
optsm.dyadprobs_row=sum(sum(optsm.pblocks,3),4);
optsm.dyadprobs_col=squeeze(sum(sum(optsm.pblocks,2),4));
optsm.dyadprobs=optsm.dyadprobs_row; %for backward compatibility
%
ard=[256 256];
if (nargin >= 2)
   if (length(area) >= 2)
      ard=area;
   else
      ard=[area(1) area(1)];
   end
end
upx=ard;
if (~isfield(optsm,'show')) optsm.show=0; end
if (~isfield(optsm,'statsub')) optsm.statsub=2; end
if (~isfield(optsm,'firstcol')) optsm.firstcol=[]; end
if (~isfield(optsm,'firstrow')) optsm.firstrow=[]; end
if (length(optsm.firstcol)~=upx(1)) optsm.firstcol=[]; end
if (length(optsm.firstrow)~=upx(2)) optsm.firstrow=[]; end
%
if (0==0); %recursive method in this case
   if (length(optsm.firstrow)==upx(2))
      imgy=reshape(round(optsm.firstrow),[1,upx(2)]);
   else
      %imgy=(rand(1,upx(2))<=(1+optsm.bias_lum_init)/2);
      %this allows for correlated first row
      imgy=genmrfmg_1d(upx(2),optsm.pixelprobs,optsm.dyadprobs_row);
	end
   if (length(optsm.firstcol)==upx(1))
      imgx=reshape(round(optsm.firstcol),[upx(1),1]);
   else
      %this allows for correlated first column
      imgx=genmrfmg_1d(upx(1),optsm.pixelprobs,optsm.dyadprobs_col,imgy(1,1)); %imgy added 22-Sep-12
      imgx=imgx';
   end
   ifix=0;
   if ((length(optsm.firstrow)>0) & (length(optsm.firstcol)==0)) imgx(1)=imgy(1); end
   if ((length(optsm.firstrow)==0) & (length(optsm.firstcol)>0)) imgy(1)=imgx(1); end
   if ((length(optsm.firstrow)>0) & (length(optsm.firstcol)>0)) 
      if (imgx(1)~=imgy(1))
         error('specified first row and first column are incompatible.')
         return
      end
   end
   img=rand(upx);
   img(1,:)=imgy;
   img(:,1)=imgx;
   for ic=2:upx(1)
      for jc=2:upx(2)
         probs=squeeze(optsm.pblocks(img(ic-1,jc-1)+1,img(ic-1,jc)+1,img(ic,jc-1)+1,:));
         %calculate conditional probabilities
         if sum(probs)==0
            img(ic,jc)=sum(rand(1)>cumsum(optsm.pixelprobs));
         else
            img(ic,jc)=sum(rand(1)>(cumsum(probs)/sum(probs)));
         end
         img(ic,jc)=min(img(ic,jc),ng-1);
      end
   end
end
%calculate stats
stats.pixelprobs=getstats(img,ng);
stats.pixelprobs_sub=cell(optsm.statsub,optsm.statsub);
for isub=1:optsm.statsub
   ilo=max(1,round(upx(1)*(isub-1)/optsm.statsub));
   ihi=round(upx(1)*isub/optsm.statsub);
   for jsub=1:optsm.statsub
	   jlo=max(1,round(upx(2)*(jsub-1)/optsm.statsub));
   	  jhi=round(upx(2)*jsub/optsm.statsub);
      stats.pixelprobs_sub{isub,jsub}=getstats(img(ilo:ihi,jlo:jhi),ng);
   end
end
%plot if requested
if (optsm.show>0)
   figure;
   imagesc(img,[0 ng-1]);
   axis equal;axis off;colormap('gray');
end
%
optsused=optsm;
return
%
function pixelprobs=getstats(img,ng)
for ig=1:ng
    pixelprobs(1,ig)=sum(sum(img==(ig-1)))/prod(size(img));
end
return
