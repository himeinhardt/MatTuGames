function [Seff bmg]=sortsets(effij,n)
% SORTSETS sorts a set of coalitions (coal,n) that reflects the
% unique integer representation of sets used in Matlab
% according to their power set representations, that is,
% the coalitions will be sorted according to their cardinality.
% Smallest coalitions coming first.
%
% Example:
%  A collection of sets based on their unique integer representation
%   S=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]
%  will be sorted according to their power set representation given by
%   S=[1 2 4 8 3 5 6 9 10 12 7 11 13 14 15]
%
% Usage: [Seff bmg]=sortsets(effij,n)
% Define variables:
%  output:
%  Seff      -- Power set/Coalitions sorted with respect to their size.
%  bmg       -- Sorted sets/coalitions given by their binary representation.
%
%  input:
%  coal      -- Sets/Coalitions sorted with respect to their unique integer
%               representation.
%   n        -- Maximal number of players. This number cannot be greater than 35.
%

%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/29/2010        0.1 alpha       hme
%   04/09/2011        0.1 beta        hme
%   10/14/2012        0.3 beta        hme
%

if n>35
   error('sortsets: The power set is too big to be sorted correctly! The number of players n must be less than or equal to 35. Sorry!')
 else
end


   bd=length(effij);
   Seff=zeros(1,bd);
   slck=cell(n,1);
   ov=ones(n,1);
   for k=1:n, PlyMat(:,k)=bitget(effij,k);end
   clsize=PlyMat*ov;
   J=1:n;
   ism=ismember(J,clsize);
   clm=J(ism); % Information needed to sort subsets, otherwise n would be sufficient.
   lcl=length(clm);
   J=J(ones(bd,1),:);
   M=PlyMat.*J;
   M1=cell(lcl,1);
   slck=cell(lcl,1);
   sco=cell(lcl,1);
   ix=cell(lcl,1);
   ii=1;
   for kk=1:lcl
       k=clm(kk);
       slck{k}=ismember(clsize,k);
       tkc{k}=effij(slck{k});
       lc=length(tkc{k});
       M1{k}=sort(M(slck{k},:),2);
       sM=size(M1{k},2);
       M2=M1{k}(:,[sM-k+1:sM]);
       sM2=size(M2);
% We cannot use a base larger than 17, otherwise the set will not be sorted correctly.
       if n<=16
       strmat{k}=base2dec(reshape(dec2base(M2,17),sM2(1),sM2(2)),17);
       else % will not be sorted correctly for base 36!!!
       strmat{k}=base2dec(reshape(dec2base(M2,36),sM2(1),sM2(2)),36);
       end
       [sco{k} ix{k}]=sort(strmat{k});
      id=ii+lc-1;
      Seff(:,[ii:id])=tkc{k}(ix{k});
      if nargout==2
        bmg([ii:id],:)=dec2bin(Seff(:,[ii:id]),n);
      end
      clix=[];
      M2=[];
      ii=id+1;
    end
