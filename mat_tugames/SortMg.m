function Seff=SortMg(effij,n);
% Sorting the set of most effective
% coalitions with respect to their
% cardinality. Ascent ordering.
% Smallest coalitions are coming first.
   bd=length(effij);
   Seff=zeros(1,bd);
   slck=cell(n,1);
   it=0:-1:1-n;
   PlyMat=rem(floor(effij(:)*pow2(it)),2)==1;
   clsize=(sum(PlyMat,2))';
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
       if n==16
       strmat{k}=base2dec(reshape(dec2base(M2,17),sM2(1),sM2(2)),17);
       else % will not be sorted correctly for base 36!!!
       strmat{k}=base2dec(reshape(dec2base(M2,36),sM2(1),sM2(2)),36);
       end
       [sco{k} ix{k}]=sort(strmat{k});
      id=ii+lc-1;
      Seff(:,[ii:id])=tkc{k}(ix{k});
      clix=[];
      M2=[];
      ii=id+1;
    end
