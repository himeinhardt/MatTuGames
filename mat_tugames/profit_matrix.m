function [prfm buy_mat sel_mat]=profit_matrix(val_buy,val_sel)
% PROFIT_MATRIX computes from a vector/matrix of buyers and sellers valuation
% a symmetric profit (assignment) matrix for a symmetric assignment problem.
% If the vectors/matrices  are not symmetric, they will be transformed into 
% symmetric ones.
%
% Usage: [prfm buy_mat sel_mat]=profit_matrix(val_buy,val_sel)
% Define variables:
%  output:
%  prfm     -- A profit matrix for a symmetric assignment problem.
%  buy_mat  -- The symmetric buyers valuation matrix.
%  sel_mat  -- The symmetric sellers valuation matrix.
%
%  input:
%  val_sl   -- A vector/matrix of buyers valuation.
%  val_sl   -- A vector/matrix of sellers valuation.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   08/17/2010        0.1 beta        hme
%                



ls=size(val_sel);
lb=size(val_buy);
dm=max(ls(2),lb(2));

Zmat=zeros(dm);

if isvector(val_sel)
  if ls(2)<dm
    val_sel(1,dm)=0;
  end
  sel_mat=repmat(val_sel,dm,1)';
 else
   if ls(2)<dm
     val_sel(dm,dm)=0;
   end
   sel_mat=val_sel;
end

if isvector(val_buy)
   if lb(2)<dm
     val_buy(1,dm)=0;
   end
   buy_mat=repmat(val_buy,dm,1);
 else
   if lb(2)<dm
     val_buy(dm,dm)=0;
   end
   buy_mat=val_buy;
end

prfm1=buy_mat-sel_mat;
prfm=max(prfm1,Zmat);

