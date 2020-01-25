function [bg_coord, bg]=p_basis_coordinates(clv)
% P_BASIS_COORDINATES computes a new basis and its corresponding
% coefficients of game v using Matlab's PCT.
% 
% Usage: [bg_coord bg]=p_basis_coordinates(clv)
%

% Define variables:
%  output:
%  bg_coord -- The coefficients of the basis. 
%  bg       -- A basis of an n-person game (Novak/Radzik 1994 IJGT). 
%
%  input:
%  clv      -- TuGame class object.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/26/2013        0.4             hme
%                

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;

int=1-n:1:0;
bg=zeros(N);


parfor k=1:N
   sS=SubSets(k,n);
   ls=length(sS);
   b=zeros(1,ls);
   mat=rem(floor(sS(:)*pow2(int)),2);
   cS=mat*ones(n,1);
   ec=cS(end);
   for jj=1:ls
       b(jj)=(factorial(ec))/(factorial(cS(jj))*factorial(ec-cS(jj)));
   end   
   temp=zeros(1,N);
   temp(sS)=b.^(-1); 
   bg(k,:)=temp;
end

sutm=sparse(bg);
bg_coord=(sutm\v')';
%bg_coord=(inv(sutm)*v')';


