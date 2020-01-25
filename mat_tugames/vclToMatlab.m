function [v sS]=vclToMatlab(clm,vlm)
% VCLTOMATLAB computes from a cell, which contains the coalition
% representation in generic power set representation (e.g. Mathematica), 
% and a valuation vector for these coalition, a game v as well as 
% the corresponding unique integer representation of coalitions in an
% unsorted order (Matlab format). 
%
% Usage: [v sS]=vclToMatlab(clm,vlm);
%
% Define variables:
%  output:
%  v        -- The generated TU game of length 2^n-1.    
%  sS       -- The list of coalitions in unique integer representation (unsorted).
%  input:
%  clm      -- The cell or matrix, which contains the coalition information
%              in generic power set representation (e.g. Mathematica format).   
%  vlm      -- The valuation vector of coalitions in the order as presented by
%              the cell input clm. 
%    
%
%    
% Examples:
% A set of coalitions must be coded as a cell, since the grand coalition
% must be listed. It is a absolutely necessity to enlist the grand coalition. 
% clm={[1 2 5 6] [1 3 4], [1 6 7], [2 3 6 7], [2 4 7], [2 5 7], [4 5 6
% 7], [1 2 3 4 5 6 7]};
% 
% The valuation of the coalitions must be given in the same order as for
% the coalitions above. Hence,
% vlm=[10 5 4 12 4 6 15 25];  The imput format must be a vector.
% 
%  [v sS]=vclToMatlab(clm,vlm);
% which returns 
% the game v of size 127.
% and sS =
%
%    51    13    97   102    74    82   120   127
%
% in Matlab's coalition representation (unsorted).    
%     
   
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   06/02/2012        0.2 beta        hme
%                    

vsz=length(vlm);    
csz=numel(clm);
if vsz~=csz
   error('The size of the input variables are not identical!')
end

if iscell(clm)
   clm2=cell(1,csz);
   for k=1:csz
       sS(k)=sum(2.^(clm{k}-1));
   end
else
  error('The input format must be a cell.!')
end

n=numel(clm{end});
N=2^n-1;
S=1:N;
v=zeros(1,N);
v(sS)=vlm;
