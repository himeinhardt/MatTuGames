function [sS,sidx]=clToMatlab(clm)
% CLTOMATLAB computes from a cell or matrix, which contains the coalition
% representation in generic power set representation (e.g. Mathematica), 
% into the corresponding unique integer representation of coalitions 
% (Matlab format). 
%
% Usage: sS=clToMatlab(clm);
%
% Define variables:
%  output:
%  sS       -- The list of coalitions in unique integer representation.
%  sidx     -- Inital position of coalitions. 
%  input:
%  clm      -- The cell or matrix, which contains the coalition information
%              in generic power set representation (e.g. Mathematica format).   
%
%    
% Examples:
% Coalitions of unequal size must be represented as a cell input:    
% clm={[1 2 5 6] [1 3 4], [1 6 7], [2 3 6 7], [2 4 7], [2 5 7], [4 5 6 7]};
%
%   sS=clToMatlab(clm)
% which returns sS =
%
%    13    51    74    82    97   102   120    
%
% in Matlab's coalition representation.    
%     
% Coalitions of equal size can be represented as a matrix
% clmat=[1 2 5; 1 3 4; 1 6 7; 2 3 6; 2 4 7; 2 5 7; 4 5 6;]
%  
%    sS2=clToMatlab(clmat);
% which returns sS2 =
%
%    13   19   38   56   74   82   97    
%
% in Matlab's coalition representation.
%
% Even a vector of players can be inputed, e.g. clv=[1   2   5   6]; 
%
% sS3=clToMatlab(clv);
% which returns sS3 = 
% 
%  51. 
% 
% Hence, cl51=dec2bin(51,6)
% cl51 = 110011.    
% 
%    
    
    
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/30/2012        0.2 beta        hme
%   01/30/2017        0.9             hme
%                    
    
if iscell(clm)
   csz=numel(clm);
   clm2=cell(1,csz);
   for k=1:csz
       sS(k)=sum(2.^(clm{k}-1));
   end
   [sS,sidx]=sort(sS);
 elseif ismatrix(clm) 
  clm=clm-1;
  clm=2.^clm;
  cl=sum(clm,2)';
  [sS,sidx]=sort(cl);
else
  error('The input format must either be a cell or matrix!')
end
