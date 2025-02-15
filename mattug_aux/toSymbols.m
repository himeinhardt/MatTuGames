function [numb dig]=toSymbols(vec,n);
% TOSYMBOLS converts a vector (vec,n) into digits.
%
%
% Define variables:
%  output:
%  numb    -- A vector vec will be converted into a number of digits or
%             a hexadecimal number. 
%             Examples:
%             A vector [2 4 6] will be converted into the number 246.
%             A vector like [2 4 11] will be converted into the 
%             hexadecimal number 587, since n=11>9.
%  dig     -- Returns a digital number or the hexadecimal value
%             For the above examples:
%             We get the digital number 246
%             and the hex-number 24B respectively.
%
%  input:
%  vec     -- A vector of maximal length of 15. A component cannot be greater
%             than 15. Hence, a vector like [1 2 4 5 16] is not allowed.
%   n      -- Maximal number of players. This number cannot be greater than 15.
%

% Record of revisions:
%   Date         Programmer                Version
%   ==========   ===================       =======
%   08/02/2010   Holger I. Meinhardt       0.1 alpha
%                University of Karlsruhe
%   E-Mail:      Holger.Meinhardt@wiwi.uni-karlsruhe.de


if nargin<1
   error('toSymbols: At least a vector must be given!');
 elseif nargin==1
   n=max(vec);
   if n>35
   error('toSymbols: The size/largest component cannot be larger than 35!');
   end
 else
   n1=max(vec);
   if n1>n
   error('toSymbols: The vector is not consistent with its size!');
   end
end

bd=length(vec);


symbols='123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ';

if n<=9
  for k=1:bd
      dig(k)=symbols(vec(k));
  end 
 else 
     for k=1:bd
       basedig(k)=symbols(vec(k));
       dig(k)=basedig(k);
     end
end

if n<=9
   numb=str2num(dig);
 else
   numb=base2dec(dig,36);
end
