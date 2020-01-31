function [imp_vert crst]=AntiImputationVertices(v,method,tol)
% ANTIIMPUTATIONVERTICES computes all vertices of the anti imputation set of game v, 
% whenever the anti-imputation set is  essential. The cdd-library by 
% Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% Note: Windows users must port the shell script 'corevert' to get full 
% operationality.
%
% Usage: [imp_vert crst]=AntiImputationVertices(v,method)
% Define variables:
%  output:
%  imp_vert   -- Matrix of anti-imputation vertices. Output is numeric or a string.
%  crst       -- The anti-imputation constraints.
%  input:
%  v          -- A Tu-Game v of length 2^n-1. 
%  method     -- A string to call a method from the cdd-library.
%                Permissible methods are: 
%                'float' that is, results are given by real numbers.
%                Default is 'float'.
%                'gmp' that is, results are given by rational numbers.
%                Choose this method whenever the result with 'float'
%                is not as expected. This method needs more time to complete. 
%  tol        -- A positive tolerance value. Its default value is set to 10^9*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   07/16/2015        0.7             hme
%   05/15/2019        1.1             hme
%                

N=length(v);
[~, n]=log2(N);
if (2^n-1)~=N
   error('Game has not the correct size!');
end


if nargin<2
  method='float';
  tol=10^9*eps;
 elseif nargin<3
  tol=10^9*eps;
end


if strcmp(method,'gmp')
elseif strcmp(method,'float')
 else
   error('Method not recognized. Permissible methods are: float or gmp!');
end

J=1:n;
S=bitset(0,J);
Nk=N-2.^(J-1);
vi=v(S);
s_vi=vi*ones(n,1);
if s_vi<=v(N)
  error('Game is inessential!');
end
v=v+tol;
Nk(end+1)=N;
lS=length(Nk);
v1=v(Nk);
it=0:-1:1-n;
PlyMat=rem(floor(Nk(:)*pow2(it)),2);
%
% Defining imputation set.
%
PlyMat(lS+1,:)=-PlyMat(end,:);
v1(lS+1)=-v1(lS);

% Defining format
crst=zeros(lS+1,n+1);
crst(:,2:n+1)=-PlyMat;
crst(:,1)=v1';
sz=size(crst);  
if strcmp(method,'gmp')
  crst=rats(crst);
  if strcmp(crst(1,2:8),'Columns')
   crst(1:2,:)=[];
   crst=num2str(crst);
  end
  rsz=size(crst);
end


% Defining auxiliary files.
str='impvert_aux';
inefile=strcat(str,'.ine');
datafile=strcat(str,'.dat');
if strcmp(method,'gmp')
 format=repmat('%3s ',1,rsz(2));
 format=strcat(format,'%3s\r');
else
 format=repmat('%3g ',1,sz(2)-1);
 format=strcat(format,' %3g\n');
end

%
% Writing into and reading from files.
%
[fid,~]=fopen(inefile,'w');
if fid>0
fprintf(fid,'%s\n','H-representation');
fprintf(fid,'%s\n','begin');
if strcmp(method,'gmp')
  fprintf(fid,'%5d %5d %5s\n',sz(1),sz(2),'rational');
elseif strcmp(method,'float')
  fprintf(fid,'%5d %5d %5s\n',sz(1),sz(2),'real');
 else
  fprintf(fid,'%5d %5d %5s\n',sz(1),sz(2),'rational');
end
fprintf(fid,format,crst');
fprintf(fid,'%s\n','end');
fclose(fid);
 else
   disp('Output file open failed');
end

% Calling the external library cdd.
if strcmp(method,'gmp')
!corevert -m gmp impvert_aux.ine
elseif strcmp(method,'float')
!corevert -m float impvert_aux.ine
else
!corevert -m gmp impvert_aux.ine
end

%
% Starting of section reading from file and reformatting of the results.
% Increase the 'et' value below under method 'gmp' if the formatting was 
% not correct. In case that you encounter a memory problem decrease 
% this value. But then the output might not be anymore formatted correctly,
% and some output will be lost.
%

% Reading from file
fk=factorial(n);
[fid,~]=fopen(datafile,'r');
if fid>0
   if strcmp(method,'gmp')
      line=fgetl(fid);
      ls=size(line);
      et=ls(2)+1;
      [imp_vert, ~]=fscanf(fid,'%c',[et*fk (n+1)]);
       imp_vert=imp_vert';
   else
      [imp_vert count]=fscanf(fid,'%f',[n+1 fk]);
      if count==0
        error('No vertex computed! Use method gmp instead.');
       else
      end
   end
  fclose(fid);
 else
  error('Output file open failed');
end
imp_vert=fliplr(imp_vert);
% Reformating the output
if strcmp(method,'gmp')
  if isempty(imp_vert)
%   imp_vert=line;
   imp_vert=str2num(imp_vert);
  else
   line=str2num(line);
   imp_vert=str2num(imp_vert);
   imp_vert=[line;imp_vert];
% Converting string into real numbers in order to reformat output. 
%  imp_vert=str2num(imp_vert);
% Converting real numbers into a rational number approximation represented 
% as a string.
  imp_vert=rats(imp_vert); 
  imp_vert(:,1:11)=[]; % deleting the first column vector of ones.
  end
else
  imp_vert=imp_vert';
  imp_vert(:,1)=[]; % deleting the first column vector of ones.
end
