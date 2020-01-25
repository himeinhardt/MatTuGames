function [core_vert crst]=AntiCoreVertices(clv,method,tol)
% ANTICOREVERTICES computes all anti-core vertices of game v, 
% whenever the anti-core exits. The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% Note: Windows users must port the shell script 'corevert' to get full 
% operationality.
%
% Usage: [acore_vert acrst]=sclv.AntiCoreVertices(method,tol)
% Define variables:
%  output:
%  core_vert  -- Matrix of anti-core vertices. Output is numeric or a string.
%  crst       -- The anti-core constraints.
%
%  input:
%  v          -- TuGame class object. 
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
%   07/15/2015        0.7             hme
%                

% Here we assume that the user has represented the game correctly.

v=clv.tuvalues;
N=clv.tusize;
n=clv.tuplayers;
gt=clv.tutype;

if nargin<2
  method='float';
  tol=10^9*eps;
 elseif nargin<3
  tol=10^9*eps; 
 else 
end

if clv.CddAntiCoreQ(tol)==0
  error('Anti-Core is empty!');
 else
end

if strcmp(method,'gmp')
elseif strcmp(method,'float')
 else
   error('Method not recognized. Permissible methods are: float or gmp!');
end

S=1:N;
it=0:-1:1-n;
PlyMat=rem(floor(S(:)*pow2(it)),2);
%
% Defining core constraints.
%
PlyMat(N+1,:)=-PlyMat(end,:);
v(N+1)=-v(N);

% Defining format
crst=zeros(N+1,n+1);
crst(:,2:n+1)=-PlyMat;
crst(:,1)=v';
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
str='corevert_aux';
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
!corevert -m gmp corevert_aux.ine
elseif strcmp(method,'float')
!corevert -m float corevert_aux.ine
else
!corevert -m gmp corevert_aux.ine
end


% Reading from file
fk=factorial(n);
[fid,~]=fopen(datafile,'r');
if fid>0
   if strcmp(method,'gmp')
      jj=1;
      while 1
         tline = fgetl(fid);
         if ~ischar(tline)
              break
         else
            core_vert(jj,:)=str2num(tline);
            jj=jj+1;
          end
      end
   else
      [core_vert count]=fscanf(fid,'%f',[n+1 fk]);
      if count==0
        error('No vertex computed! Use method gmp instead.');
       else
      end
   end
  fclose(fid);
 else
  error('Output file open failed');
end

% Reformating the output
if strcmp(method,'gmp')
  if isempty(core_vert)
   core_vert=tline;
   core_vert=str2num(core_vert);
  else
   core_vert(:,1)=[]; % deleting the first column vector of ones.
  end
% Converting string into real numbers in order to reformat output. 
%  core_vert=str2num(core_vert);
% Converting real numbers into a rational number approximation represented 
% as a string.
  core_vert=rats(core_vert); 
else
  core_vert=core_vert';
  core_vert(:,1)=[]; % deleting the first column vector of ones.
end
