function AntiCorePlot(varargin)
% ANTICOREPLOT plots at most a two (three) dimensional anti-core in case of a 
% three (four) person game v whenever the anti-core exits. 
% The cdd-library by Komei Fukuda is needed.
% http://www.cs.mcgill.ca/~fukuda/download/cdd
% Note: Windows users must port the shell script 'corevert' to get full 
% operationality.
%
% Usage: AntiCorePlot(varargin)
% Define variables:
%  output:
%             -- A plot of the anti-core of game v.
%
%  input:     At most five input arguments are admissible. Any order of
%             the input arguments below is allowed. Nevertheless, at
%             least a game v must be specified. 
%
%
%  v          -- A Tu-Game v of length 2^n-1.
%  method     -- A string to call a method from the cdd-library.
%                Permissible methods are: 
%                'float' that is, results are given by real numbers.
%                Default is 'float'.
%                'gmp' that is, results are given by rational numbers.
%                Choose this method whenever the result with 'float'
%                is not as expected. This method needs more time to complete. 
%  imp_set    -- An integer to draw the anti-core in connection with the imputation
%                set. Permissible values are 1 (true) or 0 (false).
%                Default is 1 (true). This option can not be used for
%                a game having less than 4 players.
%  add_sol    -- A string to invoke additional solutions into the plot.
%                Permissible solutions are:
%                'none', this is the default value.
%                'prk', a pre-kernel element will be incorporated.
%                'prn', the pre-nucleolus will be incorporated.
%                'shap', the Shapley value will be incorporated.
%                'all', all three solutions above will be incorporated.
%  tol        -- Positive tolerance value. Its default value is set to 10^9*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   11/24/2012        0.3             hme
%                



narginchk(1,5);

v='';
method='';
add_sol='';
imp_set='';
tol='';

for i=1:nargin

  if ischar(varargin{i})
    if strcmp(varargin{i},'gmp')
       method='gmp';
    elseif strcmp(varargin{i},'float')
       method='float';
    elseif strcmp(varargin{i},'none')
       add_sol='none';
    elseif strcmp(varargin{i},'prk')
       add_sol='prk';
    elseif strcmp(varargin{i},'prn')
       add_sol='prn';
    elseif strcmp(varargin{i},'shap')
       add_sol='shap';
    elseif strcmp(varargin{i},'all')
       add_sol='all';
    else
      error(['Input argument at position ' int2str(i) ' not recognized']);
    end
  else 
    if varargin{i}==1
      imp_set=varargin{i};     
    elseif varargin{i}==0
      imp_set=varargin{i};
    elseif varargin{i}>0 & varargin{i}<1
      tol=varargin{i};     
    elseif length(varargin{i})==7
      v=varargin{i};
    elseif length(varargin{i})==15
      v=varargin{i};
    else
     N=15;
     if varargin{i}>1
      if length(varargin{i})<N
        error('Game has not the correct size!'); 
      elseif length(varargin{i})>N
        error('Game has dimension larger than four!');
      else
        error(['Input argument at position ' int2str(i) ' not recognized']);
      end
     else
       error(['Input argument at position ' int2str(i) ' not recognized']);
     end
    end
  end
end


if isempty(v)
   error('At least the game must be given!');
elseif isempty(method)
  method='gmp';
elseif isempty(add_sol)
  add_sol='none';
elseif isempty(imp_set)
  imp_set=1;
elseif isempty(tol)
  tol=10^9*eps;
else 
end


N=length(v);
[~, n]=log2(N);


if CddAntiCoreQ(v,tol)==0
  error('Anti-Core is empty!');
 else
end


%zov=ZeroOne_Normalization(v);
zov=v;
if zov==0
  warning('CrPlt:noN','No normalization possible! Using original game data for plotting.');
  zov=v;
 else 
end
v_crv=AntiCoreVertices(zov,method);
sz=size(v_crv);

if strcmp(method,'gmp')
  v_crv=str2num(v_crv);
 else
end
if isempty(v_crv)
  error('Cannot convert string to numbers! Use method float.');
elseif sz(1)==1
  error('Anti-Core is a unique point only. No plot possible!');
else
end

if n==4
%  y=iqr(v_crv);
  y=range(v_crv);
  [~, ind]=min(y);
  X=v_crv;
  X(:,ind)=[]; % Deleting the column with smallest range.
  [rws,~]=size(X);

  if rws>3
    C=convhulln(X);
   elseif rws==3
    cr_x=sum(X)/3;
    X(end+1,:)=cr_x;
    C=convhulln(X);
    iX(end+1,:)=sum(iX)/3;
    iC=convhulln(iX); 
   else
  end 

  impv=AntiImputationVertices(zov);
  iX=impv;
  iX(:,ind)=[];
  iC=convhulln(iX); 

elseif n==3
  [rws,~]=size(v_crv);
  [X1 X2]=toSimplex(v_crv);

  if rws>2
   C=convhull(X1,X2); 
   tri=delaunay(X1,X2);
  end

  impv=AntiImputationVertices(zov);
  [isr,~]=size(impv);

  if isr<=2
   error('Imputation set is a line segment. No plot possible!');
  end

  [imX1 imX2]=toSimplex(impv);
  im_tri=delaunay(imX1,imX2);
else
  X=v_crv;
  C=convhulln(X);
end

clf;
if n==4
 figure(1); hold on
 if rws>2

  if imp_set==1
   h=trimesh(C,X(:,1),X(:,2),X(:,3));

   if strcmp(add_sol,'none')
     set(h,'FaceColor',[0 .5 1],'EdgeColor','k','LineWidth',.7);
   else
    set(h,'EdgeColor','k','LineWidth',.7,'FaceColor',[0 .5 1],'FaceAlpha',0.5);
   end

   trimesh(iC,iX(:,1),iX(:,2),iX(:,3),'FaceAlpha',0);
   text(iX(1,1)+.02,iX(1,2),iX(1,3),'Player 1');
   text(iX(2,1)-.15,iX(2,2),iX(2,3)+.09,'Player 2');
   text(iX(3,1),iX(3,2),iX(3,3)+.04,'Player 3');
   text(iX(4,1),iX(4,2)-.06,iX(4,3),'Player 4');
  else
   h=trimesh(C,X(:,1),X(:,2),X(:,3));

%   if strcmp(add_sol,'none')
%     set(h,'FaceColor',[0 .5 1],'EdgeColor','k','LineWidth',.7);
%   else
%    set(h,'EdgeColor','k','LineWidth',.7,'FaceColor',[0 .3 1],'FaceAlpha',0.5);
%   end

%  d=[ 1 2 3 1];
%  for i=1:size(C,1)
%  j=C(i,d);
%  h(i)=patch(X(j,1),X(j,2),X(j,3),i,'FaceAlpha',0.8);
%  end

   if strcmp(add_sol,'none')   
     set(h,'AmbientStrength',.2,'FaceColor',[1 0 0]);
     set(h,'FaceLighting','phong','EdgeColor','k','LineWidth',1.7);
   else
     set(h,'EdgeColor','k','LineWidth',1.7,'FaceColor',[0 .3 1],'FaceAlpha',0.5);
     set(h,'FaceLighting','phong'); %,'EdgeColor','k','LineWidth',1.7);
   end

  end




     if strcmp(add_sol,'none')
      elseif strcmp(add_sol,'prk')
         v_prk=PreKernel(zov);
         v_prk(:,ind)=[];
         h2=plot3(v_prk(1),v_prk(2),v_prk(3));
         set(h2,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
       elseif strcmp(add_sol,'prn')
         v_prn=PreNucl(zov);
         v_prn(:,ind)=[];
         h2=plot3(v_prn(1),v_prn(2),v_prn(3));
         set(h2,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
       elseif strcmp(add_sol,'shap')
         v_sh=ShapleyValue(zov);
         v_sh(:,ind)=[];
         h2=plot3(v_sh(1),v_sh(2),v_sh(3));
         set(h2,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
       elseif strcmp(add_sol,'all')
         v_prk=PreKernel(zov);
         v_prk(:,ind)=[];
         h2=plot3(v_prk(1),v_prk(2),v_prk(3));
         set(h2,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
         v_prn=PreNucl(zov);
         v_prn(:,ind)=[];
         h3=plot3(v_prn(1),v_prn(2),v_prn(3));
         set(h3,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
         v_sh=ShapleyValue(zov);
         v_sh(:,ind)=[];
         h4=plot3(v_sh(1),v_sh(2),v_sh(3));
         set(h4,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
       else
     end
   elseif rws==2
     h=line(X(:,1),X(:,2),X(:,3));
     if imp_set==1
       trimesh(iC,iX(:,1),iX(:,2),iX(:,3),'FaceAlpha',0);
       set(h,'LineWidth',3,'Color',[0 .5 1]);
       text(iX(1,1)+.02,iX(1,2),iX(1,3),'Player 1');
       text(iX(2,1)-.15,iX(2,2),iX(2,3)+.09,'Player 2');
       text(iX(3,1),iX(3,2),iX(3,3)+.04,'Player 3');
       text(iX(4,1),iX(4,2)-.06,iX(4,3),'Player 4');
     else 
       set(h,'LineWidth',1.5,'Color',[0 .5 1]);
     end

     if strcmp(add_sol,'none')
      elseif strcmp(add_sol,'prk')
         v_prk=PreKernel(zov);
         v_prk(:,ind)=[];
         h2=plot3(v_prk(1),v_prk(2),v_prk(3));
         set(h2,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
       elseif strcmp(add_sol,'prn')
         v_prn=PreNucl(zov);
         v_prn(:,ind)=[];
         h2=plot3(v_prn(1),v_prn(2),v_prn(3));
         set(h2,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
       elseif strcmp(add_sol,'shap')
         v_sh=ShapleyValue(zov);
         v_sh(:,ind)=[];
         h2=plot3(v_sh(1),v_sh(2),v_sh(3));
         set(h2,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
       elseif strcmp(add_sol,'all')
         v_prk=PreKernel(zov);
         v_prk(:,ind)=[];
         h2=plot3(v_prk(1),v_prk(2),v_prk(3));
         set(h2,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
         v_prn=PreNucl(zov);
         v_prn(:,ind)=[];
         h3=plot3(v_prn(1),v_prn(2),v_prn(3));
         set(h3,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
         v_sh=ShapleyValue(zov);
         v_sh(:,ind)=[];
         h4=plot3(v_sh(1),v_sh(2),v_sh(3));
         set(h4,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
       else
     end
   else
      trimesh(iC,iX(:,1),iX(:,2),iX(:,3),'FaceAlpha',0);

      h=plot3(X(1),X(2),X(3));
      set(h,'Marker','p','MarkerSize',5,'MarkerFaceColor',[0 .5 1]);

      if strcmp(add_sol,'none')
       elseif strcmp(add_sol,'prk')
         v_prk=PreKernel(zov);
         v_prk(:,ind)=[];
         h2=plot3(v_prk(1),v_prk(2),v_prk(3));
         set(h2,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
       elseif strcmp(add_sol,'prn')
         v_prn=PreNucl(zov);
         v_prn(:,ind)=[];
         h2=plot3(v_prn(1),v_prn(2),v_prn(3));
         set(h2,'Marker','^','MarkerSize',8,'MarkerFaceColor','c');
       elseif strcmp(add_sol,'shap')
         v_sh=ShapleyValue(zov);
         v_sh(:,ind)=[];
         h2=plot3(v_sh(1),v_sh(2),v_sh(3));
         set(h2,'Marker','o','MarkerSize',6,'MarkerFaceColor','y');
       elseif strcmp(add_sol,'all')
         v_prk=PreKernel(zov);
         v_prk(:,ind)=[];
         h2=plot3(v_prk(1),v_prk(2),v_prk(3));
         set(h2,'Marker','s','MarkerSize',6,'MarkerFaceColor','r');
         v_prn=PreNucl(zov);
         v_prn(:,ind)=[];
         h3=plot3(v_prn(1),v_prn(2),v_prn(3));
         set(h3,'Marker','^','MarkerSize',8,'MarkerFaceColor','w');
         v_sh=ShapleyValue(zov);
         v_sh(:,ind)=[];
         h4=plot3(v_sh(1),v_sh(2),v_sh(3));
         set(h4,'Marker','o','MarkerSize',6,'MarkerFaceColor','m');
	else
       end
  end
  hold off
if imp_set==1
%   view(-37.875, 31.5475), axis tight, axis off
   view(116, 30), axis tight, axis off
 else
    view(-7.5,10), axis equal, axis off
end
%    view(-7.5,10), axis equal, axis off
%    view(-17.5,20), axis equal, axis off
%    view(-37.5,30), axis equal, axis off
%  view(3)
%  box off
  if imp_set==1
  else  camorbit(90,-5);
  end
  camlight;
 lighting phong
 material shiny
% shading interp;

  title('The anti-core of the game')
elseif n==3
  figure(1); hold on

  if rws>2
   h1=patch('Faces',im_tri,'Vertices',[imX1 imX2]);
   set(h1,'FaceColor',[1 1 1],'EdgeColor',[0.5 0.5 0.5]);
   text(imX1(1)-.25,imX2(1),'Player 1');
   text(imX1(2)+.05,imX2(2),'Player 2');
   text(imX1(3)+.05,imX2(3),'Player 3');
   h=patch('Faces',tri,'Vertices',[X1 X2]);
    if strcmp(add_sol,'none')
      set(h,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
     else
     set(h,'FaceColor','c','EdgeColor','c');
    end

     if strcmp(add_sol,'none')
      elseif strcmp(add_sol,'prk')
         v_prk=PreKernel(zov);
         [x1 x2]=toSimplex(v_prk);
         h2=plot(x1,x2);
         set(h2,'Marker','s','MarkerSize',5,'MarkerFaceColor','r');
       elseif strcmp(add_sol,'prn')
         v_prn=PreNucl(zov);
         [y1 y2]=toSimplex(v_prn);
         h2=plot(y1,y2);
         set(h2,'Marker','^','MarkerSize',5,'MarkerFaceColor','b');
       elseif strcmp(add_sol,'shap')
         v_sh=ShapleyValue(zov);
         [z1 z2]=toSimplex(v_sh);
         h2=plot(z1,z2);
         set(h2,'Marker','o','MarkerSize',5,'MarkerFaceColor','y');
       elseif strcmp(add_sol,'all')
         v_prk=PreKernel(zov);
         [x1 x2]=toSimplex(v_prk);
         h2=plot(x1,x2);
         set(h2,'Marker','s','MarkerSize',5,'MarkerFaceColor','r');
         v_prn=PreNucl(zov);
         [y1 y2]=toSimplex(v_prn);
         h3=plot(y1,y2);
         set(h3,'Marker','^','MarkerSize',5,'MarkerFaceColor','c');
         v_sh=ShapleyValue(zov);
         [z1 z2]=toSimplex(v_sh);
         h4=plot(z1,z2);
         set(h4,'Marker','o','MarkerSize',5,'MarkerFaceColor','y');
       else
      end




  elseif rws==2
    h1=patch('Faces',im_tri,'Vertices',[imX1 imX2]);
    set(h1,'FaceColor',[1 1 1],'EdgeColor',[0.5 0.5 0.5]);
    text(imX1(1)-.25,imX2(1),'Player 1');
    text(imX1(2)+.05,imX2(2),'Player 2');
    text(imX1(3)+.05,imX2(3),'Player 3');
    h=line(X1,X2);
    set(h,'LineWidth',3,'Color',[1 0.5 0]);


     if strcmp(add_sol,'none')
      elseif strcmp(add_sol,'prk')
         v_prk=PreKernel(zov);
         [x1 x2]=toSimplex(v_prk);
         h2=plot(x1,x2);
         set(h2,'Marker','s','MarkerSize',5,'MarkerFaceColor','r');
       elseif strcmp(add_sol,'prn')
         v_prn=PreNucl(zov);
         [y1 y2]=toSimplex(v_prn);
         h2=plot(y1,y2);
         set(h2,'Marker','^','MarkerSize',5,'MarkerFaceColor','b');
       elseif strcmp(add_sol,'shap')
         v_sh=ShapleyValue(zov);
         [z1 z2]=toSimplex(v_sh);
         h2=plot(z1,z2);
         set(h2,'Marker','o','MarkerSize',5,'MarkerFaceColor','y');
       elseif strcmp(add_sol,'all')
         v_prk=PreKernel(zov);
         [x1 x2]=toSimplex(v_prk);
         h2=plot(x1,x2);
         set(h2,'Marker','s','MarkerSize',5,'MarkerFaceColor','r');
         v_prn=PreNucl(zov);
         [y1 y2]=toSimplex(v_prn);
         h3=plot(y1,y2);
         set(h3,'Marker','^','MarkerSize',5,'MarkerFaceColor','c');
         v_sh=ShapleyValue(zov);
         [z1 z2]=toSimplex(v_sh);
         h4=plot(z1,z2);
         set(h4,'Marker','o','MarkerSize',5,'MarkerFaceColor','y');
       else
      end
   else
    h=plot(X1,X2);
    set(h,'Marker','d','MarkerSize',25,'Color',[1 0.5 0]);
    h1=patch('Faces',im_tri,'Vertices',[imX1 imX2]);
    set(h1,'FaceColor',[1 1 1],'EdgeColor',[0.5 0.5 0.5]);
     
    if strcmp(add_sol,'none')
      elseif strcmp(add_sol,'prk')
         v_prk=PreKernel(zov);
         [x1 x2]=toSimplex(v_prk);
         h2=plot(x1,x2);
         set(h2,'Marker','s','MarkerSize',5,'MarkerFaceColor','r');
       elseif strcmp(add_sol,'prn')
         v_prn=PreNucl(zov);
         [y1 y2]=toSimplex(v_prn);
         h2=plot(y1,y2);
         set(h2,'Marker','^','MarkerSize',5,'MarkerFaceColor','b');
       elseif strcmp(add_sol,'shap')
         v_sh=ShapleyValue(zov);
         [z1 z2]=toSimplex(v_sh);
         h2=plot(z1,z2);
         set(h2,'Marker','o','MarkerSize',5,'MarkerFaceColor','y');
       elseif strcmp(add_sol,'all')
         v_prk=PreKernel(zov);
         [x1 x2]=toSimplex(v_prk);
         h2=plot(x1,x2);
         set(h2,'Marker','s','MarkerSize',5,'MarkerFaceColor','r');
         v_prn=PreNucl(zov);
         [y1 y2]=toSimplex(v_prn);
         h3=plot(y1,y2);
         set(h3,'Marker','^','MarkerSize',5,'MarkerFaceColor','c');
         v_sh=ShapleyValue(zov);
         [z1 z2]=toSimplex(v_sh);
         h4=plot(z1,z2);
         set(h4,'Marker','o','MarkerSize',5,'MarkerFaceColor','y');
       else
    end
  end
  hold off
  view(2)
  axis([-1 1 -1 1]), axis off
  title('The anti-core of the game');
 else
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X1 X2]=toSimplex(X);

[sr sc]=size(X);

if sr>1
 X1=[(X(:,2)-X(:,1))*sqrt(3)/2];
 X2=[ X(:,3)-(X(:,2)+X(:,1))/2];
else
 X1=[(X(2)-X(1))*sqrt(3)/2];
 X2=[ X(3)-(X(2)+X(1))/2];
end

