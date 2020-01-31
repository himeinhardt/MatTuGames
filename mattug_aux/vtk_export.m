function vtk_export(v,filename,crit_val)
% VTK_EXPORT exportes the graphical raw data to VTK legacy format.
% Requires VTKWRITE from 
% http://www.mathworks.com/matlabcentral/fileexchange/47814-export-3d-data-to-paraview-in-vtk-legacy-file-format
%  
%   Usage: vtk_export(v,'Game02')
%          It exports the imputation set, the core, a strong epsilon core, the core core, 
%          the Weber set, the upper, and lower set as well as their faces. The vtk data
%          can be visualized, for instance, with ParaView.  
%
%  output:
%               -- Writes the graphical raw data to vtk-files.  
%
%  input:
%  v            -- A Tu-Game v of length 2^n-1.
%  filename     -- A string to specify the basename of the vtk-files.
%  crit_val     -- To determine the epsilon value of a strong core
%                  that should be drawn.
%
%
%


%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/22/2015        0.6             hme
%



if nargin <2
  s1='ExpGame';
  crit_val='';
else
  s1=filename;
  crit_val='';
end

% Exporting Imputation Set
[ip_vert,~,ip_vol,Pip]=CddImputationVertices(v);
Pip.minHRep;
ip_facet=Pip.getFacet;
nfi=numel(ip_facet);
for jj=1:nfi
   ipM=ip_facet(jj).incidenceMap;
   ifm=ipM.incVH;
   isf=size(ifm,2);
   bvi=[];for k=1:isf bvi=[bvi;ip_facet(jj).V(ifm(:,k),:)];end
   x=bvi(:,1);
   y=bvi(:,2);
   z=bvi(:,3);
   s3='ImpFacet%.03d.vtk';
   fn = strcat(s1,s3);
   fname=sprintf(fn,jj);
   vtkwrite(fname,'polydata','lines',x,y,z);
end

xi = Pip.V(:,1);
yi = Pip.V(:,2);
zi = Pip.V(:,3);
Ki = convhulln (Pip.V);
s2='Imputation.vtk';
fn = strcat(s1,s2);
vtkwrite(fn,'polydata','triangle', xi,yi,zi,Ki);


% Exporting Core
rg=range(ip_vert);
[~,idx]=min(rg);
[cr_vert,~,cr_vol,Pcr]=CddCoreVertices(v,idx);
Pcr.minHRep;
cr_facet=Pcr.getFacet;
nf=numel(cr_facet);
if nf>2
   for ii=1:nf
       iM=cr_facet(ii).incidenceMap;
       fm=iM.incVH;
       sf=size(fm,2);
       bv=[];for k=1:sf bv=[bv;cr_facet(ii).V(fm(:,k),:)];end
       x=bv(:,1);
       y=bv(:,2);
       z=bv(:,3);
       s4='CoreFacet%.03d.vtk';
       fn=strcat(s1,s4);
       fname=sprintf(fn,ii);
%   fname=sprintf('core_facet%.03d.vtk',ii);
       vtkwrite(fname,'polydata','lines',x,y,z);
   end
   try
     xc = Pcr.V(:,1);
     yc = Pcr.V(:,2);
     zc = Pcr.V(:,3);
     Kc = convhulln (Pcr.V);
     s5='Core.vtk';
     fn = strcat(s1,s5);
     vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);
   catch
     cM=Pcr.incidenceMap;
     fm=cM.incVH;
     sf=size(fm,2);
     bv=[];for k=1:sf bv=[bv;Pcr.V(fm(:,k),:)];end
     x=bv(:,1);
     y=bv(:,2);
     z=bv(:,3);
     s4='CoreFacet.vtk';
     fn=strcat(s1,s4);
     vtkwrite(fn,'polydata','lines',x,y,z); 
   end
elseif nf==2
       x=Pcr.V(:,1);
       y=Pcr.V(:,2);
       z=Pcr.V(:,3);
       s4='CoreLine.vtk';
       fn=strcat(s1,s4);
       vtkwrite(fn,'polydata','lines',x,y,z,'Precision',5);
else
      [x1,y1,z1] = sphere(75);
      crpt=Pcr.V;
      x0=crpt(1);
      y0=crpt(2);
      z0=crpt(3);
      x = x1*r + x0;
      y = y1*r + y0;
      z = z1*r + z0;
      x = [x(:);0];
      y = [y(:);0];
      z = [z(:);0];
      DTC = delaunayTriangulation(x,y,z);
      sL='CorePoint.vtk';
      fn=strcat(s1,sL);
      vtkwrite(fn,'polydata','tetrahedron',x,y,z,DTC.ConnectivityList);
end

% Radius for Point Solutions
tol=10^6*eps;
if abs(cr_vol)<tol
 r=ip_vol/50;
elseif abs(cr_vol)>tol
 r=cr_vol/ip_vol;
else
 r=1/2;
end

% Exporting Pre-Nucleolus, Pre-Kernel, and Shapley value
v_prk=PreKernel(v);
v_prn=CddPrenucl(v);
v_sh=ShapleyValue(v);
ResQ_v_pk=MatQ(v,v_prk,1);
rk=rank(ResQ_v_pk.matQ);
n=length(v_prk);
v_prk(:,idx)=[];
if rk<n
   P_v1=Polyhedron(ResQ_v_pk.matQ,ResQ_v_pk.vec_b);
   P_v1.minVRep;
   pk2=P_v1.V;
   pk2(:,idx)=[];
end
v_prn(:,idx)=[];
v_sh(:,idx)=[];
[x2,y2,z2] = sphere(50);
% Pre-Kernel;
if rk==n
  x0=v_prk(1);
  y0=v_prk(2);
  z0=v_prk(3);
  xK = x2*r + x0;                                              
  yK = y2*r + y0;                                              
  zK = z2*r + z0;
  x2K = [xK(:);0];
  y2K = [yK(:);0];
  z2K = [zK(:);0];
  DTK = delaunayTriangulation(x2K,y2K,z2K);
  sL='PreKr.vtk';
  fn=strcat(s1,sL);
  vtkwrite(fn,'polydata','tetrahedron',x2K,y2K,z2K,DTK.ConnectivityList);
else
  xkr=[v_prk(1);pk2(1)];
  ykr=[v_prk(2);pk2(2)];
  zkr=[v_prk(3);pk2(3)];
  sL='PreKrLine.vtk';
  fn=strcat(s1,sL);
  vtkwrite(fn,'polydata','lines',xkr,ykr,zkr,'Precision',5);
end
% Pre-Nucleolus
x0=v_prn(1);
y0=v_prn(2);
z0=v_prn(3);
xN = x2*r + x0;
yN = y2*r + y0; 
zN = z2*r + z0;
x2N = [xN(:);0];
y2N = [yN(:);0];
z2N = [zN(:);0];
DTN = delaunayTriangulation(x2N,y2N,z2N);
sL='PreNuc.vtk';
fn=strcat(s1,sL);
vtkwrite(fn,'polydata','tetrahedron',x2N,y2N,z2N,DTN.ConnectivityList);
% Shapley Value
x0=v_sh(1);
y0=v_sh(2);
z0=v_sh(3);
xS = x2*r + x0;
yS = y2*r + y0;
zS = z2*r + z0;
x2S = [xS(:);0];
y2S = [yS(:);0];
z2S = [zS(:);0];
DTS = delaunayTriangulation(x2S,y2S,z2S);
sL='ShapVal.vtk';
fn=strcat(s1,sL);
vtkwrite(fn,'polydata','tetrahedron',x2S,y2S,z2S,DTS.ConnectivityList);

% Exporting Strong Core
if isempty(crit_val)
   ctv1=critical_value1(v);
   ctv2=critical_value2(v);
   ctv3=critical_value_star(v);
   vc=[ctv1,ctv2,ctv3];
   ctv=max(vc);
   if ctv<=0
      ctv=2;
   end
else
   ctv=crit_val;
end

v_eps=streps_value(v,ctv);
[~,~,~,Pstc]=CddCoreVertices(v_eps,idx);
Pstc.minHRep;
stc_facet=Pstc.getFacet;
nf=numel(stc_facet);
if nf>2
   for ii=1:nf
       iM=stc_facet(ii).incidenceMap;
       fm=iM.incVH;
       sf=size(fm,2);
       bv=[];for k=1:sf bv=[bv;stc_facet(ii).V(fm(:,k),:)];end
       x=bv(:,1);
       y=bv(:,2);
       z=bv(:,3);
       s4='StrongCoreFacet%.03d.vtk';
       fn=strcat(s1,s4);
       fname=sprintf(fn,ii);
       vtkwrite(fname,'polydata','lines',x,y,z);
   end
   xc = Pstc.V(:,1);
   yc = Pstc.V(:,2);
   zc = Pstc.V(:,3);
   Kc = convhulln (Pstc.V);
   s6='StrongCore.vtk';
   fn = strcat(s1,s6);
   vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);
elseif nf==2
       x=Pstc.V(:,1);
       y=Pstc.V(:,2);
       z=Pstc.V(:,3);
       s4='StrongCoreLine.vtk';
       fn=strcat(s1,s4);
       vtkwrite(fn,'polydata','lines',x,y,z);
else
      [x1,y1,z1] = sphere(75);
      crpt=Pstc.V;
      x0=crpt(1);
      y0=crpt(2);
      z0=crpt(3);
      x = x1*r + x0;
      y = y1*r + y0;
      z = z1*r + z0;
      x = [x(:);0];
      y = [y(:);0];
      z = [z(:);0];
      DTC = delaunayTriangulation(x,y,z);
      sL='StrongCorePoint.vtk';
      fn=strcat(s1,sL);
      vtkwrite(fn,'polydata','tetrahedron',x,y,z,DTC.ConnectivityList);
end

% Exporting Least Core
LCV=CddLeastCoreVertices(v,idx);
LCV.P.minHRep;
lc_facet=LCV.P.getFacet;
nf=numel(lc_facet);
if nf>2
   for ii=1:nf
       iM=lc_facet(ii).incidenceMap;
       fm=iM.incVH;
       sf=size(fm,2);
       bv=[];for k=1:sf bv=[bv;lc_facet(ii).V(fm(:,k),:)];end
       x=bv(:,1);
       y=bv(:,2);
       z=bv(:,3);
       s4='LeastCoreFacet%.03d.vtk';
       fn=strcat(s1,s4);
       fname=sprintf(fn,ii);
       vtkwrite(fname,'polydata','lines',x,y,z);
   end
   xc = LCV.P.V(:,1);
   yc = LCV.P.V(:,2);
   zc = LCV.P.V(:,3);
   Kc = convhulln (LCV.P.V);
   s6='LeastCore.vtk';
   fn = strcat(s1,s6);
   vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);
elseif nf==2
       x=LCV.P.V(:,1);
       y=LCV.P.V(:,2);
       z=LCV.P.V(:,3);
       s4='LeastCoreLine.vtk';
       fn=strcat(s1,s4);
       vtkwrite(fn,'polydata','lines',x,y,z);
else
      [x1,y1,z1] = sphere(75);
      crpt=LCV.P.V;
      x0=crpt(1);
      y0=crpt(2);
      z0=crpt(3);
      x = x1*r + x0;
      y = y1*r + y0;
      z = z1*r + z0;
      x = [x(:);0];
      y = [y(:);0];
      z = [z(:);0];
      DTC = delaunayTriangulation(x,y,z);
      sL='LeastCorePoint.vtk';
      fn=strcat(s1,sL);
      vtkwrite(fn,'polydata','tetrahedron',x,y,z,DTC.ConnectivityList);
end


% Exporting Weber Set
[~,~,~,Pweb]=CddWeberSet(v,idx);
Pweb.minHRep;
wbs_facet=Pweb.getFacet;
nf=numel(wbs_facet);
if nf>2
   for ii=1:nf
       iM=wbs_facet(ii).incidenceMap;
       fm=iM.incVH;
       sf=size(fm,2);
       bv=[];for k=1:sf bv=[bv;wbs_facet(ii).V(fm(:,k),:)];end
       x=bv(:,1);
       y=bv(:,2);
       z=bv(:,3);
       s4='WeberSetFacet%.03d.vtk';
       fn=strcat(s1,s4);
       fname=sprintf(fn,ii);
       vtkwrite(fname,'polydata','lines',x,y,z);
   end
   try
     xc = Pweb.V(:,1);
     yc = Pweb.V(:,2);
     zc = Pweb.V(:,3);
     Kc = convhulln (Pweb.V);
     s7='WeberSet.vtk';
     fn = strcat(s1,s7);
     vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);
   catch
     cM=Pweb.incidenceMap;
     fm=cM.incVH;
     sf=size(fm,2);
     bv=[];for k=1:sf bv=[bv;Pweb.V(fm(:,k),:)];end
     x=bv(:,1);
     y=bv(:,2);
     z=bv(:,3);
     s4='WeberSetFacet.vtk';
     fn=strcat(s1,s4);
     vtkwrite(fn,'polydata','lines',x,y,z);
   end
elseif nf==2
       x=Pweb.V(:,1);
       y=Pweb.V(:,2);
       z=Pweb.V(:,3);
       s4='WeberSetLine.vtk';
       fn=strcat(s1,s4);
       vtkwrite(fn,'polydata','lines',x,y,z);
else
      [x1,y1,z1] = sphere(75);
      crpt=Pweb.V;
      x0=crpt(1);
      y0=crpt(2);
      z0=crpt(3);
      x = x1*r + x0;
      y = y1*r + y0;
      z = z1*r + z0;
      x = [x(:);0];
      y = [y(:);0];
      z = [z(:);0];
      DTC = delaunayTriangulation(x,y,z);
      sL='CorePoint.vtk';
      fn=strcat(s1,sL);
      vtkwrite(fn,'polydata','tetrahedron',x,y,z,DTC.ConnectivityList);
end

%Exporting Core Cover
[~,~,~,Pcc]=CddCoreCoverVertices(v,idx);
Pcc.minHRep;
cc_facet=Pcc.getFacet;
nf=numel(cc_facet);
if nf>2
   for ii=1:nf
       iM=cc_facet(ii).incidenceMap;
       fm=iM.incVH;
       sf=size(fm,2);
       bv=[];for k=1:sf bv=[bv;cc_facet(ii).V(fm(:,k),:)];end
       x=bv(:,1);
       y=bv(:,2);
       z=bv(:,3);
       s4='CoreCoverFacet%.03d.vtk';
       fn=strcat(s1,s4);
       fname=sprintf(fn,ii);
       vtkwrite(fname,'polydata','lines',x,y,z);
   end
   xc = Pcc.V(:,1);
   yc = Pcc.V(:,2);
   zc = Pcc.V(:,3);
   Kc = convhulln (Pcc.V);
   s8='CoreCover.vtk';
   fn = strcat(s1,s8);
   vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);
elseif nf==2
       x=Pcc.V(:,1);
       y=Pcc.V(:,2);
       z=Pcc.V(:,3);
       s4='CoreCoverLine.vtk';
       fn=strcat(s1,s4);
       vtkwrite(fn,'polydata','lines',x,y,z);
else
      [x1,y1,z1] = sphere(75);
      crpt=Pcc.V;
      x0=crpt(1);
      y0=crpt(2);
      z0=crpt(3);
      x = x1*r + x0;
      y = y1*r + y0;
      z = z1*r + z0;
      x = [x(:);0];
      y = [y(:);0];
      z = [z(:);0];
      DTC = delaunayTriangulation(x,y,z);
      sL='CoreCoverPoint.vtk';
      fn=strcat(s1,sL);
      vtkwrite(fn,'polydata','tetrahedron',x,y,z,DTC.ConnectivityList);
end

% Exporting Lower Set
[~,~,~,Plws]=CddLowerSetVertices(v,idx);
Plws.minHRep;
lws_facet=Plws.getFacet;
nf=numel(lws_facet);
for ii=1:nf
   iM=lws_facet(ii).incidenceMap;
   fm=iM.incVH;
   sf=size(fm,2);
   bv=[];for k=1:sf bv=[bv;lws_facet(ii).V(fm(:,k),:)];end
   x=bv(:,1);
   y=bv(:,2);
   z=bv(:,3);
   s4='LowerSetFacet%.03d.vtk';
   fn=strcat(s1,s4);
   fname=sprintf(fn,ii);
   vtkwrite(fname,'polydata','lines',x,y,z);
end

xc = Plws.V(:,1);
yc = Plws.V(:,2);
zc = Plws.V(:,3);
Kc = convhulln (Plws.V);
s9='LowerSet.vtk';
fn = strcat(s1,s9);
vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);

%Exporting Upper Set
[~,~,~,Pups]=CddUpperSetVertices(v,idx);
Pups.minHRep;
ups_facet=Pups.getFacet;
nf=numel(ups_facet);
for ii=1:nf
   iM=ups_facet(ii).incidenceMap;
   fm=iM.incVH;
   sf=size(fm,2);
   bv=[];for k=1:sf bv=[bv;ups_facet(ii).V(fm(:,k),:)];end
   x=bv(:,1);
   y=bv(:,2);
   z=bv(:,3);
   s4='UpperSetFacet%.03d.vtk';
   fn=strcat(s1,s4);
   fname=sprintf(fn,ii);
   vtkwrite(fname,'polydata','lines',x,y,z);
end

xc = Pups.V(:,1);
yc = Pups.V(:,2);
zc = Pups.V(:,3);
Kc = convhulln (Pups.V);
s10='UpperSet.vtk';
fn = strcat(s1,s10);
vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);

% Exporting Reasonable Set
[~,~,~,Preas]=CddReasonableSetVertices(v,idx);
Preas.minHRep;
reas_facet=Preas.getFacet;
nf=numel(reas_facet);
for ii=1:nf
   iM=reas_facet(ii).incidenceMap;
   fm=iM.incVH;
   sf=size(fm,2);
   bv=[];for k=1:sf bv=[bv;reas_facet(ii).V(fm(:,k),:)];end
   x=bv(:,1);
   y=bv(:,2);
   z=bv(:,3);
   s4='ReasSetFacet%.03d.vtk';
   fn=strcat(s1,s4);
   fname=sprintf(fn,ii);
   vtkwrite(fname,'polydata','lines',x,y,z);
end

xc = Preas.V(:,1);
yc = Preas.V(:,2);
zc = Preas.V(:,3);
Kc = convhulln (Preas.V);
s11='ReasSet.vtk';
fn = strcat(s1,s11);
vtkwrite(fn,'polydata','triangle', xc,yc,zc,Kc);
