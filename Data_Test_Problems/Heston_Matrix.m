function [Y,Dvm,Dsm,X,Dss,Iv,Ds,I,Dvv,Dv,Is,Dsmt,Dvmt,Xt,Yt,Dsst,Dvvt,Dst,Dvt,Ivt,u0,S,V,G,s,ds,v] = Heston_Matrix(m1_1,m2_1,K,T)
B = 0;           % lower barrier (must be <K)

%%%%%%%%%%%%%%%%%%%%%%%
% discretization data %
%%%%%%%%%%%%%%%%%%%%%%%

% s-direction
Smax = 8*K;
%m1=50;
%m1=200;
m1 = m1_1;
M1 = m1+1; 

% v-direction
Vmax = 5;
%m2=25;
%m2=50;
%m2=100;
m2 = m2_1;
M2 = m2+1;

%%%%%%%%%%%%%%%%
% spatial grid %
%%%%%%%%%%%%%%%%

% nonuniform s-mesh
c = K/10;
pointl = max(B,max(0.5,exp(-0.1*T))*K);
pointr = min(1.5,exp(0.1*T))*K;
xmin = asinh((B-pointl)/c);
xint = (pointr-pointl)/c;
xmax = xint + asinh((Smax-pointr)/c);
x  = linspace(xmin,xmax,M1);
x1 = x(x<=0);
x2 = x(0<x & x<xint);
x3 = x(x>=xint);
s1 = pointl + sinh(x1)*c;
s2 = pointl + x2*c;
s3 = pointr + sinh(x3-xint)*c;
s  = [s1 s2 s3];
s  = max(s,B);

% nonuniform v-mesh
v_in=0.0025;
d = Vmax/500;
W = asinh(Vmax/d);
w = linspace(v_in,W,M2);
v = d*sinh(w);
v = max(v,0);

% uniform mesh (not recommended)
% s = B+(Smax-B)*linspace(0,1,M1);
% v = Vmax*linspace(0,1,M2);

% quadratic mesh (as a paper)
% alpha=3/8; s=[]; v=[];
% for i=1:M1
%     if i==0
%     s(1)=(((i-1)/(alpha*M1)-1)*abs((i-1)/(alpha*M1)-1)+1)*K;
%     else
%         s(i)=(((i-1)/(alpha*M1)-1)*abs((i-1)/(alpha*M1)-1)+1)*K;
%     end
% end
% % s(2)=[];
% for j=1:M2
%     v(j)=(((j-1)/M2)^2)*Vmax;
% end
% v(2)=[];



i1 = find(s<2*K,1,'last');     % s(i1+1) >= 2K, s(i1) < 2K
j2 = find(v<1,1,'last');       % v(j2+1) >= 1,  v(j2) < 1

ii = 2:M1;
jj = 1:m2;
% s = s(ii);
% v= v(jj);
mi = max(size(ii));                                                         %GN =m1?
mj = max(size(jj));                                                         %GN =m2?

Xt = spdiags(s',0,M1,M1);
Yt = spdiags(v',0,M2,M2);
X  = Xt(ii,ii);
Y  = Yt(jj,jj);

[S,V] = ndgrid(s,v);                                                        %GN as mesh?


% initial function (payoff)
% an=1/2;
% u = 0.5*(S-K)+0.5*sqrt(((S-K).^2)+2*(an^2)*(sqrt(2)-1));
% initial function (payoff)
u = max(0,S-K);

% cell average near the strike                                               %GN regolarizzazione del dato iniziale 
% ind = find(s>K,1,'first');
% if abs(s(ind)-K) > abs(s(ind-1)-K)                                           %Why?
%     ind = ind-1;                   
% end
% sl = (s(ind-1)+s(ind))/2;
% sr = (s(ind)+s(ind+1))/2;
% u(ind,:) = 0.5*(sr-K)^2/(sr-sl);

u = u(ii,jj);
u = u(:); u0=u;

% Dirichlet at s=0 and v=Vmax
G = zeros(M1,M2); G(:,M2) = (s-B)';

% FD coefficients s-mesh
ds = (s(2:m1+1)-s(1:m1))';
hh = ds(1:m1-1)+ds(2:m1);
% convection central
betsl = -ds(2:m1)./(ds(1:m1-1).*hh);
betsp =  ds(1:m1-1)./(ds(2:m1).*hh);
betsm = -betsl-betsp;
% diffusion central
delsl = 2./(ds(1:m1-1).*hh);
delsp = 2./(ds(2:m1).*hh);
delsm = -delsl-delsp;

betsl = [betsl;0;0];
betsm = [0;betsm;0];
betsp = [0;0;betsp];
delsl = [delsl;0;0];
delsm = [0;delsm;0];
delsp = [0;0;delsp];

% FD coefficients v-mesh
dv = (v(2:m2+1)-v(1:m2))';
hh = dv(1:m2-1)+dv(2:m2);
% convection central
betvl = -dv(2:m2)./(dv(1:m2-1).*hh);
betvp = dv(1:m2-1)./(dv(2:m2).*hh);
betvm = -betvl-betvp;
% convection backward
alpvk = -betvl;                                                             %GN why?
alpvl = -hh./(dv(1:m2-1).*dv(2:m2));
alpvm = -alpvk-alpvl;
% diffusion central
delvl = 2./(dv(1:m2-1).*hh);
delvp = 2./(dv(2:m2).*hh);
delvm = -delvl-delvp;

alpvk = [alpvk;0;0];
alpvl = [0;alpvl;0];
alpvm = [0;0;alpvm];
betvl = [betvl;0;0];
betvm = [0;betvm;0];
betvp = [0;0;betvp];
delvl = [delvl;0;0];
delvm = [0;delvm;0];
delvp = [0;0;delvp];

% CONVECTION
% central
Dst = spdiags([betsl betsm betsp],[-1 0 1],M1,M1);
% Neumann at s=Smax
Dst(M1,:) = 0;
% central
Dvc = spdiags([betvl betvm betvp],[-1 0 1],M2,M2);
% backward
Dvb = spdiags([alpvk alpvl alpvm],[-2 -1 0],M2,M2);
% combination
Dvt = Dvb; Dvt(1:j2,:) = Dvc(1:j2,:);                                       %GN fino a v<1 apply central then backward
% forward at v=0
Dvt(1,:) = 0;
Dvt(1,1) = -(2*dv(1)+dv(2))/(dv(1)*(dv(1)+dv(2)));
Dvt(1,2) = (dv(1)+dv(2))/(dv(1)*dv(2));
Dvt(1,3) = -dv(1)/(dv(2)*(dv(1)+dv(2)));

% DIFFUSION
% central
Dsst = spdiags([delsl delsm delsp],[-1 0 1],M1,M1);
% Neumann at s=Smax
Dsst(M1,:) = 0;
Dsst(M1,M1-1) = 2/(ds(m1)^2);                                               %GN detla at neumann
Dsst(M1,M1)  = -2/(ds(m1)^2);                                               %GN delta at neumann
% central
Dvvt = spdiags([delvl delvm delvp],[-1 0 1],M2,M2);

% MIXED DERIVATIVE
% central
Dsmt = spdiags([betsl betsm betsp],[-1 0 1],M1,M1);
Dvmt = spdiags([betvl betvm betvp],[-1 0 1],M2,M2);
% Neumann at s=Smax
Dsmt(M1,:) = 0;

% restriction to actual grid points
Ds  = Dst(ii,ii);
Dv  = Dvt(jj,jj);
Dss = Dsst(ii,ii);
Dvv = Dvvt(jj,jj);
Dsm = Dsmt(ii,ii);
Dvm = Dvmt(jj,jj);

Is = speye(mi,mi);
Iv = speye(mj,mj);
Ivt = speye(M2,M2);

mij = mi*mj; I = speye(mij,mij);

% S = S(1:i1,1:j2);
% V = V(1:i1,1:j2);


end

