
% results are perfect for  general field point and general centre, 
% but fibre aligned along x direction

% sign error in 3,1 ; 1,3 ; 2,3 and 3,2 entries when aligned along y direction

% some zv dependence, suggesting rotations aren't being done correctly

clear all

% arbitrary point
x0=[2;3.5;-8];
% arbitrary tangent angles
phi=(pi/180)*115;
th=(pi/180)*34;
% segment half-length
h=0.05;
% regularization parameter
epsilon=0.01;

% local frame
zv=[0;1;1];
tv=[cos(phi)*sin(th); sin(phi)*sin(th); cos(th)];
nv=cross(zv,tv);nv=nv/norm(nv);
bv=cross(tv,nv);bv=bv/norm(bv);

%%

%
% check frame
norm(tv)
norm(nv)
norm(bv)
norm(bv-cross(tv,nv))
norm(tv-cross(nv,bv))
norm(nv-cross(bv,tv))
%%
% numerical integral

x00=x0+[0.05;-0.2;0.5];

% line along which to integrate
xv=@(s) kron(x0,s*0+1)+kron(tv,s);

% quadrature points - midpoint rule
Q=501;
w=2*h/Q;
s=linspace(-h,h,Q+1)'+h/Q;s(end)=[];
xg=xv(s);
[xg1,xg2,xg3]=ExtractComponents(xg);
figure(1);clf;plot3(xg1,xg2,xg3,'.');axis equal;hold on;plot3(x00(1),x00(2),x00(3),'ro');

N=kron(eye(3),ones(Q,1));
S=RegStokeslet(x00,xg,epsilon);
ANumeric=S*N*w

%%
% analytic integral

R=[tv nv bv];
AAnalyticGeneral=RegStokesletAnalyticIntegrals(x00,x0,h,R,epsilon)

%%

