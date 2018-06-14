%testSimpleCiliaProblemV3.m

% tests analytic integrals for cilium. All line integrals analytic

clear all

epsilon=0.01;
domain='i';
blockSize=0.2;

m=5;
q=12;
a=3;

% stationary points
[xb,Xb]=GenerateSpherePoints(m,q,a);
vb=0*xb;

% cilia force points
Nc=10;
s=linspace(0,1,Nc+1)+1/2/Nc;s(end)=[];s=s(:);
xc=[0*s;0*s;s-a];
vc=[s;0*s;0*s];

Qc=20;
S=linspace(0,1,Qc+1)+1/2/Qc;S(end)=[];S=S(:);
Xc=[0*S;0*S;S-a];

%%
% Full nearest neighbour implementation

x=MergeVectorGrids(xb,xc);
v=MergeVectorGrids(vb,vc);
X=MergeVectorGrids(Xb,Xc);

[A,NN]=AssembleStokesletMatrix(x,X,x,epsilon,domain,blockSize);
f=A\v;

% report velocity at centre of sphere
xTest=[0*S;0*S;-a+2*a*S];
[A,NN]=AssembleStokesletMatrix(xTest,X,x,epsilon,domain,blockSize);

vTest=A*f;
[xT1,xT2,xT3]=ExtractComponents(xTest);
[vT1,vT2,vT3]=ExtractComponents(vTest);

figure(1);clf;hold on;
quiver3(xT1,xT2,xT3,vT1,vT2,vT3);
[x1,x2,x3]=ExtractComponents(x);
plot3(x1,x2,x3,'.');
view([-4 20]);

figure(3);clf;hold on;title('centreline velocity in x1 direction through cilium')
plot(xT3,vT1,'b');
xlabel('x3');
ylabel('v1');

%%
% All line integrals analytic
%
% note: matrix structure is slightly different from the above! all boundary points first, then
% all cilia points

h=1/2/Nc;

% bb interactions
[Abb,Nbb]=AssembleStokesletMatrix(xb,Xb,xb,epsilon,domain,blockSize);

% cb interaction
[Acb,Ncb]=AssembleStokesletMatrix(xc,Xb,xb,epsilon,domain,blockSize);

% bc interactions
R=kron(ones(1,Qc),[[0;0;1] [1;0;0] [0;1;0]]);
Abc=RegStokesletAnalyticIntegrals(xb,xc,h,R,epsilon);

% cc interactions
Acc=RegStokesletAnalyticIntegrals(xc,xc,h,R,epsilon);

% assemble
A=[Abb Abc ; Acb Acc];
v=[vb;vc];
f=A\v;

% report velocity at centre of sphere
xTest=[0*S;0*S;-a+2*a*S];
[Atb,Ntb]=AssembleStokesletMatrix(xTest,Xb,xb,epsilon,domain,blockSize);
Atc=RegStokesletAnalyticIntegrals(xTest,xc,h,R,epsilon);

vTest=[Atb Atc]*f;
[xT1,xT2,xT3]=ExtractComponents(xTest);
[vT1,vT2,vT3]=ExtractComponents(vTest);

figure(2);clf;hold on;
quiver3(xT1,xT2,xT3,vT1,vT2,vT3);
[x1,x2,x3]=ExtractComponents(x);
plot3(x1,x2,x3,'.');
view([-4 20]);

%%
figure(3);hold on;
plot(xT3,vT1,'r');
legend('full nearest-neighbour','analytic integral for cilium');
