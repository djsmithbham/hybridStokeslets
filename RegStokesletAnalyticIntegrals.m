function D=RegStokesletAnalyticIntegrals(x,X,h,R,epsilon)

% x is 3Nx1 field points
% X is 3Qx1 centres of intervals of integration
% h is half-length of integration
% R is a 3x3Q rotation matrix into the local coordinates of each interval
% of integration

x=x(:);
X=X(:);
[x1,x2,x3]=ExtractComponents(x);
[X1,X2,X3]=ExtractComponents(X);

N=length(x1);
Q=length(X1);

RM=kron(R,ones(N,1)); % RM is 3N x 3Q. Repeat for every field point.

RM11=repmat(R(1,      1:Q),N,1);
RM12=repmat(R(1,  Q+1:2*Q),N,1);
RM13=repmat(R(1,2*Q+1:3*Q),N,1);

RM21=repmat(R(2,      1:Q),N,1);
RM22=repmat(R(2,  Q+1:2*Q),N,1);
RM23=repmat(R(2,2*Q+1:3*Q),N,1);

RM31=repmat(R(3,      1:Q),N,1);
RM32=repmat(R(3,  Q+1:2*Q),N,1);
RM33=repmat(R(3,2*Q+1:3*Q),N,1);

xi1=x1*ones(1,Q)-ones(N,1)*X1';
xi2=x2*ones(1,Q)-ones(N,1)*X2';
xi3=x3*ones(1,Q)-ones(N,1)*X3';

% rotate into local coordinates.
xiL1=RM11.*xi1 + RM21.*xi2 + RM31.*xi3;
xiL2=RM12.*xi1 + RM22.*xi2 + RM32.*xi3;
xiL3=RM13.*xi1 + RM23.*xi2 + RM33.*xi3;

denom=xiL2.^2+xiL3.^2+epsilon^2;
idenom=1./denom;
repsM=sqrt((xiL1-h).^2+denom);
repsP=sqrt((xiL1+h).^2+denom);
irepsM=1./repsM;
irepsP=1./repsP;

I11=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*(epsilon^2.*idenom-1)...
    +2*log((h-xiL1+repsM)./(-h-xiL1+repsP));
I22=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*(epsilon^2+xiL2.^2).*idenom...
      +log((h-xiL1+repsM)./(-h-xiL1+repsP));
I33=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*(epsilon^2+xiL3.^2).*idenom...
      +log((h-xiL1+repsM)./(-h-xiL1+repsP));
I12=xiL2.*(irepsM-irepsP);
I21=I12;
I13=xiL3.*(irepsM-irepsP);
I31=I13;
I23=(-(xiL1-h).*irepsM+(xiL1+h).*irepsP).*xiL2.*xiL3.*idenom;
I32=I23;

% transform 2nd rank tensor - might be neater and safer with cell arrays
D11=(RM11.*I11+RM12.*I21+RM13.*I31).*RM11...
   +(RM11.*I12+RM12.*I22+RM13.*I32).*RM12...
   +(RM11.*I13+RM12.*I23+RM13.*I33).*RM13;

D12=(RM11.*I11+RM12.*I21+RM13.*I31).*RM21...
   +(RM11.*I12+RM12.*I22+RM13.*I32).*RM22...
   +(RM11.*I13+RM12.*I23+RM13.*I33).*RM23;

D13=(RM11.*I11+RM12.*I21+RM13.*I31).*RM31...
   +(RM11.*I12+RM12.*I22+RM13.*I32).*RM32...
   +(RM11.*I13+RM12.*I23+RM13.*I33).*RM33;

D21=(RM21.*I11+RM22.*I21+RM23.*I31).*RM11...
   +(RM21.*I12+RM22.*I22+RM23.*I32).*RM12...
   +(RM21.*I13+RM22.*I23+RM23.*I33).*RM13;

D22=(RM21.*I11+RM22.*I21+RM23.*I31).*RM21...
   +(RM21.*I12+RM22.*I22+RM23.*I32).*RM22...
   +(RM21.*I13+RM22.*I23+RM23.*I33).*RM23;

D23=(RM21.*I11+RM22.*I21+RM23.*I31).*RM31...
   +(RM21.*I12+RM22.*I22+RM23.*I32).*RM32...
   +(RM21.*I13+RM22.*I23+RM23.*I33).*RM33;

D31=(RM31.*I11+RM32.*I21+RM33.*I31).*RM11...
   +(RM31.*I12+RM32.*I22+RM33.*I32).*RM12...
   +(RM31.*I13+RM32.*I23+RM33.*I33).*RM13;

D32=(RM31.*I11+RM32.*I21+RM33.*I31).*RM21...
   +(RM31.*I12+RM32.*I22+RM33.*I32).*RM22...
   +(RM31.*I13+RM32.*I23+RM33.*I33).*RM23;

D33=(RM31.*I11+RM32.*I21+RM33.*I31).*RM31...
   +(RM31.*I12+RM32.*I22+RM33.*I32).*RM32...
   +(RM31.*I13+RM32.*I23+RM33.*I33).*RM33;

D=[D11 D12 D13; D21 D22 D23; D31 D32 D33]/8/pi;
