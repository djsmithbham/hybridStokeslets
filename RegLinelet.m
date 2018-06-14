function L=RegLinelet(h,epsilon)

% inner line integral of regularized stokeslet in local coordinates
% corrected from appendix B of D.J. Smith (2009) Proc. R. Soc. Lond. A 465, 3605-3626.

L=zeros(3);
alpha=1./sqrt(1+(epsilon/h)^2);
L(1,1)=4*atanh(alpha);
L(2,2)=2*(alpha+atanh(alpha));
L(3,3)=L(2,2);
L=1/8/pi*L;