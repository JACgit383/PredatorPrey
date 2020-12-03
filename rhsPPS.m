function F=rhsPPS(t,u)
%use this
%u(1)=x
%u(2)=y

global a;
global r;
global alpha;
global c;
global gamma;

%our 2 differentials
F=zeros(length(u),1);
F(1)=(a.*u(1)-r.*(u(1).^2))-alpha.*u(1).*u(2);
F(2)=-c.*u(2)+gamma.*u(1).*u(2);
