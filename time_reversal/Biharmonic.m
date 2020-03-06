close all;
clear;
range=[320];
iter=0;
for m=range;
    iter=iter+1;
    h=1/(m+1);
    f=zeros(m+2);
    gx=zeros(m+2,2);
    gy=zeros(2,m+2);
    x=0:h:1;
    [X,Y]=meshgrid(x,x);
    
    u=exp(X.^2+Y.^2);
    du=4*(X.^2+Y.^2+1).*u;
    f=(X.^4+Y.^4+2*X.^2.*Y.^2+4*X.^2+4*Y.^2+2)*16.*u;
    df=64*(X.^6+3*X.^4.*Y.^2+9*X.^4+3*X.^2.*Y.^4+18*X.^2.*Y.^2+18*X.^2+Y.^6+9*Y.^4+18*Y.^2+6).*u;
    gx(:,2)=2*exp(1+(0:m+1).^2*h^2);
    gy(2,:)=2*exp(1+(0:m+1).^2*h^2);


%     u=exp(X).*sin(Y)+1;
%     du=zeros(m+2);
%     df=zeros(m+2);
%     gx(:,1)=-sin(0:h:1);
%     gx(:,2)=exp(1)*sin(0:h:1);
%     gy(1,:)=-exp(0:h:1);
%     gy(2,:)=exp(0:h:1)*cos(1);
    
    [U,U13,U17]=biharmonic(u,f,df(2:m+1,2:m+1),gx,gy,m,h);
    e(iter)=mean(mean(abs((U-u(2:m+1,2:m+1))./u(2:m+1,2:m+1))));
    e2(iter)=mean(mean(abs((U13-u(2:m+1,2:m+1))./u(2:m+1,2:m+1))));
    e3(iter)=mean(mean(abs((U17-u(2:m+1,2:m+1))./u(2:m+1,2:m+1))));
    dU=Laplacian(u,U,f,gx,gy,m,h);
    e1(iter)=mean(mean(abs((dU-du)./du)));
end
% figure(1)
% subplot(1,3,1)
% loglog(range,e)
% subplot(1,3,2)
% loglog(range,e2)
% subplot(1,3,3)
% loglog(range,e3)
% set(gcf,'position',[700,700,1600,400])
% figure(2)
% loglog(range,e1)
% beta=[ones(length(range),1),log(1./range)']\log(e')
% beta1=[ones(length(range),1),log(1./range)']\log(e1')
% figure(2)
% subplot(2,3,1)
% UU=u;
% UU(2:m+1,2:m+1)=U;
% UU13=u;
% UU(2:m+1,2:m+1)=U13;
% UU17=u;
% UU(2:m+1,2:m+1)=U17;
% imagesc(u)
% subplot(2,3,2)
% imagesc(UU)
% subplot(2,3,3)
% imagesc(u-UU)
% subplot(2,3,4)
% imagesc(du)
% subplot(2,3,5)
% imagesc(dU)
% subplot(2,3,6)
% imagesc(du-dU)
% set(gcf,'position',[500,500,2000,700])


%% Plot
function Plot(u1,g1)
m=length(u1)/4-1;
h=2/(m+1);
u=zeros(m+2);
gx=zeros(m+2,2);
gy=zeros(2,m+2);
f=zeros(m+2);
df=zeros(m+2);
u(1,1:m+1)=u1(1:m+1);
u(1:m+1,m+2)=u1(m+2:2*m+2);
u(m+2,m+2:-1:2)=u1(2*m+3:3*m+3);
u(m+2:-1:2,1)=u1(3*m+4:4*m+4);
gx(:,1)=g1(4*m+8:-1:3*m+7);
gx(:,2)=g1(m+3:2*m+4);
gy(1,:)=g1(1:m+2);
gy(2,:)=g1(3*m+6:-1:2*m+5);
U=biharmonic(u,f,df,gx,gy,m,h);
dU=Laplacian(u,U,f,gx,gy,m,h);
figure(1)
u(2:m+1,2:m+1)=U;
imagesc(u)
figure(2)
imagesc(dU)
end

%% 5-Stencil Laplacian
function du=Laplacian(u,tildeu,f,gx,gy,m,h)
i=-1:7;
b=[0;1;0;0;0;0;0;0;0];
A=fliplr(vander(i))';
a=linsolve(A,b)';
U=zeros(m+4);
U(2:m+3,2:m+3)=u;
U(3:m+2,3:m+2)=tildeu;
U(1,2:m+3)=-(h*gy(1,:)+a(2:9)*U(2:9,2:m+3))/a(1);
U(m+4,2:m+3)=-(h*gy(2,:)+a(9:-1:2)*U(m-4:m+3,2:m+3))/a(1);
U(2:m+3,1)=-(h*gx(:,1)+U(2:m+3,2:9)*a(2:9)')/a(1);
U(2:m+3,m+4)=-(h*gx(:,2)+U(2:m+3,m-4:m+3)*a(9:-1:2)')/a(1);
d1=diag(U);
d2=diag(fliplr(U));
U(1,1)=(-(gx(1,1)+gy(1,1))*h-a(2:9)*d1(2:9))/a(1);
U(m+4,m+4)=(-(gx(m+2,2)+gy(2,m+2))*h-a(9:-1:2)*d1(m-4:m+3))/a(1);
U(m+4,1)=(-(gx(m+2,1)+gy(2,1))*h-a(9:-1:2)*d2(m-4:m+3))/a(1);
U(1,m+4)=(-(gx(1,2)+gy(1,m+2))*h-a(2:9)*d2(2:9))/a(1);
e=ones(m+4,1);
T1=spdiags([4*e -20*e 4*e],-1:1,m+4,m+4);
T2=spdiags([e 4*e e],-1:1,m+4,m+4);
E1=spdiags(e,0,m+4,m+4);
E2=spdiags([e e],[-1 1],m+4,m+4);
A=kron(E1,T1)+kron(E2,T2);
d=reshape(A*reshape(U',[],1),[],m+4)'/(6*h^2);
du=d(2:m+3,2:m+3)-f*h^2/12;
end

%% Fourth order Biharmonic
function [U,U13,U17]=biharmonic(u,f,df,gx,gy,m,h)
    [A1,F1]=biharmonic13(u,f,gx,gy,m,h);
    [A2,F2]=biharmonic17(u,f,gx,gy,m,h);
    U=reshape((A1+A2)\reshape((F1+F2+h^6*df/2)',[],1),[],m)';
    U13=reshape(A1\reshape(F1',[],1),[],m)';
    U17=reshape(A2\reshape(F2',[],1),[],m)';
end

%% 17-Stencil Biharmonic
function [A,F]=biharmonic17(u,f,gx,gy,m,h)
i=-1:7;
b=[0;1;0;0;0;0;0;0;0];
A=fliplr(vander(i))';
a=linsolve(A,b)';
e=ones(m,1);
T1=spdiags([-2*e 16*e -2*e],-1:1,m,m);
T2=spdiags([e -4*e -2*e -4*e e],-2:2,m,m);
T3=spdiags([e e],[-1,1],m,m);
E1=spdiags(e,0,m,m);
E2=spdiags([e e],[-1 1],m,m);
E3=spdiags([e e],[-2 2],m,m);
A=kron(E1,T1)+kron(E2,T2)+kron(E3,T3);
D=zeros(m,m);
D(1,1:7)=a(3:9)/a(1);
D(m,m-6:m)=a(9:-1:3)/a(1);
A=A-kron(E2,D);
A(1:m,1:7*m)=A(1:m,1:7*m)-kron(a(3:9)/a(1),E2);
A(m^2-m+1:m^2,m^2-7*m+1:m^2)=A(m^2-m+1:m^2,m^2-7*m+1:m^2)-kron(a(9:-1:3)/a(1),E2);
F=f(2:m+1,2:m+1)*h^4*2;
F(1,2:m-1)=F(1,2:m-1)+(4+a(2)/a(1))*u(1,2:m-1)+(4+a(2)/a(1))*u(1,4:m+1)+2*u(1,3:m)-u(1,1:m-2)-u(1,5:m+2);
F(m,2:m-1)=F(m,2:m-1)+(4+a(2)/a(1))*u(m+2,2:m-1)+(4+a(2)/a(1))*u(m+2,4:m+1)+2*u(m+2,3:m)-u(m+2,1:m-2)-u(m+2,5:m+2);
F(2:m-1,1)=F(2:m-1,1)+(4+a(2)/a(1))*u(2:m-1,1)+(4+a(2)/a(1))*u(4:m+1,1)+2*u(3:m,1)-u(1:m-2,1)-u(5:m+2,1);
F(2:m-1,m)=F(2:m-1,m)+(4+a(2)/a(1))*u(2:m-1,m+2)+(4+a(2)/a(1))*u(4:m+1,m+2)+2*u(3:m,m+2)-u(1:m-2,m+2)-u(5:m+2,m+2);
F(1,1)=F(1,1)+4*u(1,1)+2*u(1,2)+2*u(2,1)+(4+a(2)/a(1))*u(1,3)+(4+a(2)/a(1))*u(3,1)-u(1,4)-u(4,1)+a(2:9)/a(1)*u(1,1:8)'+a(2:9)/a(1)*u(1:8,1);
F(1,m)=F(1,m)+4*u(1,m+2)+2*u(1,m+1)+2*u(2,m+2)+(4+a(2)/a(1))*u(1,m)+(4+a(2)/a(1))*u(3,m+2)-u(1,m-1)-u(4,m+2)+a(9:-1:2)/a(1)*u(1,m-5:m+2)'+a(2:9)/a(1)*u(1:8,m+2);
F(m,1)=F(m,1)+4*u(m+2,1)+2*u(m+2,2)+2*u(m+1,1)+(4+a(2)/a(1))*u(m+2,3)+(4+a(2)/a(1))*u(m,1)-u(m+2,4)-u(m-1,1)+a(2:9)/a(1)*u(m+2,1:8)'+a(9:-1:2)/a(1)*u(m-5:m+2,1);
F(m,m)=F(m,m)+4*u(m+2,m+2)+2*u(m+2,m+1)+2*u(m+1,m+2)+(4+a(2)/a(1))*u(m+2,m)+(4+a(2)/a(1))*u(m,m+2)-u(m+2,m-1)-u(m-1,m+2)+a(9:-1:2)/a(1)*u(m+2,m-5:m+2)'+a(9:-1:2)/a(1)*u(m-5:m+2,m+2);
F(2,:)=F(2,:)-u(1,1:m)-u(1,3:m+2);
F(:,2)=F(:,2)-u(1:m,1)-u(3:m+2,1);
F(m-1,:)=F(m-1,:)-u(m+2,1:m)-u(m+2,3:m+2);
F(:,m-1)=F(:,m-1)-u(1:m,m+2)-u(3:m+2,m+2);
F(2,1)=F(2,1)+u(1,1);
F(2,m)=F(2,m)+u(1,m+2);
F(m-1,1)=F(m-1,1)+u(m+2,1);
F(m-1,m)=F(m-1,m)+u(m+2,m+2);
F(1,2)=F(1,2)+u(1,1);
F(m,2)=F(m,2)+u(m+2,1);
F(1,m-1)=F(1,m-1)+u(1,m+2);
F(m,m-1)=F(m,m-1)+u(m+2,m+2);
F([1,m],:)=F([1,m],:)+(gy(:,1:m)+gy(:,3:m+2))*h/a(1);
F(:,[1,m])=F(:,[1,m])+(gx(1:m,:)+gx(3:m+2,:))*h/a(1);
end
%% 13-Stencil Biharmonic
function [A,F]=biharmonic13(u,f,gx,gy,m,h)
i=-1:7;
b=[0;1;0;0;0;0;0;0;0];
A=fliplr(vander(i))';
a=linsolve(A,b)';
e=ones(m,1);
T1=spdiags([e -8*e 20*e -8*e e],-2:2,m,m);
T1(1,1:7)=T1(1,1:7)-a(3:9)/a(1);
T1(m,m-6:m)=T1(m,m-6:m)-a(9:-1:3)/a(1);
T2=spdiags([2*e -8*e 2*e],-1:1,m,m);
T3=spdiags(e,0,m,m);
E1=spdiags(e,0,m,m);
E2=spdiags([e e],[-1 1],m,m);
E3=spdiags([e e],[-2 2],m,m);
A=kron(E1,T1)+kron(E2,T2)+kron(E3,T3);
A(1:m,1:7*m)=A(1:m,1:7*m)-kron(a(3:9)/a(1),E1);
A(m^2-m+1:m^2,m^2-7*m+1:m^2)=A(m^2-m+1:m^2,m^2-7*m+1:m^2)-kron(a(9:-1:3)/a(1),E1);
F=f(2:m+1,2:m+1)*h^4;
F(1,:)=F(1,:)-2*u(1,1:m)-2*u(1,3:m+2)+(8+a(2)/a(1))*u(1,2:m+1);
F(m,:)=F(m,:)-2*u(m+2,1:m)-2*u(m+2,3:m+2)+(8+a(2)/a(1))*u(m+2,2:m+1);
F(:,1)=F(:,1)-2*u(1:m,1)-2*u(3:m+2,1)+(8+a(2)/a(1))*u(2:m+1,1);
F(:,m)=F(:,m)-2*u(1:m,m+2)-2*u(3:m+2,m+2)+(8+a(2)/a(1))*u(2:m+1,m+2);
F([1,m],[1,m])=F([1,m],[1,m])+2*u([1,m+2],[1,m+2]);
F(2,:)=F(2,:)-u(1,2:m+1);
F(m-1,:)=F(m-1,:)-u(m+2,2:m+1);
F(:,2)=F(:,2)-u(2:m+1,1);
F(:,m-1)=F(:,m-1)-u(2:m+1,m+2);
F([1,m],:)=F([1,m],:)+gy(:,2:m+1)*h/a(1);
F(:,[1,m])=F(:,[1,m])+gx(2:m+1,:)*h/a(1);
end