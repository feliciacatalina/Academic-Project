function [x,y,n,b,e] = RK4model(M,awal,akhir,xa,ya,na,ba,ea,b, a, r, K,L,d,beta,lambda,v,alpha,s,s1,s0,delta,Q,delta0)
% Fungsi Runge-Kutta 4 untuk solusi persamaan pada tugas pendahuluan
% a and b adalah batas kiri dan kanan 
% xa, ya nilai awal x(t0) = x0 dan y(t0)=y0
% M banyaknya subinterval

% sistem persamaan orde 1
Xdot = @(X,Y,N,B,E) (b-a*r*(N/K))*N - (d+(1-a)*r*N/K)*X - (beta*Y+lambda*beta)*X + v*Y;
Ydot = @(X,Y,N,B,E) (beta*Y+lambda*B)*X - (v+alpha+d+(1-a)*r*N/K)*Y;
Ndot = @(X,Y,N,B,E) r*(1-N/K)*N - alpha*Y;
Bdot = @(X,Y,N,B,E) s*B*(1-B/L)+s1*Y-s0*B+delta*B*E;
Edot = @(X,Y,N,B,E) Q*N-delta0*E;

h=(akhir-awal)/M; %step size
x=zeros(1,M+1);
y=zeros(1,M+1);
n=zeros(1,M+1);
b=zeros(1,M+1);
e=zeros(1,M+1);

x(1) = xa; %X(0) karena di matlab indexing dari 1
y(1) = ya; %Y(0) karena di matlab indexing dari 1
n(1) = na;
b(1) = ba;
e(1) = ea;

for j=1:M
    % mencari xdot1,ydot1
    xdot1=h*Xdot(x(j),y(j),n(j),b(j),e(j));
    ydot1=h*Ydot(x(j),y(j),n(j),b(j),e(j));
    ndot1=h*Ndot(x(j),y(j),n(j),b(j),e(j));
    bdot1=h*Bdot(x(j),y(j),n(j),b(j),e(j));
    edot1=h*Edot(x(j),y(j),n(j),b(j),e(j));

    %mencari xdot2,ydot2
    xdot2=h*Xdot(x(j)+xdot1/2,y(j)+ydot1/2,n(j)+ndot1/2,b(j)+bdot1/2,e(j)+edot1/2);
    ydot2=h*Ydot(x(j)+xdot1/2,y(j)+ydot1/2,n(j)+ndot1/2,b(j)+bdot1/2,e(j)+edot1/2);
    ndot2=h*Ndot(x(j)+xdot1/2,y(j)+ydot1/2,n(j)+ndot1/2,b(j)+bdot1/2,e(j)+edot1/2);
    bdot2=h*Bdot(x(j)+xdot1/2,y(j)+ydot1/2,n(j)+ndot1/2,b(j)+bdot1/2,e(j)+edot1/2);
    edot2=h*Edot(x(j)+xdot1/2,y(j)+ydot1/2,n(j)+ndot1/2,b(j)+bdot1/2,e(j)+edot1/2);

    %mencari xdot3,ydot3
    xdot3=h*Xdot(x(j)+xdot2/2,y(j)+ydot2/2,n(j)+ndot2/2,b(j)+bdot2/2,e(j)+edot2/2);
    ydot3=h*Ydot(x(j)+xdot2/2,y(j)+ydot2/2,n(j)+ndot2/2,b(j)+bdot2/2,e(j)+edot2/2);
    ndot3=h*Ndot(x(j)+xdot2/2,y(j)+ydot2/2,n(j)+ndot2/2,b(j)+bdot2/2,e(j)+edot2/2);
    bdot3=h*Bdot(x(j)+xdot2/2,y(j)+ydot2/2,n(j)+ndot2/2,b(j)+bdot2/2,e(j)+edot2/2);
    edot3=h*Edot(x(j)+xdot2/2,y(j)+ydot2/2,n(j)+ndot2/2,b(j)+bdot2/2,e(j)+edot2/2);
    %mencari xdot4,ydot4
    xdot4=h*Xdot(x(j)+xdot3,y(j)+ydot3,n(j)+ndot3,b(j)+bdot3,e(j)+edot3);
    ydot4=h*Ydot(x(j)+xdot3,y(j)+ydot3,n(j)+ndot3,b(j)+bdot3,e(j)+edot3);
    ndot4=h*Ndot(x(j)+xdot3,y(j)+ydot3,n(j)+ndot3,b(j)+bdot3,e(j)+edot3);
    bdot4=h*Bdot(x(j)+xdot3,y(j)+ydot3,n(j)+ndot3,b(j)+bdot3,e(j)+edot3);
    edot4=h*Edot(x(j)+xdot3,y(j)+ydot3,n(j)+ndot3,b(j)+bdot3,e(j)+edot3);

    %memasukan hasil x dan y dalam sebuah vektor
    x(j+1)=x(j)+(xdot1+2*xdot2+2*xdot3+xdot4)/6;
    y(j+1)=y(j)+(ydot1+2*ydot2+2*ydot3+ydot4)/6;
    n(j+1)=n(j)+(ndot1+2*ndot2+2*ndot3+ndot4)/6;
    b(j+1)=b(j)+(bdot1+2*bdot2+2*bdot3+bdot4)/6;
    e(j+1)=e(j)+(edot1+2*edot2+2*edot3+edot4)/6;
end
end