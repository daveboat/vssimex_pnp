function [x,dx]=makemesh(N1,N2, N3, L1, L2, L3)
% makemesh,m
%
% The interval [0,L] is divided into 3 subintervals of length L1, L2, and
% L3 respectively.
% 
% The mesh x has N1+N2+N3 nodes, x(1) = 0, x(N1+N2+N3) = L.  The first
% subinterval has a uniform mesh of width L1/(N1-1), the second has width
% L2/N2, and the third has L3/N3.
%
% returns dx, the vector of mesh widths and x, the mesh vector.  Note that
%
% For a uniform mesh on [0,L], we take L1=L2=L3 and N1-1=N2=N3.  This will 
% result in 3*N2 intervals and 3*N2+1 meshpoints.  The first subinterval 
% [0,L1] has meshpoints at 0 and at L1: x(1) = 0 and x(N1) = L1, the second 
% subinterval (L1,L2] has meshpoint x(N1+1) is to the right of L1 and 
% x(N1+N2) = L1+L2, the third subinterval (L1+L2,L1+L2+L3] has meshpoint 
% x(N1+N2+1) to the right of L1+L2 and x(N1+N2+N3) = L1+L2+L3.


t1=linspace(1,N1,N1);
t2=linspace(N1+1,N1+N2,N2);
t3=linspace(N1+N2+1,N1+N2+N3,N3);

%piecewise linear mesh
x=zeros(N1+N2+N3,1);
x(1:N1) = L1*(t1-1)/(N1-1);
x(N1+1:N1+N2) = L2*(t2-N1)/N2 + L1;
x(N1+N2+1:N1+N2+N3) = L3*(t3-N1-N2)/N3+L1+L2;

%straight up linear, debug
% N=N1+N2+N3;
% t=linspace (1,N,N);
% x=t/N;

dx=diff(x);
end