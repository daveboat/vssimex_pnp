function [ phi ] = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp,cm,z_cp,z_cm, phi_x_right,v )
% poisson1d.m 
%
% 1d poisson solver
%
% solves epsilon^2 d^2 phi/dx^2 = -1/2 sum_i(c_i)
%
% there are internal boundary conditions at x=L1 and x=L1+L2, specifically 
%
%     epsilon_1^2 phi_x_L(L1) = epsilon_2^2 phi_x_R(L1)
%     epsilon_2^2 phi_x_L(L1+L2) = epsilon_3^2 phi_x_R(L1+L2)
%
% where phi_x_L is the derivative approximated from the left and phi_x_R
% is the derivative approximate from the right.
%
% The boundary condition at x=0 is 
%
%        epsilon_s_left^2 (phi(0) - 0)/lambda_s = epsilon_1^2 phi_x(0)
%
% and the boundary condition at x=L=L1+L2+L3 is one of two options:
% current boundary conditions
%
%   epsilon_3^2 phi_x(L) = epsilon_s_right^2 phi_x_right
%
% where phi_x_right is found from solving the ODE 
%
% or voltage boundary conditions
%
% epsilon_3^2 phi_x(L) = epsilon_s_right^2 (v - phi(L))/(lambda_s)
%
% where v(t) is imposed (given by the voltage.m script)
%
% takes geometry vectors and constants, some polarization parameters, the
% cp and cm vectors, and phi_x_right and v for boundary conditions (but only
% uses one of the two depending on the 'bc' flag)

N=N1+N2+N3;

A = zeros(N,N);
b = zeros(N,1);

% the boundary condition at x=0:
%        epsilon_s_left^2 (phi(0) - 0)/lambda_s = epsilon_1^2 phi_x(0)
%
A(1,1) = -(2*dx(1)+dx(2))/(dx(1)*(dx(1)+dx(2))) - epsilon_s_left^2/(lambda_s*epsilon_1^2);
A(1,2) = (dx(1)+dx(2))/(dx(1)*dx(2));
A(1,3) = -dx(1)/(dx(2)*(dx(1)+dx(2)));
b(1) = 0;
% the PDE at internal nodes
%  epsilon_1^2 d^2 phi/dx^2 = -1/2 sum_i(c_i)
for i=2:N1-1
    A(i,i-1)=2/(dx(i-1)*(dx(i)+dx(i-1)));
    A(i,i)=-2/(dx(i)*dx(i-1));
    A(i,i+1)=2/(dx(i)*(dx(i)+dx(i-1)));
    b(i)=-1/(2*epsilon_1^2)*(z_cp*cp(i) + z_cm* cm(i));
end
%
% internal boundary condition
%     epsilon_1^2 phi_x_L(L1) = epsilon_2^2 phi_x_R(L1)
%     epsilon_1^2 phi_x_L(L1) - epsilon_2^2 phi_x_R(L1) = 0
% where phi_x_L is the derivative approximated from the left and phi_x_R
% is the derivative approximate from the right.
%
A(N1, N1-2) = epsilon_1^2*dx(N1-1)/(dx(N1-2)*(dx(N1-2)+dx(N1-1)));
A(N1, N1-1) = -epsilon_1^2*(dx(N1-2)+dx(N1-1))/(dx(N1-2)*dx(N1-1));
A(N1, N1) = epsilon_1^2*(2*dx(N1-1)+dx(N1-2))/(dx(N1-1)*(dx(N1-2)+dx(N1-1))) ...
    + epsilon_2^2*(2*dx(N1)+dx(N1+1))/(dx(N1)*(dx(N1)+dx(N1+1)));
A(N1, N1+1) = - epsilon_2^2*(dx(N1)+dx(N1+1))/(dx(N1)*dx(N1+1));
A(N1, N1+2) = epsilon_2^2 * (dx(N1))/(dx(N1+1)*(dx(N1)+dx(N1+1)));
b(N1)=0;
% the PDE at internal nodes
%  epsilon_2^2 d^2 phi/dx^2 = -1/2 sum_i(c_i)
for i=N1+1:N1+N2-1
    A(i,i-1)=2/(dx(i-1)*(dx(i)+dx(i-1)));
    A(i,i)=-2/(dx(i)*dx(i-1));
    A(i,i+1)=2/(dx(i)*(dx(i)+dx(i-1)));
    b(i)=-1/(2*epsilon_2^2)*(z_cp*cp(i) + z_cm*cm(i));
end
%
% internal boundary condition
%     epsilon_2^2 phi_x_L(L1+L2) = epsilon_3^2 phi_x_R(L1+L2)
%     epsilon_2^2 phi_x_L(L1+L2) - epsilon_3^2 phi_x_R(L1+L2) = 0
% where phi_x_L is the derivative approximated from the left and phi_x_R
% is the derivative approximate from the right.
%
A(N1+N2,N1+N2-2) = epsilon_2^2*(dx(N1+N2-1))/(dx(N1+N2-2)*(dx(N1+N2-2)+dx(N1+N2-1)));
A(N1+N2,N1+N2-1) = -epsilon_2^2*(dx(N1+N2-2)+dx(N1+N2-1))/(dx(N1+N2-2)*dx(N1+N2-1));
A(N1+N2,N1+N2) = epsilon_2^2*(2*dx(N1+N2-1)+dx(N1+N2-2))/(dx(N1+N2-1)*(dx(N1+N2-2) + dx(N1+N2-1))) ...
    + epsilon_3^2*(2*dx(N1+N2) + dx(N1+N2+1))/(dx(N1+N2)*(dx(N1+N2)+dx(N1+N2+1)));
A(N1+N2,N1+N2+1) =-epsilon_3^2*(dx(N1+N2) + dx(N1+N2+1))/(dx(N1+N2)*dx(N1+N2+1));
A(N1+N2,N1+N2+2) = epsilon_3^2*(dx(N1+N2))/(dx(N1+N2+1)*(dx(N1+N2)+dx(N1+N2+1)));
b(N1+N2)=0;
%
% the PDE at internal nodes
%  epsilon_3^2 d^2 phi/dx^2 = -1/2 sum_i(c_i)
%
for i=N1+N2+1:N-1
    A(i,i-1)=2/(dx(i-1)*(dx(i)+dx(i-1)));
    A(i,i)=-2/(dx(i)*dx(i-1));
    A(i,i+1)=2/(dx(i)*(dx(i)+dx(i-1)));
    b(i)=-1/(2*epsilon_3^2)*(z_cp*cp(i) + z_cm*cm(i));
end

% the boundary condition at x=L1+L2+L3
if bc==0 %current boundary condition
    %
    %  epsilon_3^2 phi_x(L1+L2+L3) = epsilon_s_right^2 phi_x_right
    %  phi_x(L1+L2+L3) = (epsilon_s_right/epsilon_3)^2 phi_x_right
    %
    A(N,N) = (2*dx(N-1)+dx(N-2))/(dx(N-1)*(dx(N-1)+dx(N-2)));
    A(N,N-1) = -(dx(N-2)+dx(N-1))/(dx(N-2)*dx(N-1));
    A(N,N-2) = dx(N-1)/(dx(N-2)*(dx(N-1)+dx(N-2)));
    b(N) = epsilon_s_right^2/epsilon_3^2*phi_x_right;
else %voltage boundary condition
    %
    % epsilon_3^2 phi_x(L) = epsilon_s_right^2 (v - phi(L)/(lambda_s)
    % phi_x(L) = (epsilon_s_right/epsilon_3)^2 (v - phi(L))/(lambda_s)
    % phi_x(L) + (epsilon_s_right/epsilon_3)^2 phi(L)/(lambda_s)
    %                      = (epsilon_s_right/epsilon_3)^2 v/(lambda_s)
    %
    A(N,N) = (2*dx(N-1)+dx(N-2))/(dx(N-1)*(dx(N-1)+dx(N-2))) + epsilon_s_right^2/(epsilon_3^2*lambda_s);
    A(N,N-1) = -(dx(N-2)+dx(N-1))/(dx(N-2)*dx(N-1));
    A(N,N-2) = dx(N-1)/(dx(N-2)*(dx(N-1)+dx(N-2)));
    b(N) = epsilon_s_right^2/(epsilon_3^2*lambda_s)*v;
end

%compute phi
phi = A \ b;


end
