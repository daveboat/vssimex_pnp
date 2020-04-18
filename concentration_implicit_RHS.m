function [ RHS ] = concentration_implicit_RHS( z,dx,N1,N2,N3,D_0,D,j_left,j_right,c_new,phi)
%concentration.m
%creates c_new depending on value of z (+1 for cp, -1 for cm, etc)
%works for both types of concentrations
%j_left_(old/now) and j_right_(old/now) are calculated in step.m and passed in. They
%represent the values of the fluxes at the previous/current time.
%
% This is simulating the RHS of the PDEs
% 
% cp_t = (D_p/D_0)[cp_x +  cp phi_x]_x   where D_p is passed in as D
% cm_t = (D_m/D_0)[cm_x  - cm phi_x]_x   where D_m is passed in as D  
% 
% it's being done fully implicitly with phi having been extrapolated
% forward to time t_new.  

N=N1+N2+N3;
A=zeros(N,N);
RHS = zeros(N,1);

% set up A and b:
% at the internal nodes we have the PDE dc/dt = D/D0 [c_x + z c phi_x]_x
% I'm coding this somewhat kloodgily in that I write the matrix A for the
% operator on the RHS and I then subsequently use it in defining the matrix
% that needs to be inverted.
for i=2:N-1    
    prefactor = (D/D_0)*(2/(dx(i-1)+dx(i)));
    phi_x_a = (z/2)*(phi(i)-phi(i-1))/dx(i-1);
    phi_x_b = (z/2)*(phi(i+1)-phi(i))/dx(i);
    
    A(i,i-1) = prefactor * (1/dx(i-1) - phi_x_a);
    A(i,i) = prefactor * (-1/dx(i-1) - 1/dx(i) + phi_x_b - phi_x_a);
    A(i,i+1) = prefactor * (1/dx(i) + phi_x_b);
   
    RHS(i,1) = A(i,i-1)*c_new(i-1) + A(i,i)*c_new(i) + A(i,i+1)*c_new(i+1);
end



%handle boundary conditions via ghost points.

%lhs
A(1,1) = -1/dx(1) + (z/2)*(phi(2)-phi(1))/dx(1);
A(1,2) = 1/dx(1) + (z/2)*(phi(2)-phi(1))/dx(1);
A(1,:) = (D/D_0)*(2/dx(1))*A(1,:);
RHS(1,1) = A(1,:)*c_new;
RHS(1,1) = RHS(1) + 1/(dx(1)/2)*D/D_0*(j_left);


%rhs
A(N,N-1)= 1/dx(N-1) - (z/2)*(phi(N)-phi(N-1))/dx(N-1);
A(N,N)= -1/dx(N-1) - (z/2)*(phi(N)-phi(N-1))/dx(N-1);
A(N,:) = (D/D_0)*(2/dx(N-1))*A(N,:);
RHS(N,1) = A(N,:)*c_new;
RHS(N,1)= RHS(N,1) - 1/(dx(N-1)/2)*D/D_0*(j_right);

end
