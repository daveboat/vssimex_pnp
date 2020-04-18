function [ RHS ] = concentration_RHS( z,dx,N1,N2,N3,dt_now,dt_old,D_0,D,lambda_s,j_left_now,j_right_now,j_left_old,j_right_old,c_new,c_now,phi_now,c_old,phi_old,cp_now,cm_now,cp_old,cm_old )
%concentration_RHS.m
%creates c_new depending on value of z (+1 for cp, -1 for cm, etc)
%works for both types of concentrations
%j_left_(old/now) and j_right_(old/now) are calculated in step.m and passed in. They
%represent the values of the fluxes at the previous/current time.
%
% This is simulating the RHS of the PDEs
% 
% cp_t = (D/D_0)[cp_x +  cp phi_x ]_x
% cm_t = (D/D_0)[cm_x  - cm phi_x ]_x
% 
% I.e. it's the usual flux of c + z c phi_x
% but with an additional term of
% 
% - cp d/dx in the cp_t equation and
% - cm d/dx in the cp_t equation 
%
% verified that it's correct using test_concentration_RHS.m

w=dt_now/dt_old;
N=N1+N2+N3;

%set up A and b
% at the internal nodes we have the PDE dc/dt = D/D0 [c_x + z c phi_x - c (ln(1-nu(cp+cm)))_x]_x
for i=2:N-1
    RHS(i)=D/D_0*(2/(dx(i-1)*(dx(i)+dx(i-1))))*c_new(i-1) - D/D_0*(2/(dx(i)*dx(i-1)))*c_new(i)...
            + D/D_0*(2/(dx(i)*(dx(i)+dx(i-1))))*c_new(i+1)...
            + (1+w)*D/D_0*( ( z*(c_now(i)+c_now(i+1))/(1+1)*(phi_now(i+1)-phi_now(i))/dx(i)...
               - z*(c_now(i)+c_now(i-1))/(1+1)*(phi_now(i)-phi_now(i-1))/dx(i-1) ) / ((dx(i)+dx(i-1))/2)) ...
              - w*D/D_0*( ( z*(c_old(i)+c_old(i+1))/(1+1)*(phi_old(i+1)-phi_old(i))/dx(i) ...
                 - z*(c_old(i)+c_old(i-1))/(1+1 )*(phi_old(i)-phi_old(i-1))/dx(i-1) ) / ((dx(i)+dx(i-1))/2));
end

%handle boundary conditions via ghost points.
RHS(1)= -(1/(dx(1)*(dx(1)/2)))*c_new(1) + (1/(dx(1)*(dx(1)/2)))*c_new(2) ...
    + (1+w)*1/(dx(1)/2)*D/D_0*( z*(c_now(2)+c_now(1))/(1+1 ) * (phi_now(2)-phi_now(1))/dx(1) + j_left_now ) ...
    - w*1/(dx(1)/2)*D/D_0*( z*(c_old(2)+c_old(1))/(1 +1 ) * (phi_old(2)-phi_old(1))/dx(1) + j_left_old );

%rhs
RHS(N)=-(1/(dx(N-1)*(dx(N-1)/2)))*c_new(N) +(1/(dx(N-1)*(dx(N-1)/2)))*c_new(N-1) ...
    - (1+w)*1/(dx(N-1)/2)*D/D_0*( z*(c_now(N-1)+c_now(N))/( 1+1 ) * (phi_now(N)-phi_now(N-1))/dx(N-1) + j_right_now ) ...
    + w*1/(dx(N-1)/2)*D/D_0*( z*(c_old(N-1)+c_old(N))/(1 + 1 ) * (phi_old(N)-phi_old(N-1))/dx(N-1) + j_right_old );

end
