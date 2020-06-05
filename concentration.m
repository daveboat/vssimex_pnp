function [ c_new ] = concentration( z,dx,N1,N2,N3,dt_now,dt_old,D_0,D,lambda_s,j_left_now,j_right_now,j_left_old,j_right_old,c_now,phi_now,c_old,phi_old,cp_now,cm_now,cp_old,cm_old )
%concentration.m
%creates c_new depending on value of z (+1 for cp, -1 for cm, etc)
%works for both types of concentrations
%j_left_(old/now) and j_right_(old/now) are calculated in step.m and passed in. They
%represent the values of the fluxes at the previous/current time.
%
% This is simulating the PDEs
% 
% cp_t = (D/D_0)[cp_x +  cp phi_x + nu cp (cp_x + cm_x)/(1-nu(cp+cm)) ]_x
% cm_t = (D/D_0)[cm_x  - cm phi_x + nu cm (cp_x + cm_x)/(1-nu(cp+cm)) ]_x
% 
% I.e. it's the usual flux of c + z c phi_x
% but with an additional term of
% 
% - cp d/dx ln( 1 - nu(cp+cm) ) in the cp_t equation and
% - cm d/dx ln( 1 - nu(cp+cm) ) in the cp_t equation 

w=dt_now/dt_old;
N=N1+N2+N3;

A=zeros(N,N);
b=zeros(N,1);

%set up A and b
% at the internal nodes we have the PDE dc/dt = D/D0 [c_x + z c phi_x - c (ln(1-nu(cp+cm)))_x]_x
for i=2:N-1
    prefactor = dt_now*(D/D_0)*(2/(dx(i)+dx(i-1)));
    
    A(i,i-1)=-prefactor/dx(i-1);
    A(i,i)=(1+2*w)/(1+w) + dt_now*D/D_0*2/(dx(i)*dx(i-1));
    A(i,i+1)=-prefactor/dx(i);

    b(i)=(1+w)*c_now(i) - w^2/(1+w)*c_old(i) + prefactor*(1+w)*( ( z*(c_now(i)+c_now(i+1))/2 * (phi_now(i+1)-phi_now(i))/dx(i)- z*(c_now(i)+c_now(i-1))/(1+1)*(phi_now(i)-phi_now(i-1))/dx(i-1) ) ) - prefactor * w * ( ( z*(c_old(i)+c_old(i+1))/(1+1)*(phi_old(i+1)-phi_old(i))/dx(i) - z*(c_old(i)+c_old(i-1))/2 *(phi_old(i)-phi_old(i-1))/dx(i-1) ));
end

%lhs
A(1,1)=(1+2*w)/(1+w) + dt_now/(dx(1)*(dx(1)/2));
A(1,2)=-dt_now/(dx(1)*(dx(1)/2));
b(1)=(1+w)*c_now(1) - w^2/(1+w)*c_old(1) ...
    + (1+w)*dt_now/(dx(1)/2)*D/D_0*( z*(c_now(2)+c_now(1))/(1+1 ) * (phi_now(2)-phi_now(1))/dx(1) + j_left_now ) ...
    - w*dt_now/(dx(1)/2)*D/D_0*( z*(c_old(2)+c_old(1))/(1 +1 ) * (phi_old(2)-phi_old(1))/dx(1) + j_left_old );

%rhs
A(N,N)=(1+2*w)/(1+w) + dt_now/(dx(N-1)*(dx(N-1)/2));
A(N,N-1)=-dt_now/(dx(N-1)*(dx(N-1)/2));
b(N)=(1+w)*c_now(N) - w^2/(1+w)*c_old(N) ...
    - (1+w)*dt_now/(dx(N-1)/2)*D/D_0*( z*(c_now(N-1)+c_now(N))/( 1+1 ) * (phi_now(N)-phi_now(N-1))/dx(N-1)  + j_right_now ) ...
    + w*dt_now/(dx(N-1)/2)*D/D_0*( z*(c_old(N-1)+c_old(N))/(1 + 1 ) * (phi_old(N)-phi_old(N-1))/dx(N-1)  + j_right_old );


% Using the direct method rather than the ghost point method.  Have 
% flux = BC, i.e. D/D_0( - c_x - z c phi_x) = j 
% and so                   -c_x - z c phi_x = D_0/D j
% and so                               -c_x = z c phi_x + D_0/D j = RHS
% the RHS is extrapolated forward in time with (1+w) RHS_new - w RHS_old
% and the LHS is handled implicitly with -c_new_x.  The first derivatives
% are approximated using a two-point stencil
% A(1,2) = -1/dx(1);
% A(1,1) = 1/dx(1);
% b(1) = (1+w)*( z*((c_now(1)+c_now(2))/2)*(phi_now(2)-phi_now(1))/dx(1) + (D_0/D)*j_left_now)...
%     - w*( z*((c_old(1)+c_old(2))/2)*(phi_old(2)-phi_old(1))/dx(1) + (D_0/D)*j_left_old);
% 
% A(N,N) = - 1/dx(N-1);
% A(N,N-1) = 1/dx(N-1);
% b(N) = (1+w)*( z*((c_now(N-1)+c_now(N))/2)*(phi_now(N)-phi_now(N-1))/dx(N-1) + (D_0/D)*j_right_now)...
%     - w*( z*((c_old(N-1)+c_old(N))/2)*(phi_old(N)-phi_old(N-1))/dx(N-1) + (D_0/D)*j_right_old);


% % another approach to the direct method but using 3 point stencils...
% a0 = -(2*dx(1)+dx(2))/(dx(1)*(dx(1)+dx(2)));
% a1 = 1/dx(1) + 1/dx(2);
% a2 = - dx(1)/(dx(2)*(dx(1)+dx(2)));
% % the three point stencil approximation of f'(x) is a0 f(x0) + a1 f(x1) + a2 f(x2)
% % we want -c_new_x
% A(1,1) = - a0;
% A(1,2) = - a1;
% A(1,3) = - a2;
% phi_now_x = a0*phi_now(1) + a1*phi_now(2) + a2*phi_now(3);
% phi_old_x = a0*phi_old(1) + a1*phi_old(2) + a2*phi_old(3);
% b(1) = (1+w)*( z*((c_now(1)+c_now(2))/2)*phi_now_x + (D_0/D)*j_left_now)...
%     - w*( z*((c_old(1)+c_old(2))/2)*phi_old_x + (D_0/D)*j_left_old);
% %
% a0 = (2*dx(N-1)+dx(N-2))/(dx(N-1)*(dx(N-1)+dx(N-2)));
% am1 = -(1/dx(N-1) + 1/dx(N-2));
% am2 = dx(N-1)/(dx(N-2)*(dx(N-1)+dx(N-2)));
% % the three point stencil approximation of f'(x) is a0 f(x(N)) + am1 f(x(N-1)) + am2 f(x(N-2))
% % we want - c_new_x
% A(N,N-2) = -am2;
% A(N,N-1) = -am1;
% A(N,N) = -a0;
% phi_now_x = am2*phi_now(N-2)+am1*phi_now(N-1)+a0*phi_now(N);
% phi_old_x = am2*phi_old(N-2)+am1*phi_old(N-1)+a0*phi_old(N);
% b(N) = (1+w)*( z*((c_now(N-1)+c_now(N))/2)*phi_now_x + (D_0/D)*j_right_now)...
%     - w*( z*((c_old(N-1)+c_old(N))/2)*phi_old_x + (D_0/D)*j_right_old);

c_new=A\b;

end
