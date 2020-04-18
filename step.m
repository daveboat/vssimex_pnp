function [ cp_new,cm_new,phi_x_right_new,debug1 ] = step(bc,t,dt_now,dt_old,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_old,cm_old,phi_old,phi_x_right_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
%
% Steps the vectors cp, cm, phi, and phi_x_right to the next time step.
% Uses a variable step size, nonuniform mesh implicit-explicit BDF time
% stepper for the concentrations, temperature and phi_x, and an O(dx^2)
% Poisson solver for phi.
%
% bc: 1 for voltage boundary conditions, 0 for current boundary conditions
% t: current time (scalar)
% dt_now: current time step (i.e. t[n+1] = t[n] + dt_now)
% dt_old: previous time step (i.e. t[n] = t[n-1] + dt_old)
% dx, N1,N2,N3: N-1 length vector of mesh widths (vector) and numbers of
% meshes in each region
% <vector>_now: N length vectors (or scalar for phi_x_right) holding current
% values in time
% <vector>_old: N length vectors (or scalar for phi_x_right) holding
% previous values in time (that is, at t[n] - dt_old)


w = dt_now/dt_old; %omega_(n)
N=N1+N2+N3; % total number of meshpoints
debug1=zeros(N,1);

% define the currents at the left and right electrodes; these are needed
% for time-stepping the concentrations and are also needed for
% time-stepping phi_x_right when the current is being imposed.

BV_left_now = 4*kc_left*cp_now(1)*exp(phi_now(1)/(2*1 )) ...
    - 4*jr_left*exp(-phi_now(1)/(2*1 ));
j_p_left_now= - BV_left_now;
j_m_left_now=0;
j_m_right_now=0;

BV_left_old = 4*kc_left*cp_old(1)*exp(phi_old(1)/(2*1 )) ...
    - 4*jr_left*exp(-phi_old(1)/(2*1 ));
j_p_left_old= - BV_left_old;
j_m_left_old=0;
j_m_right_old=0;


% time-step phi_x_right
% d/dt phi_x_right = -2/epsilon_s_right^2 [j - k_cr*cp*exp(- lambda_s phi_x_right/2) + jrr*exp(lambda_s phi_x_right/2) ]
%
% d/dt phi_x_right = RHS then
%
% 1/dt_now [ (1+2w)/(1+w) phi_x_right_new - (1+w) phi_x_right_now + w^2/(1+w) phi_x_right_old ] = RHS
%
% and so
%
% phi_x_right_new = (1+w)^2/(1+2*w) phi_x_right_now - w^2/(1+2*w) phi_x_right_old
%                      + dt_now (1+w)/(1+2*w) RHS
% where
% RHS = -2/epsilon_s_right^2( j_new + (1+w)(- k_cr*cp_now*exp(- lambda_s phi_x_right_now/2) + jrr*exp(lambda_s phi_x_right_now/2) )
%           - w (- k_cr*cp_old*exp(- lambda_s phi_x_right_old/2) + jrr*exp(lambda_s phi_x_right_old/2) ))
%     = -2/epsilon_s_right^2 (j_new - (1+w) (1/4) j_p_right_now + w (1/4)  j_p_right_old)
%
if bc==0 %for current boundary conditions
    BV_right_now = 4*kc_right*cp_now(N)*exp(-lambda_s*phi_x_right_now/(2*1 )) ...
        - 4*jr_right*exp(lambda_s*phi_x_right_now/(2*1 ));
    j_p_right_now = BV_right_now;
    BV_right_old = 4*kc_right*cp_old(N)*exp(-lambda_s*phi_x_right_old/(2*1 )) ...
        - 4*jr_right*exp(lambda_s*phi_x_right_old/(2*1 ));
    j_p_right_old = BV_right_old;
    
    j_new = cur(t+dt_now);
    RHS = -(2/epsilon_s_right^2)*(j_new - (1+w)*(1/4)*z_cp*j_p_right_now + w*(1/4)*z_cp*j_p_right_old);
    phi_x_right_new =((1+w)^2/(1+2*w))*phi_x_right_now - (w^2/(1+2*w))*phi_x_right_old + dt_now*((1+w)/(1+2*w))*RHS;
    
elseif bc==1
    BV_right_now = 4*kc_right*cp_now(N)*exp(-(voltage(t)-phi_now(N))/(2*1 )) ...
        - 4*jr_right*exp((voltage(t)-phi_now(N))/(2*1 ));
    j_p_right_now = BV_right_now;
    BV_right_old = 4*kc_right*cp_old(N)*exp(-(voltage(t-dt_old)-phi_old(N))/(2*1 )) ...
        - 4*jr_right*exp((voltage(t-dt_old)-phi_old(N))/(2*1 ));
    j_p_right_old = BV_right_old;
    
    phi_x_right_new=0;
end

%step cp,cm by calling their respective functions

cp_new=concentration( z_cp,dx,N1,N2,N3,dt_now,dt_old,D_0,D_p,lambda_s,j_p_left_now,j_p_right_now,j_p_left_old,j_p_right_old,cp_now,phi_now,cp_old,phi_old,cp_now,cm_now,cp_old,cm_old );
cm_new=concentration( z_cm,dx,N1,N2,N3,dt_now,dt_old,D_0,D_m,lambda_s,j_m_left_now,j_m_right_now,j_m_left_old,j_m_right_old,cm_now,phi_now,cm_old,phi_old,cp_now,cm_now,cp_old,cm_old );


v_now=voltage(t);
if dt_old==Inf
    v_old=0;
else
    v_old=voltage(t-dt_old);
end


end

