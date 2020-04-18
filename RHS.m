function b = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_new,cp_now,cp_old,cm_new,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
% the output of this function is a 2*N vector b = [cp_RHS;cm_RHS].  This
% woudl be writing the PDE as X_t = [cp_t;cm_t] = [cp_RHS;cm_RHS] = b


w = dt_now/dt_old; %omega_(n)
N=N1+N2+N3; % total number of meshpoints

% step 1: given cm_now, cp_now, and v_now/phi_x_right_now compute phi_now.
%         given cm_old, cp_old, and v_old/phi_x_right_old compute phi_old. 

phi_now = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_now,cm_now,z_cp,z_cm, phi_x_right_now,v_now);
phi_old = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_old,cm_old,z_cp,z_cm, phi_x_right_old,v_old);

% step 2: define the currents that'll be needed for the concentration
% equations

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
if bc==0 %for current boundary conditions
    BV_right_now = 4*kc_right*cp_now(N)*exp(-lambda_s*phi_x_right_now/(2*1 )) ...
        - 4*jr_right*exp(lambda_s*phi_x_right_now/(2*1 ));
    j_p_right_now = BV_right_now;
    BV_right_old = 4*kc_right*cp_old(N)*exp(-lambda_s*phi_x_right_old/(2*1 )) ...
        - 4*jr_right*exp(lambda_s*phi_x_right_old/(2*1 ));
    j_p_right_old = BV_right_old;
elseif bc==1
    BV_right_now = 4*kc_right*cp_now(N)*exp(-(v_now-phi_now(N))/(2*1 )) ...
        - 4*jr_right*exp((v_now-phi_now(N))/(2*1 ));
    j_p_right_now = BV_right_now;
    BV_right_old = 4*kc_right*cp_old(N)*exp(-(v_old-phi_old(N))/(2*1 )) ...
        - 4*jr_right*exp((v_old-phi_old(N))/(2*1 ));
    j_p_right_old = BV_right_old;
end

% step 3: now compute cp_RHS and cm_RHS
cp_RHS = concentration_RHS(z_cp,dx,N1,N2,N3,dt_now,dt_old,D_0,D_p,lambda_s,j_p_left_now,j_p_right_now,j_p_left_old,j_p_right_old,cp_new,cp_now,phi_now,cp_old,phi_old,cp_now,cm_now,cp_old,cm_old );
cm_RHS = concentration_RHS(z_cm,dx,N1,N2,N3,dt_now,dt_old,D_0,D_m,lambda_s,j_m_left_now,j_m_right_now,j_m_left_old,j_m_right_old,cm_new,cm_now,phi_now,cm_old,phi_old,cp_now,cm_now,cp_old,cm_old );

b(1:N,1) = cp_RHS;
b(N+1:2*N,1) = cm_RHS;
