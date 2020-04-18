function NF = Nonlinear_F_c_one_step(dt_now,phi,cp,cp_now,cm,cm_now,j_m_left,j_m_right,j_p_left,j_p_right, dx,N1,N2,N3,D_0,D_p,D_m,z_cp,z_cm)
% the output of this function is a 2*N vector b = [cp_RHS-cp_LHS;cm_RHS-cm_LHS].  This
% woudl be writing the PDE as X_t = [cp_t;cm_t] = [cp_RHS;cm_RHS] = b -->
% NF = [cp_RHS-cp_t; cm_RHS-cm_t];


N=N1+N2+N3; % total number of meshpoints

% step 1: given cm_now, cp_now, and v_now/phi_x_right_now compute phi_now.
%         given cm_old, cp_old, and v_old/phi_x_right_old compute phi_old. 
% 
% phi_now = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_now,cm_now,z_cp,z_cm, phi_x_right_now,v_now);
% phi_old = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_old,cm_old,z_cp,z_cm, phi_x_right_old,v_old);

% step 2: define the currents that'll be needed for the concentration
% equations

% BV_left_now = 4*kc_left*cp_now(1)*exp(phi_now(1)/(2*1 )) ...
%     - 4*jr_left*exp(-phi_now(1)/(2*1 ));
% j_p_left_now= - BV_left_now;
% j_m_left_now=0;
% j_m_right_now=0;
% 
% BV_left_old = 4*kc_left*cp_old(1)*exp(phi_old(1)/(2*1 )) ...
%     - 4*jr_left*exp(-phi_old(1)/(2*1 ));
% j_p_left_old= - BV_left_old;
% j_m_left_old=0;
% j_m_right_old=0;
% if bc==0 %for current boundary conditions
%     BV_right_now = 4*kc_right*cp_now(N)*exp(-lambda_s*phi_x_right_now/(2*1 )) ...
%         - 4*jr_right*exp(lambda_s*phi_x_right_now/(2*1 ));
%     j_p_right_now = BV_right_now;
%     BV_right_old = 4*kc_right*cp_old(N)*exp(-lambda_s*phi_x_right_old/(2*1 )) ...
%         - 4*jr_right*exp(lambda_s*phi_x_right_old/(2*1 ));
%     j_p_right_old = BV_right_old;
% elseif bc==1
%     BV_right_now = 4*kc_right*cp_now(N)*exp(-(v_now-phi_now(N))/(2*1 )) ...
%         - 4*jr_right*exp((v_now-phi_now(N))/(2*1 ));
%     j_p_right_now = BV_right_now;
%     BV_right_old = 4*kc_right*cp_old(N)*exp(-(v_old-phi_old(N))/(2*1 )) ...
%         - 4*jr_right*exp((v_old-phi_old(N))/(2*1 ));
%     j_p_right_old = BV_right_old;
% end

% step 3: now compute cp_RHS and cm_RHS.  I can use the previous routine by
% simply taking phi_old = phi_now = phi, and so forth.  
% cp_RHS = concentration_implicit_RHS(z_cp,dx,N1,N2,N3,dt_now,dt_old,D_0,D_p,lambda_s,nu,j_p_left_now,j_p_right_now,j_p_left_old,j_p_right_old,cp_new,cp_now,phi_now,cp_old,phi_old,cp_now,cm_now,cp_old,cm_old );
% cm_RHS = concentration_implicit_RHS(z_cm,dx,N1,N2,N3,dt_now,dt_old,D_0,D_m,lambda_s,nu,j_m_left_now,j_m_right_now,j_m_left_old,j_m_right_old,cm_new,cm_now,phi_now,cm_old,phi_old,cp_now,cm_now,cp_old,cm_old );
cp_RHS = concentration_implicit_RHS(z_cp,dx,N1,N2,N3,D_0,D_p,j_p_left,j_p_right,cp,phi);
cm_RHS = concentration_implicit_RHS(z_cm,dx,N1,N2,N3,D_0,D_m,j_m_left,j_m_right,cm,phi);

% we have backward Euler time-stepping: (c_new-c_now)/dt_now = RHS
cm_LHS = (cm - cm_now)/dt_now;
cp_LHS = (cp - cp_now)/dt_now;

NF(1:N,1) = cp_RHS-cp_LHS;
NF(N+1:2*N,1) = cm_RHS-cm_LHS;
