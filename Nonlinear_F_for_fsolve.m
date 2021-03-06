function NF = Nonlinear_F_for_fsolve(y,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
% the input is y = [phi;cp;cm;phi_x_right] (if bc==0) and it's
% y = [phi,cp,cm] if bc==1.
%
% the output of this function is a 3*N+1 vector NF = [cp_RHS-cp_LHS;cm_RHS-cm_LHS;phi_RHS-phi_LHS,phi_x_right_RHS-phi_x_right_LHS].
%
% if bc==0 and we have current boundary conditions then NF is a 3N+1 vector
% because the +1 comes from the ODE that we need to solve for phi_x_right.
% If bc==1 and we have voltage boundary conditions then NF is a 3N vector.
%
% the free variables are cp,cm,phi,phi_x_right (playing the roles of
% cp_new,cm_new,phi_new,phi_x_right_new) and parameters include cp_now,
% cp_old, cm_now, cm_old, phi_x_right_now, and phi_x_right_old but these
% are *only* used in the _LHS portions which are the time derivatives.

w = dt_now/dt_old; %omega_(n)
N=N1+N2+N3; % total number of meshpoints
if bc==0
    phi = y(1:N);
    cp = y(N+1:2*N);
    cm = y(2*N+1:3*N);
    phi_x_right = y(3*N+1);
else
    phi = y(1:N);
    cp = y(N+1:2*N);
    cm = y(2*N+1:3*N);
    phi_x_right = rand; % it doesn't matter what value phi_x_right has and I need to put some value here because it's needed by Nonilnear_F_phi
end




% % step 1: given cm_now, cp_now, and v_now/phi_x_right_now compute phi_now.
% %         given cm_old, cp_old, and v_old/phi_x_right_old compute phi_old.
%
% phi_now = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_now,cm_now,z_cp,z_cm, phi_x_right_now,v_now);
% phi_old = poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_old,cm_old,z_cp,z_cm, phi_x_right_old,v_old);

% step 2: define the currents that'll be needed for the concentration
% equations
%
% voltage at time t_now+dt_now is needed in the Poisson problem and also in
% the BV boundary conditions...
v = voltage(t_now+dt_now);

NF(1:N,1) = Nonlinear_F_phi(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp,cm,z_cp,z_cm,phi, phi_x_right,v);

% Now create the fluxes that will be used for the concentration equations
j_m_left=0;
j_m_right=0;
BV_left = 4*kc_left*cp(1)*exp(phi(1)/(2*1 )) ...
    - 4*jr_left*exp(-phi(1)/(2*1 ));
j_p_left= - BV_left;
if bc==0 %for current boundary conditions
    BV_right = 4*kc_right*cp(N)*exp(-lambda_s*phi_x_right/(2*1 )) ...
        - 4*jr_right*exp(lambda_s*phi_x_right/(2*1 ));
    j_p_right = BV_right;
    
    % here's the ODE for phi_x_right when we have current boundary
    % conditions
    j_new = cur(t_now+dt_now);
    RHS = -(2/epsilon_s_right^2)*(j_new - (1+w)*(1/4)*z_cp*j_p_right + w*(1/4)*z_cp*j_p_right);
    LHS = (((1+2*w)/(1+w))*phi_x_right - (1+w)*phi_x_right_now + (w^2/(1+w))*phi_x_right_old)/dt_now;
    NF(3*N+1,1) = RHS-LHS;
    % NF(3*N+1,1) = ((1+w)^2/(1+2*w))*phi_x_right_now - (w^2/(1+2*w))*phi_x_right_old + dt_now*((1+w)/(1+2*w))*RHS - phi_x_right;
elseif bc==1
    BV_right = 4*kc_right*cp(N)*exp(-(v-phi(N))/(2*1 )) ...
        - 4*jr_right*exp((v-phi(N))/(2*1 ));
    j_p_right = BV_right;
    
    % there's no ODE for phi_x_right when there're voltage boundary
    % conditions
    % NF(3*N+1,1) = 0;
end
% use these fluxes in the PDEs for cp and cm
NF(N+1:3*N,1) = Nonlinear_F_c(dt_now,dt_old,phi,cp,cp_now,cp_old,cm,cm_now,cm_old,j_m_left,j_m_right,j_p_left,j_p_right, dx,N1,N2,N3,D_0,D_p,D_m,z_cp,z_cm);


