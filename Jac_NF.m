function JacNF = Jac_NF(h,phi,cp,cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
% NF is in R^(3*N+1).  we want Jac_NF which is (3*N+1)x(3*N+1).
% first compute derivatives wrt phi then wrt cp then wrt cm then wrt
% phi_x_right.

N = N1+N2+N3;
% Now that I have a steady state, I can linearize about it.
A = eye(N); % use these for the perturbations of phi, cp, and cm
% h = .1;
% For current boundary conditions, we need to solve an ODE for phi_x_right
% and so JacNF is (3N+1)x(3N+1).
if bc == 0
    JacNF = zeros(3*N+1,3*N+1);
else
    JacNF = zeros(3*N,3*N);
end
for ii=1:N
    % derivatives wrt phi
    f1 = Nonlinear_F(phi-h*A(:,ii),cp,cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    f2 = Nonlinear_F(phi+h*A(:,ii),cp,cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    a1 = (f2-f1)/(2*h);
    JacNF(:,ii)=a1;
    
    % derivatives wrt cp
    f1 = Nonlinear_F(phi,cp-h*A(:,ii),cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    f2 = Nonlinear_F(phi,cp+h*A(:,ii),cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    a1 = (f2-f1)/(2*h);
    JacNF(:,N+ii)=a1;
    
    % derivatives wrt cm
    f1 = Nonlinear_F(phi,cp,cm-h*A(:,ii),phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    f2 = Nonlinear_F(phi,cp,cm+h*A(:,ii),phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    a1 = (f2-f1)/(2*h);
    JacNF(:,2*N+ii)=a1;
    
end

if bc == 0
    % derivative wrt phi_x_right
    f1 = Nonlinear_F(phi,cp,cm,phi_x_right-h,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    f2 = Nonlinear_F(phi,cp,cm,phi_x_right+h,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    a1 = (f2-f1)/(2*h);
    JacNF(:,3*N+1)=a1;
end








