function [phi_new,cp_new,cm_new,phi_x_right_new,NF] = Fsolve_one_step(t_now,dt_now,phi_x_right_now,cp_now,cm_now,phi_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)

h = .01; % used by the Jacobian approximator (f(x+h)-f(x-h))/(2h)
N=N1+N2+N3; % total number of meshpoints

tol = 1.e-10;
it_max = 10;

% define the initial guess:
y0(1:N,1) = phi_now;
y0(N+1:2*N,1)= cp_now;
y0(2*N+1:3*N,1)= cm_now;
if bc==0
    y0(3*N+1,1) = phi_x_right_now;
end
[y,NF] = fsolve(@(y) Nonlinear_F_one_step_for_fsolve(y,t_now,dt_now,phi_x_right_now,cp_now,cm_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm),y0,optimset('Display','off'));
% [y,NF] = fsolve(@(y) Nonlinear_F_one_step_for_fsolve(y,t_now,dt_now,phi_x_right_now,cp_now,cm_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm),y0,optimset('Display','iter'));

% once we've converged:
phi_new = y(1:N,1);
cp_new = y(N+1:2*N,1);
cm_new = y(2*N+1:3*N,1);
if bc==0
    phi_x_right_new = y(3*N+1,1);
else
    phi_x_right_new = phi_x_right_now;
end