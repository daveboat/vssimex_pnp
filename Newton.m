function [phi_new,cp_new,cm_new,phi_x_right_new,NF] = Newton(t_now,dt_now,dt_old,phi_now,phi_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
% Hand-coded Newton-Raphson method for solving the nonlinear system at each
% implicit time step


h = .01; % used by the Jacobian approximator (f(x+h)-f(x-h))/(2h)
N=N1+N2+N3; % total number of meshpoints

tol = 1.e-10;
it_max = 10;

% first guess:
w = dt_now/dt_old;
cp = (1+w)*cp_now-w*cp_old;
cm = (1+w)*cm_now-w*cm_old;
phi = (1+w)*phi_now-w*phi_old;
phi_x_right = (1+w)*phi_x_right_now-w*phi_x_right_old;
NF = Nonlinear_F(phi,cp,cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
max(abs(NF));

it = 1;
while ((max(abs(NF))) > tol)&&(it < it_max)
    % the output of this function is a 3*N vector if bc=1 and a 3*N+1 vector if
    % bc=0
    % now I need to compute the Jacobian
    JacNF = Jac_NF(h,phi,cp,cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);

    % take one Newton step
    y(1:N,1) = phi;
    y(N+1:2*N,1) = cp;
    y(2*N+1:3*N,1) = cm;
    if bc == 0
        y(3*N+1,1) = phi_x_right;
    end
    % newton step is x_new = x_now - (2*x_now)^(-1)*(x_now^2-4);
    % y_new = y - JacNF\NF;
    y = y-JacNF\NF;
    phi = y(1:N,1);
    cp = y(N+1:2*N,1);
    cm = y(2*N+1:3*N,1);
    if bc == 0
        phi_x_right = y(3*N+1,1);
    end
    NF = Nonlinear_F(phi,cp,cm,phi_x_right,t_now,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_now,cp_old,cm_now,cm_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    max(abs(NF));
    it = it+1;
end

% once we've converged:
phi_new = phi;
cp_new = cp;
cm_new = cm;
phi_x_right_new = phi_x_right;