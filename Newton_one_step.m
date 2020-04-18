function [phi_new,cp_new,cm_new,phi_x_right_new,NF] = Newton_one_step(t_now,dt_now,phi_x_right_now,cp_now,cm_now,phi_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)

h = .01; % used by the Jacobian approximator (f(x+h)-f(x-h))/(2h)
N=N1+N2+N3; % total number of meshpoints

tol = 1.e-10;
it_max = 10;

% first guess:
cp = cp_now;
cm = cm_now;
phi = phi_now;
phi_x_right = phi_x_right_now;
NF = Nonlinear_F_one_step(phi,cp,cm,phi_x_right,t_now,dt_now,phi_x_right_now,cp_now,cm_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
max(abs(NF));

it = 1;
while ((max(abs(NF))) > tol)&&(it < it_max)
    % the output of this function is a 3*N vector if bc=1 and a 3*N+1 vector if
    % bc=0
    % now I need to compute the Jacobian
    JacNF = Jac_NF_one_step(h,phi,cp,cm,phi_x_right,t_now,dt_now,phi_x_right_now,cp_now,cm_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
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
    
    NF = Nonlinear_F_one_step(phi,cp,cm,phi_x_right,t_now,dt_now,phi_x_right_now,cp_now,cm_now,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    max(abs(NF));
    it = it+1;
end


% once we've converged:
phi_new = phi;
cp_new = cp;
cm_new = cm;
phi_x_right_new = phi_x_right;