function [ j_f, j_d, j, v ] = postprocess( bc,dx,N1,N2,N3,t,t_old,cp,cm,phi,phi_old,phi_x_right,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)
%postprocess for either j or v depending on the boundary conditions.
%remember that dx and t aren't uniform.

N=N1+N2+N3;
if bc==0
        v=phi(N) + epsilon_3^2/epsilon_s_right^2*lambda_s*phi_x_right;
        j=cur(t);
        j_f=0;
        j_d=0;
else

    jpr_now = 4*kc_right*cp(N)*exp(-(voltage(t)-phi(N))/(2*1))...
            - 4*jr_right*exp((voltage(t)-phi(N))/(2*1));
    jmr_now = 0;
    v=voltage(t);
    j_f=(z_cp*jpr_now + z_cm*jmr_now)/4;
    j_d=- epsilon_s_right^2/(2*(t-t_old))*((voltage(t) - phi(N))/lambda_s - (voltage(t_old) - phi_old(N))/lambda_s);
    j=j_f+j_d;
end

end
