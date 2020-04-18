function [ err ] = errcomp( step,dt_now,dt_old,cp_coarse,cm_coarse,phi_x_right_coarse,cp_fine,cm_fine,phi_x_right_fine)
%errcomp.m
%computes error value based on richardson extrapolation
%step = 0 means first time step, step = 1 means subsequent steps
%make sure to set phi_x_right_coarse and phi_x_right_fine to 0 if using
%voltage boundary conditions

if step==0
    err=4/3*norm(cp_coarse-cp_fine);
else

% err=norm(cp_coarse-cp_fine)*8*(dt_old+dt_now)/(7*dt_old+5*dt_now);
err=norm(cat(1,cp_coarse,cm_coarse,phi_x_right_coarse) ...
    -cat(1,cp_fine,cm_fine,phi_x_right_fine))*8*(dt_old+dt_now)/(7*dt_old+5*dt_now);
end

end

