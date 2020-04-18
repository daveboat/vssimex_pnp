function [B_new,B_now,B_old] = auto_linearize_RHS(h,test,v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm)


N = N1+N2+N3;
dt_now=rand;
dt_old=dt_now;
t=rand;
% Now that I have a steady state, I can linearize about it.
A = eye(N);
% h = .1;

for ii=1:N
    
    
    f1 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f2 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new+h*A(:,ii),cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f3 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new+2*h*A(:,ii),cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f4 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new+3*h*A(:,ii),cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    if isequal(test,'linear')
        a1 = (f2-f1)/h;
    elseif isequal(test,'quadratic')
        a1 = -((3*f1 - 4*f2 + f3)/(2*h));
    elseif isequal(test,'cubic')
        a1 = -(11*f1 - 18*f2 + 9*f3 - 2*f4)/(6*h);
    else
        display('choose flag correctly')
    end
    B_new(:,ii)=a1;

    
    f2 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new+h*A(:,ii),cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f3 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new+2*h*A(:,ii),cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f4 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new+3*h*A(:,ii),cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
        
    if isequal(test,'linear')
        a1 = (f2-f1)/h;
    elseif isequal(test,'quadratic')
        a1 = -((3*f1 - 4*f2 + f3)/(2*h));
    elseif isequal(test,'cubic')
        a1 = -(11*f1 - 18*f2 + 9*f3 - 2*f4)/(6*h);
    else
        display('choose flag correctly')
    end
    B_new(:,N+ii)=a1;
    
      
    f2 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now+h*A(:,ii),cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f3 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now+2*h*A(:,ii),cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f4 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now+3*h*A(:,ii),cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    if isequal(test,'linear')
        a1 = (f2-f1)/h;
    elseif isequal(test,'quadratic')
        a1 = -((3*f1 - 4*f2 + f3)/(2*h));
    elseif isequal(test,'cubic')
        a1 = -(11*f1 - 18*f2 + 9*f3 - 2*f4)/(6*h);
    else
        display('choose flag correctly')
    end
    B_now(:,ii)=a1;

    
    f2 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now+h*A(:,ii),cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f3 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now+2*h*A(:,ii),cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f4 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now+3*h*A(:,ii),cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
        
    if isequal(test,'linear')
        a1 = (f2-f1)/h;
    elseif isequal(test,'quadratic')
        a1 = -((3*f1 - 4*f2 + f3)/(2*h));
    elseif isequal(test,'cubic')
        a1 = -(11*f1 - 18*f2 + 9*f3 - 2*f4)/(6*h);
    else
        display('choose flag correctly')
    end
    B_now(:,N+ii)=a1;
    
      
    f2 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old+h*A(:,ii),cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f3 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old+2*h*A(:,ii),cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f4 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old+3*h*A(:,ii),cm_ss_new,cm_ss_now,cm_ss_old,bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    if isequal(test,'linear')
        a1 = (f2-f1)/h;
    elseif isequal(test,'quadratic')
        a1 = -((3*f1 - 4*f2 + f3)/(2*h));
    elseif isequal(test,'cubic')
        a1 = -(11*f1 - 18*f2 + 9*f3 - 2*f4)/(6*h);
    else
        display('choose flag correctly')
    end
    B_old(:,ii)=a1;

    
    f2 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old+h*A(:,ii),bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f3 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old+2*h*A(:,ii),bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
    f4 = RHS(v_now,v_old,dt_now,dt_old,phi_x_right_now,phi_x_right_old,cp_ss_new,cp_ss_now,cp_ss_old,cm_ss_new,cm_ss_now,cm_ss_old+3*h*A(:,ii),bc, dx,N1,N2,N3, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    
        
    if isequal(test,'linear')
        a1 = (f2-f1)/h;
    elseif isequal(test,'quadratic')
        a1 = -((3*f1 - 4*f2 + f3)/(2*h));
    elseif isequal(test,'cubic')
        a1 = -(11*f1 - 18*f2 + 9*f3 - 2*f4)/(6*h);
    else
        display('choose flag correctly')
    end
    B_old(:,N+ii)=a1;
    
 
end




