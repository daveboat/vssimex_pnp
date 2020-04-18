debug1=zeros(N,1);


t_now=t(n);
cp_now=cp(:,n);
cm_now=cm(:,n);
phi_x_right_now=phi_x_right(n);
phi_now=phi(:,n);

[j_f(n),j_d(n),j(n),v(n)]=postprocess(bc,dx,N1,N2,N3,t(n),Inf,cp(:,n),cm(:,n),phi(:,n),zeros(N,1),phi_x_right(n),lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);


% flag=0;

%first time step
%coarse step
[cp_coarse,cm_coarse,phi_x_right_coarse,debug1_coarse ] = step(bc,t_now,dt,Inf,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,ones(N,1),ones(N,1),ones(N,1),0,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
%fine step
[ cp_half,cm_half,phi_x_right_half,debug1_half ] = step(bc,t_now,dt/2,Inf,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,ones(N,1),ones(N,1),ones(N,1),0,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
phi_half=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_half,cm_half,z_cp,z_cm, phi_x_right_half,voltage(t_now+dt/2) );
[ cp_fine,cm_fine,phi_x_right_fine,debug1_fine ] = step(bc,t_now+dt/2,dt/2,dt/2,dx,N1,N2,N3,cp_half,cm_half,phi_half,phi_x_right_half,cp_now,cm_now,phi_now,phi_x_right_now,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
%evaluate error
err=errcomp( 0,dt,0,cp_coarse,cm_coarse,phi_x_right_coarse,cp_fine,cm_fine,phi_x_right_fine);
while abs(err - tol)>range 
    %take a stab at the new dt
    fac=(tol/err)^(.5);
    if fac > eta_max fac=eta_max;
    elseif fac < eta_min fac=eta_min;
    end
    dt=dt*fac;
    %If my calculated dt goes over dt_max, take a step with dt_max and exit
    %the loop.
    if dt>=dt_max
        dt=dt_max;
        [ cp_coarse,cm_coarse,phi_x_right_coarse,debug1_coarse ] = step(bc,t_now,dt,Inf,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,ones(N,1),ones(N,1),ones(N,1),ones(N,1),0,0,0, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        [ cp_half,cm_half,phi_x_right_half,debug1_half ] = step(bc,t_now,dt/2,Inf,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,ones(N,1),ones(N,1),ones(N,1),ones(N,1),0,0,0, lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        phi_half=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_half,cm_half,z_cp,z_cm, phi_x_right_half,voltage(t_now+dt/2) );
        [ cp_fine,cm_fine,phi_x_right_fine,debug1_fine ] = step(bc,t_now+dt/2,dt/2,dt/2,dx,N1,N2,N3,cp_half,cm_half,phi_half,phi_x_right_half,cp_now,cm_now,phi_now,phi_x_right_now,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        break;
    end
    %coarse step
    [ cp_coarse,cm_coarse,phi_x_right_coarse,debug1_coarse ] = step(bc,t_now,dt,Inf,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,ones(N,1),ones(N,1),ones(N,1),0,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    %fine step
    [ cp_half,cm_half,phi_x_right_half,debug1_half ] = step(bc,t_now,dt/2,Inf,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,ones(N,1),ones(N,1),ones(N,1),0,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    phi_half=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_half,cm_half,z_cp,z_cm, phi_x_right_half,voltage(t_now+dt/2) );
    [ cp_fine,cm_fine,phi_x_right_fine,debug1_fine ] = step(bc,t_now+dt/2,dt/2,dt/2,dx,N1,N2,N3,cp_half,cm_half,phi_half,phi_x_right_half,cp_now,cm_now,phi_now,phi_x_right_now,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    %evaluate error
    err=errcomp( 0,dt,0,cp_coarse,cm_coarse,phi_x_right_coarse,cp_fine,cm_fine,phi_x_right_fine);
    c=c+1;


end
% display(sprintf('err=%.3g,time=%.8f,dt=%.3g,c=%d',err,t_now,dt,c));
%perform richardson extrapolation
p=1; %the first step with dt_old=Inf is equivalent to coarse/fine forward euler, p=1
% cp_new = (2^p*cp_fine-cp_coarse)/(2^p-1);
% cm_new = (2^p*cm_fine-cm_coarse)/(2^p-1);
% phi_x_right_new=(2^p*phi_x_right_fine-phi_x_right_coarse)/(2^p-1); %phi_x_right continues to be 0 for voltage boundary conditions
cp_new = cp_coarse;
cm_new = cm_coarse;
phi_x_right_new=phi_x_right_coarse; %phi_x_right continues to be 0 for voltage boundary conditions

%solve for phi at this time level
phi_new = poisson1d( bc,dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_new,cm_new,z_cp,z_cm, phi_x_right_new,voltage(t_now+dt) );

%store values at fine half-step
cp_half_old=cp_half;
cm_half_old=cm_half;
phi_x_right_half_old=phi_x_right_half;
phi_half_old=phi_half;

%advance time
cp_old=cp_now;
cm_old=cm_now;
phi_x_right_old=phi_x_right_now;
phi_old=phi_now;
dt_old=dt;

cp_now=cp_new;
cm_now=cm_new;
phi_x_right_now=phi_x_right_new;
phi_now=phi_new;

% if t_now>=t_save*n;
%     display('snapshot taken');
    t(n+1)=t_now+dt;
    cp(:,n+1)=cp_new;
    cm(:,n+1)=cm_new;
    phi_x_right(n+1)=phi_x_right_new;
    phi(:,n+1)=phi_new;
    [j_f(n+1),j_d(n+1),j(n+1),v(n+1)]=postprocess(bc,dx,N1,N2,N3,t(n+1),t_now,cp(:,n+1),cm(:,n+1),phi(:,n+1),phi_old,phi_x_right(n+1),lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    err_save(n+1)=err;
    % save(file_name,'cp','cm','phi','j_f','j_d','j','v','t','x','err_save'); %print to file
    n=n+1;
% end

t_now=t_now+dt;
c=0; %reset convergence counter




%subsequent time steps

while t_now<t_end - eps
    %coarse step (with same dt as previous step)
    [ cp_coarse,cm_coarse,phi_x_right_coarse,debug1_coarse ] = step(bc,t_now,dt,dt_old,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_old,cm_old,phi_old,phi_x_right_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    %fine step (with same dt as previous step)
    [ cp_half,cm_half,phi_x_right_half,debug1_half ] = step(bc,t_now,dt/2,(dt_old)/2,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_half_old,cm_half_old,phi_half_old,phi_x_right_half_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    phi_half=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_half,cm_half,z_cp,z_cm, phi_x_right_half,voltage(t_now+dt/2) );
    [ cp_fine,cm_fine,phi_x_right_fine,debug1_fine ] = step(bc,t_now+dt/2,dt/2,dt/2,dx,N1,N2,N3,cp_half,cm_half,phi_half,phi_x_right_half,cp_now,cm_now,phi_now,phi_x_right_now,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
    %evaluate error
    err=errcomp( 1,dt,dt_old,cp_coarse,cm_coarse,phi_x_right_coarse,cp_fine,cm_fine,phi_x_right_fine);
    while abs(err - tol)>range 
        %take a stab at the new dt
        fac=(tol/err)^(.5);
        if fac > eta_max fac=eta_max;
        elseif fac < eta_min fac=eta_min;
        end
        dt=dt*fac;
        %If my calculated dt goes over dt_max, take a step with dt_max and exit
        %the loop.
        if dt>=dt_max
            dt=dt_max;
            [ cp_coarse,cm_coarse,phi_x_right_coarse,debug1_coarse ] = step(bc,t_now,dt,dt_old,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_old,cm_old,phi_old,phi_x_right_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
            [ cp_half,cm_half,phi_x_right_half,debug1_half ] = step(bc,t_now,dt/2,(dt_old)/2,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_half_old,cm_half_old,phi_half_old,phi_x_right_half_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
            phi_half=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_half,cm_half,z_cp,z_cm, phi_x_right_half,voltage(t_now+dt/2) );
            [ cp_fine,cm_fine,phi_x_right_fine,debug1_fine ] = step(bc,t_now+dt/2,dt/2,dt/2,dx,N1,N2,N3,cp_half,cm_half,phi_half,phi_x_right_half,cp_now,cm_now,phi_now,phi_x_right_now,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
            break;
        end
        %coarse step (with same dt as previous step)
        [ cp_coarse,cm_coarse,phi_x_right_coarse,debug1_coarse ] = step(bc,t_now,dt,dt_old,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_old,cm_old,phi_old,phi_x_right_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        %fine step (with same dt as previous step)
        [ cp_half,cm_half,phi_x_right_half,debug1_half ] = step(bc,t_now,dt/2,(dt_old)/2,dx,N1,N2,N3,cp_now,cm_now,phi_now,phi_x_right_now,cp_half_old,cm_half_old,phi_half_old,phi_x_right_half_old,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        phi_half=poisson1d(bc, dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_half,cm_half,z_cp,z_cm, phi_x_right_half,voltage(t_now+dt/2) );
        [ cp_fine,cm_fine,phi_x_right_fine,debug1_fine ] = step(bc,t_now+dt/2,dt/2,dt/2,dx,N1,N2,N3,cp_half,cm_half,phi_half,phi_x_right_half,cp_now,cm_now,phi_now,phi_x_right_now,lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        %evaluate error
        err=errcomp( 1,dt,dt_old,cp_coarse,cm_coarse,phi_x_right_coarse,cp_fine,cm_fine,phi_x_right_fine);
        c=c+1;
        %display(sprintf('n=%d,err=%.3g,time=%.3g,dt=%.3g,c=%d',n,err,t(n),dt,c));
        

    end
    % display(sprintf('err=%.3g,time=%.8f,dt=%.3g,c=%d',err,t_now,dt,c));
    %extrapolate using extrapolation formula
%     cp_new = extrapolate(cp_coarse,cp_fine, dt,dt_old);
%     cm_new = extrapolate(cm_coarse,cm_fine, dt,dt_old);
%     phi_x_right_new=extrapolate(phi_x_right_coarse,phi_x_right_fine, dt,dt_old); %phi_x_right continues to be 0 for voltage boundary conditions
    cp_new = cp_coarse;
    cm_new = cm_coarse;
    phi_x_right_new = phi_x_right_coarse;
    %solve for phi at this time level
    phi_new = poisson1d( bc,dx,N1,N2,N3, lambda_s, epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,cp_new,cm_new,z_cp,z_cm, phi_x_right_new,voltage(t_now+dt) );
    
    cp_half_old=cp_half;
    cm_half_old=cm_half;
    phi_x_right_half_old=phi_x_right_half;
    phi_half_old=phi_half;
    
    %advance time
    cp_old=cp_now;
    cm_old=cm_now;
    phi_x_right_old=phi_x_right_now;
    phi_old=phi_now;
    dt_old=dt;

    cp_now=cp_new;
    cm_now=cm_new;
    phi_x_right_now=phi_x_right_new;
    phi_now=phi_new;

%     if t_now>=t_save*n;
%         display('snapshot taken');
        t(n+1)=t_now+dt;
        cp(:,n+1)=cp_new;
        cm(:,n+1)=cm_new;
        phi_x_right(n+1)=phi_x_right_new;
        phi(:,n+1)=phi_new;
        [j_f(n+1),j_d(n+1),j(n+1),v(n+1)]=postprocess(bc,dx,N1,N2,N3,t(n+1),t_now,cp(:,n+1),cm(:,n+1),phi(:,n+1),phi_old,phi_x_right(n+1),lambda_s,D_0,D_p,D_m,kc_left,jr_left,kc_right,jr_right,epsilon_1,epsilon_2,epsilon_3,epsilon_s_left,epsilon_s_right,z_cp,z_cm);
        err_save(n+1)=err;
        % save(file_name,'cp','cm','phi','j_f','j_d','j','v','t','x','err_save'); %print to file
        n=n+1;
%    end
    
    t_now=t_now+dt;
    c=0; %reset convergence counter
end



runtime=toc(start);
fprintf('Runtime: %d minutes and %f seconds\n',floor(runtime/60),rem(runtime,60));
% plot(t(1:end-1),log10(diff(t)),'k-','Linewidth',1);
% xlabel('Time','FontSize',14);ylabel('$log_{10}(dt)$','Interpreter','LaTex','FontSize',14)

% xlabel('$x$','Interpreter','LaTex','FontSize',14);ylabel('$c_+$','Interpreter','LaTex','FontSize',14)
% xlabel('$x$','Interpreter','LaTex','FontSize',14);ylabel('$\phi$','Interpreter','LaTex','FontSize',14)
