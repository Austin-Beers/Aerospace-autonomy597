
clc
clear
% close all hidden;

delta_t = .01;
end_time = 150;
l = 0:delta_t:end_time;

p = 0.04 * sin(pi * l / 15);
q = 0.05 * cos(pi * l / 6);
r = pi / 60;
z = (pi/60) * ones(1,length(l));
% plots pqr with respect to t
figure(1)
plot (l,p);
hold on


plot (l,subs(z),l,subs(q),l,subs(r));

hold off
% initial conditions
x_new = zeros(3,1);
p_ned = zeros(3,1);
x_fe_new = zeros(3,1);
p_fe_ned = zeros(3,1);
x_be_val = zeros(3,1);
v_0 = [30.8667;
      0;
      0;];

[theta_matrix_euler_rates, ...
    phi_matrix_euler_rates, ...
    si_matrix_euler_rates, ...
    theta_matrix_euler_angles, ...
    phi_matrix_euler_angles, ...
    si_matrix_euler_angles, ...
    velocity_matrix_n, ...
    velocity_matrix_e, ...
    velocity_matrix_d, ...
    position_matrix_n, ...
    position_matrix_e, ...
    position_matrix_d, ...
    ...
    tmer_fe, ...
    pmer_fe, ...
    smer_fe, ...
    tmea_fe, ...
    pmea_fe, ...
    smea_fe, ...
    vmn_fe, ...
    vme_fe, ...
    vmd_fe, ...
    pmn_fe, ...
    pme_fe, ...
    pmd_fe, ...
    ] = deal(zeros(1, (end_time/delta_t)+1));
disp(length(theta_matrix_euler_angles))


% Calculate DCM at (0,0,0)
syms cphi ctheta csi t_var
c_phi = [1 0 0;
        0 cos(cphi) sin(cphi);
        0 -sin(cphi) cos(cphi);];
c_theta = [cos(ctheta) 0 -sin(ctheta);
          0 1 0;
          sin(ctheta) 0 cos(ctheta);];
c_si = [cos(csi) sin(csi) 0;
       -sin(csi) cos(csi) 0;
       0 0 1;];
% finds dcm
dcm = c_phi * c_theta * c_si;
c_dot_dcm = c_deriv(dcm, t_var);
% solves C_dot at 0,0,0
% % [cphi, ctheta, csi] = deal(0)
% matlab didnt like when I tried to pass dcm_0 to rk4 alg. after assigning
% result of  0,0,0 to dcm so I hard coded the result below
% dcm_0 = subs(dcm)
dcm_0 = [1 0 0;
        0 1 0;
        0 0 1;];

dcm_fe_0 = [1 0 0;
        0 1 0;
        0 0 1;];


% rk4 algorithm to numerically integrate between 0 and end time
for t = 0:delta_t:end_time

    %%%RUNGE KUTTA%%%

    % update global x_new to local variable 
    x_old = x_new;
    xdot_1 = deriv(x_old,t);
    xdot_2 = deriv(x_old + xdot_1 * (delta_t/2),t + (delta_t/2));
    xdot_3 = deriv(x_old + xdot_2 * (delta_t/2),t + (delta_t/2));
    xdot_4 = deriv(x_old + xdot_3 * delta_t,t + delta_t);
    x_dot_RK4 = (1/6) * (xdot_1 + 2*xdot_2 +2*xdot_3 + xdot_4);
    x_new = x_old + (delta_t*x_dot_RK4);
    
    % pass dcm_0 into c_deriv function to update transformation
    cdot_1 = c_deriv(dcm_0,t);
    cdot_2 = c_deriv(dcm_0 + cdot_1 * (delta_t/2),t + (delta_t/2));
    cdot_3 = c_deriv(dcm_0 + cdot_2 * (delta_t/2),t + (delta_t/2));
    cdot_4 = c_deriv(dcm_0 + cdot_3 * delta_t,t + delta_t);
    c_dot_RK4 = (1/6) * (cdot_1 + 2*cdot_2 + 2*cdot_3 + cdot_4);
    % update dcm from result of rk4 algorithm
    dcm_0 = dcm_0 + (delta_t*c_dot_RK4);
    % transform velocity of body frame to local ned frame
    v_ned = dcm_0 * v_0;
    % integrate position to get velocity and add old postion to new position
    p_ned = p_ned + (delta_t * v_ned);

    %%%FORWARD EULER%%%

    % gimbal equation
    x_fe_old = x_fe_new;
    x_fe_dot = x_fe_old + (delta_t * deriv(x_fe_old, t));
    x_fe_new = x_fe_old + (delta_t * deriv(x_fe_dot, t));

    % strapdown equation
    c_fe_dot = dcm_fe_0 + (delta_t * c_deriv(dcm_fe_0, t));
    dcm_fe_0 = dcm_fe_0 + (delta_t * c_deriv(c_fe_dot, t));
    v_fe_ned = dcm_fe_0 * v_0;
    p_fe_ned = p_fe_ned + (delta_t * v_fe_ned);

    %%%BACKWARD EULER%%%
%     options = optimset('TolX', 1e-06);
%     x_be_old = x_be_val(:, end);
%     der = @deriv;
%     anon_x = @(x) x - delta_t*der(x,t+delta_t) - x_be_old;
%     new_x = fsolve(anon_x, x_be_old, options)
%     x_be_val = [x_be_val new_x]

    
    if(t<=end_time)
        theta_matrix_euler_rates(uint32((t/delta_t)+1)) = x_dot_RK4 (1);
        phi_matrix_euler_rates(uint32((t/delta_t)+1)) = x_dot_RK4 (2);
        si_matrix_euler_rates(uint32((t/delta_t)+1)) = x_dot_RK4 (3);

        theta_matrix_euler_angles(uint32((t/delta_t)+1)) = x_new (1);
        phi_matrix_euler_angles(uint32((t/delta_t)+1)) = x_new (2);
        si_matrix_euler_angles(uint32((t/delta_t)+1)) = x_new (3);

        velocity_matrix_n(uint32((t/delta_t)+1)) = v_ned (1);
        velocity_matrix_e(uint32((t/delta_t)+1)) = v_ned (2);
        velocity_matrix_d(uint32((t/delta_t)+1)) = v_ned (3);

        position_matrix_n(uint32((t/delta_t)+1)) = p_ned (1);
        position_matrix_e(uint32((t/delta_t)+1)) = p_ned (2);
        position_matrix_d(uint32((t/delta_t)+1)) = p_ned (3);

        tmer_fe(uint32((t/delta_t)+1)) = x_fe_dot (1);
        pmer_fe(uint32((t/delta_t)+1)) = x_fe_dot (2);
        smer_fe(uint32((t/delta_t)+1)) = x_fe_dot (3);
        

        tmea_fe(uint32((t/delta_t)+1)) = x_fe_new (1);
        pmea_fe(uint32((t/delta_t)+1)) = x_fe_new (2);
        smea_fe(uint32((t/delta_t)+1)) = x_fe_new (3);

        vmn_fe(uint32((t/delta_t)+1)) = v_fe_ned (1);
        vme_fe(uint32((t/delta_t)+1)) = v_fe_ned (2);
        vmd_fe(uint32((t/delta_t)+1)) = v_fe_ned (3);

        pmn_fe(uint32((t/delta_t)+1)) = p_fe_ned (1);
        pme_fe(uint32((t/delta_t)+1)) = p_fe_ned (2);
        pmd_fe(uint32((t/delta_t)+1)) = p_fe_ned (3);
    end
end
display(smer_fe)


%PLOTS FOR FORWARD EULER
figure(2);
plot(l,mod(tmer_fe-pi,2*pi)-pi);

title('Forward euler Euler rates');
hold on;
plot(l,pmer_fe);
plot(l,smer_fe);
hold off;
legend('theta_dot','phi_dot','psi_dot')
xlabel('time (sec)')
ylabel('angle (rad) per time (sec)')
figure(3);
plot(l,mod(tmea_fe-pi,2*pi)-pi);
title('Forward euler Euler angles');
hold on;
plot(l,pmea_fe);
plot(l,smer_fe);
hold off;
legend('theta','phi','psi')
xlabel('time (sec)')
ylabel('angle (rad)')

% Plot for time velocity and position using forward Euler
figure(4);
plot(l,vmn_fe);
title('Forward euler velocity in ned');
hold on;
plot(l,vme_fe);
plot(l,vmd_fe);
hold off;
legend('v_n','v_e','v_d')
xlabel('time (sec)')
ylabel('velocity')
figure(5);
plot(l,pmn_fe);
title('Forward euler position in ned');
hold on;
plot(l,pme_fe);
plot(l,pmd_fe);
hold off;
legend('p_n','p_e','p_d')
xlabel('time (sec)')
ylabel('position')

%%%PLOTS FROM RUNGE KUTTA
% Plot for time rate of change of Euler angles and Euler angles
figure(6);
plot(l,theta_matrix_euler_rates);
title('Runge Kutta Euler Rates');
hold on;
plot(l,phi_matrix_euler_rates);
plot(l,si_matrix_euler_rates);
hold off;
legend('theta_dot','phi_dot','psi_dot')
xlabel('time (sec)')
ylabel('angle (rad) per time (sec)')
figure(7);
plot(l,mod(theta_matrix_euler_angles-pi,2*pi)-pi);
title('Runge Kutta Euler Angles')
hold on;
plot(l,phi_matrix_euler_angles);
plot(l,si_matrix_euler_angles);
hold off;
legend('theta','phi','psi')
xlabel('time (sec)')
ylabel('angle (rad)')

% Plot for time velocity and position
figure(8);
plot(l,velocity_matrix_n);
title('Runge Kutta velocity ned');
hold on;
plot(l,velocity_matrix_e);
plot(l,velocity_matrix_d);
hold off;
legend('v_n','v_e','v_d')
xlabel('time (sec)')
ylabel('velocity')
figure(9);
plot(l,position_matrix_n);
title('Runge Kutta position ned');
hold on;
plot(l,position_matrix_e);
plot(l,position_matrix_d);
hold off;
legend('p_n','p_e','p_d')
xlabel('time (sec)')
ylabel('position')



function xdot = deriv(x,t)
    phi = x(1);
    theta = x(2);
    
    A = [1 tan(theta)*sin(phi) tan(theta)*cos(phi);
        0 cos(phi) -sin(phi)
        0 sec(theta)*sin(phi) sec(theta)*cos(phi)];
    B = [(pi/6);
        cos(6*t/pi) ;
        3*sin(30*t/pi)];
    % returns euler rates using gimbal equation
    xdot = A * B;
end

function cdot = c_deriv(c,t)

    p = 0.04 * sin(pi * t / 15);
    q = 0.05 * cos(pi * t / 6);
    r = pi / 60;

    % returns cdot using strapdown equation
    cdot = c * [0 -r q; 
               r 0 -p; 
               -q p 0;];
end

