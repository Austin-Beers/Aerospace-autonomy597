clc
clear
% close all hidden;

delta_t = .1;
end_time = 60;
l = 0:delta_t:end_time;

p = (pi/6);
q = cos(6*l/pi);
r = 3*sin(30*l/pi);
z = (30*pi/6) * ones(1,length(l));
% plots pqr with respect to t
figure(1)
plot (l,p);
hold on


plot (l,subs(z),l,subs(q),l,subs(r));

hold off
% initial conditions
x_new = zeros(3,1);
p_ned = zeros(3,1);
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
dcm = c_phi * c_theta * c_si
c_dot_dcm = c_deriv(dcm, t_var)
% solves C_dot at 0,0,0
% % [cphi, ctheta, csi] = deal(0)
% matlab didnt like when I tried to pass dcm_0 to rk4 alg. after assigning
% result of  0,0,0 to dcm so I hard coded the result below
% dcm_0 = subs(dcm)
dcm_0 = [1 0 0;
        0 1 0;
        0 0 1;]


% rk4 algorithm to numerically integrate between 0 and end time
for t = 0:delta_t:end_time
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
    v_ned = dcm_0 * v_0
    % integrate position to get velocity and add old postion to new position
    p_ned = p_ned + (delta_t * v_ned)

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
    end
end


%Plot for time rate of change of Euler angles and Euler angles
figure(2);
plot(l,theta_matrix_euler_rates);
hold on;
plot(l,phi_matrix_euler_rates);
plot(l,si_matrix_euler_rates);
hold off;
legend('theta_dot','phi_dot','psi_dot')
xlabel('time (sec)')
ylabel('angle (rad) per time (sec)')
figure(3);
plot(l,mod(theta_matrix_euler_angles-pi,2*pi)-pi);
hold on;
plot(l,phi_matrix_euler_angles);
plot(l,si_matrix_euler_angles);
hold off;
legend('theta','phi','psi')
xlabel('time (sec)')
ylabel('angle (rad)')

% Plot for time velocity and position
figure(4);
plot(l,velocity_matrix_n);
hold on;
plot(l,velocity_matrix_e);
plot(l,velocity_matrix_d);
hold off;
legend('v_n','v_e','v_d')
xlabel('time (sec)')
ylabel('velocity')
figure(5);
plot(l,position_matrix_n);
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

    p = pi / 6;
    q = cos(t * 6 / pi);
    r = 3 * sin(t * 30 / pi);

    % returns cdot using strapdown equation
    cdot = c * [0 -r q; 
               r 0 -p; 
               -q p 0;];
end


