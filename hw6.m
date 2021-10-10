clc
clear
close all hidden

q ="part-3" % sets homework part via conditional in rk4 options are "part-1", "part-2", or "part-3"
% initial conditions
x_0 = [0;
       0;
       60; 
       0;];
u = [0 750 750;
     60 60 60;
     0 0 180;];

% Setting up time for plots
delta_t = 0.1;
end_time = 60;
l = 0:delta_t:end_time; 

state_vector = deal(zeros((end_time/delta_t)+1,4));


K_si_x_org = 0.011;
k_si_x_update = 0;

% rk4 algorithm to numerically integrate between 0 and end time
for t = 0:delta_t:end_time

    %%%RUNGE KUTTA%%%
    limit = false; % limit values of si_dot 
    local_u = [];
    local_k_si_x = [];
    if q == "part-1" 
        local_k_si_x = K_si_x_org;
        local_u = u(:, 1);
    end
    if q == "part-2"
        local_k_si_x = K_si_x_org;
        local_u = u(:, 2);
    end
    if q == "part-3"
        local_k_si_x = k_si_x_update;
        local_u = u(:, 3);
    end

    % update global x_0 to local variable 
    x_old = x_0;
    xdot_1 = deriv_func(x_old,local_u, local_k_si_x,limit );
    xdot_2 = deriv_func(x_old + xdot_1 * (delta_t/2),local_u,local_k_si_x, limit);
    xdot_3 = deriv_func(x_old + xdot_2 * (delta_t/2),local_u, local_k_si_x, limit);
    xdot_4 = deriv_func(x_old + xdot_3 * delta_t,local_u, local_k_si_x, limit);
    x_dot_RK4 = (1/6) * (xdot_1 + 2*xdot_2 +2*xdot_3 + xdot_4);
    x_0 = x_old + (delta_t*x_dot_RK4);

    if(t<=end_time)
        state_vector(uint32((t/delta_t)+1), :) = x_0;
    end
end



%%%PLOTS FROM RUNGE KUTTA
% % Plot for time rate of change of Euler angles and Euler angles
transpose_state = state_vector.';
figure(1);
plot(l,transpose_state);
title("State vector vs Time " + q);
xlabel("Time (sec)");
ylabel("State vector");
legend("East (deg)","North (deg)","Velocity (fps)","Heading (fps)");
figure(2);
plot(transpose_state(2,:), transpose_state(1,:));
title("North vs East " + q);
xlabel("East (deg)")
ylabel("North (deg)")







