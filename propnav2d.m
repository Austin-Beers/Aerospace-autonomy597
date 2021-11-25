%% Aerospace Autonomy
% Homework 11
% proportionalnavigation2D performs intercept / rendezvous in 2D space

% Housekeeping tools
clear all; close all hidden; clc;

% Simulation Setup parameters
n = 1;                   % index counter for time steps
RANGETHRESHOLD = 50;     % distance which to consider "intercepted"
timestep = 0.01;         % size of timestep in numerical integration

MODEL_TYPE = 0;
% 0 -- Intercept
% 1 -- Rendezvous


% TARGET DYNAMICS
targetposition(n, :) = [2000, 2000]; % position in feet
targetvelocity(n) = 250;             % magnitude of velocity in fps
targetheading(n) = 0;                % heading in radians

% target changes speed and heading via randomly selected jumps
% these values scale how "big" those jumps are
TARGETVELOCITYCHANGE = 5 * (10 * timestep);            % fps
TARGETHEADINGCHANGE = 5 * pi / 180 * (10 * timestep);  % rad/second

% INTERCEPTOR DYNAMICS
interceptorposition(n, :) = [0, 0];  % position in feet
interceptorvelocity(n) = 350;        % magnitude of velocity in fps
interceptorheading(n) = 0;           % heading in radians
interceptorsideacceleration(n) = 0;  % acceleration normal to current heading, in fpsps

% simplified notation to make the math more readable
% XDIFF is just the "delta x" for each time step of line-of-sight
xdiff = targetposition(n, 1) - interceptorposition(n, 1);
% YDIFF is just the "delta y" for each time step of line-of-sight
ydiff = targetposition(n, 2) - interceptorposition(n, 2)

% RANGE is distance between the TARGET and the INTERCEPTOR
range(n) = sqrt(xdiff * xdiff + ydiff * ydiff);

time(n) = 0;

while (range(n) > RANGETHRESHOLD)
    % while the distance between the target and the interceptor is large
    % CALCULATE the angle for the line-of-sight
    sigma(n) = atan2(xdiff, ydiff) - interceptorheading(n); %3.1-11-16
    % CALCULATE the derivative of angle of line-of-sight using chain and quotient rules
    sigmadot(n) = calcSigmaDot( ...
        xdiff, ...
        ydiff, ...
        [targetvelocity(n), interceptorvelocity(n)], ...
        [targetheading(n), interceptorheading(n)] ...
        ); %3.2-11-19

    % here's where you need to figure out the interceptor's side acceleration to command for an intercept
    if MODEL_TYPE == 0
        interceptorsideacceleration(n) = -3 * interceptorvelocity(n) * sigmadot(n); %3.2-11-19
    elseif MODEL_TYPE == 1
%         interceptorsideacceleration(n) = calcAccRen( ...
%             xdiff, ...
%             ydiff, ...
%             [targetvelocity(n), interceptorvelocity(n)], ...
%             sigma, ...
%             [targetheading(n), interceptorheading(n)], ...
%             sigmadot(n) ...
%             ); %3.3-11-19
       interceptorsideacceleration(n) = -1 * interceptorvelocity(n) * (4 * sigmadot(n) + (2 * sigma(n))/(range(n)/ interceptorvelocity(n)));
    end

    % INCREMENT the simulation counter
    n = n + 1;

    % UPDATE the TARGET velocity and heading based on its dynamics
    % target changes speed via randomly selected jumps and is turning at a constant rate
    % RAND is a uniform value from [0,1], so we subtract 0.5 and multiply by 2 to make the 
    %   velocity change somewhere between [-1,1]*TARGETVELOCITYCHANGE
    targetvelocity(n) = targetvelocity(n - 1) +(rand(1) - .5) * 2 * TARGETVELOCITYCHANGE;
    % HEADING is just a constant turn
    targetheading(n) = targetheading(n - 1) + TARGETHEADINGCHANGE;
    % UNWRAP the heading of the TARGET
    if (targetheading(n) > 2 * pi) targetheading(n) = targetheading(n) - 2 * pi; end
    % INTEGRATE position using Euler method
    targetposition(n, 1) = targetposition(n - 1, 1) + targetvelocity(n) * sin(targetheading(n)) * timestep;
    targetposition(n, 2) = targetposition(n - 1, 2) + targetvelocity(n) * cos(targetheading(n)) * timestep;

    % INTERCEPTOR maintains a constant VELOCITY
    % interceptor's side acceleration changes the magnitude of its velocity and its heading
    interceptorvelocity(n) = interceptorvelocity(n - 1);
    % HEADING changes based on velocity steering (propotional gain)
    interceptorheading(n) = interceptorheading(n - 1) - atan2((interceptorsideacceleration(n - 1) * timestep), interceptorvelocity(n - 1));
    % UNWRAP the heading of the INTERCEPTOR
    if (interceptorheading(n) > 2 * pi) interceptorheading(n) = interceptorheading(n) - 2 * pi; end
    % INTEGRATE position using Euler method
    interceptorposition(n, 1) = interceptorposition(n - 1, 1) + interceptorvelocity(n) * sin(interceptorheading(n)) * timestep;
    interceptorposition(n, 2) = interceptorposition(n - 1, 2) + interceptorvelocity(n) * cos(interceptorheading(n)) * timestep;

    % CALCULATE the line-of-sight
    % simplified notation to make the math more readable
    xdiff = targetposition(n, 1) - interceptorposition(n, 1);
    ydiff = targetposition(n, 2) - interceptorposition(n, 2);
    range(n) = sqrt(xdiff * xdiff + ydiff * ydiff);

    % INCREMENT time in our simulation
    time(n) = time(n - 1) + timestep;
    % safety value in case you don't converge...
    if time(n) > 100
        break
    end

    % now to make the animation 'go'
    pause(timestep);
    plot(interceptorposition(1:n - 1, 1), interceptorposition(1:n - 1, 2), 'g', targetposition(1:n - 1, 1), targetposition(1:n - 1, 2), 'r', ...
        interceptorposition(n - 1, 1), interceptorposition(n - 1, 2), 'g+', targetposition(n - 1, 1), targetposition(n - 1, 2), 'or')
end

% print out some simulation information
fprintf('Simulation took %d steps and %.2f seconds to intercept', n, time(end))

%% MAKING PLOTS
figure
subplot(3, 1, 1)
plot(targetposition(:, 1), targetposition(:, 2), interceptorposition(:, 1), interceptorposition(:, 2))
legend('target', 'interceptor')
title('Target and interceptor positions')

subplot(3, 1, 2)
plot(time, targetheading * 180 / pi, time, interceptorheading * 180 / pi, time(1:n - 1), sigma * 180 / pi)
legend('target', 'interceptor', 'angle of line of sight')
title('Target and interceptor headings, and angle of line of sight, as a function of time')

subplot(3, 1, 3)
plot(time, range)
title('Range as a function of time')

function sigmaDot = calcSigmaDot(dx, dy, V, SI)
    arguments
        dx % xdiff (x position of target - x position of interceptor)
        dy % ydiff (y position of target - y position of interceptor)
        V % velocities relevant target,V(1), and interceptor, V(2)
        SI % headings relevant target,SI(1), and interceptor, SI(2)
    end

    dxdt = V(1) * sin(SI(1)) - V(2) * sin(SI(2));
    dydt = V(1) * cos(SI(1)) - V(2) * cos(SI(2));

    sigmaDot = (1 / (1 + (dx/dy)^2)) * ((dy * dxdt - dx * dydt) / dy^2);
end

function acc_ren = calcAccRen(dx, dy, V, SIGMA, SI, SIGMADOT)
    arguments
        dx % xdiff (x position of target - x position of interceptor)
        dy % ydiff (y position of target - y position of interceptor)
        V % velocities relevant target,V(1), and interceptor, V(2)
        SIGMA % line of sight (LOS) angle
        SI % headings relevant target,SI(1), and interceptor, SI(2)
        SIGMADOT % rate of change of LOS angle
    end
    R = sqrt(dx^2 + dy^2) * asin2((pi/2) - SIGMA);
    dxdt = V(1) * sin(SI(1)) - V(2) * sin(SI(2));
    dydt = V(1) * cos(SI(1)) - V(2) * cos(SI(2));
    Va = sqrt(dxdt^2 + dydt^2);
    
    acc_ren = -V * ( 4 * SIGMADOT + ((2 * SIGMA) / (R / Va)));
end
