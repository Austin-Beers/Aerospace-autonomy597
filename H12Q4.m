%% Aerospace Autonomy 497/597
% Homework 12

% Housekeeping tools
clear all; close all hidden; clc;


%% SETTING UP OUR SIMULATION OF REALITY USING NON-LINEAR MODEL
DT = 0.1;          % size of timestep in numerical integration (seconds)
time = 0;          % initializing the simulation time to 0 seconds
ENDTIME = 150;     % simulation ENDTIME in seconds
i = 1;             % counter for indexing arrays in relation to 'time'

% parameterize the statistics of our measures
sigmaPositionSensor = 150; % feet; 2 of these sensors, one to measure x and one to measure y
sigmaVelocitySensor = 5;   % knots
sigmaHeadingSensor = 5;    % degrees

%% SIMULATION
% set the time between measurements
mInterval = 10; % how often do we get sensor measurements? in seconds
mCntr     = 1; % index for storing the captured sensor measurement

% initialize our time history of state, with state vector containing:
% [PosX(East)_feet; PosY(North)_feet; V_knots; heading_degrees]
x(:, i) = [-750 0 60 0]';

% initialize our time history of state, with state vector containing:
% [Positionx(East)_feet; Positiony(North)_feet; V_knots; heading_degrees]
% Use a column vector to describe the state
% Notation dictates state information is stored in index "i" to allow use
% to capture "time" history using the second index
x(:, i) = [0 0 60 0]';

% initialize our controls
% We are "COMMANDING":  [x-position (feet), velocity (kts), heading (deg)]
u = [750 60 0]';

% initialize a TIME array
time_arr(1) = 0;
% and a MEASUREMENT TIME array
m_time_arr(1) = 0;

% controller GAINS
gains = [0.011, 0.0667, 1];
Kpsibyx = gains(1);
Kpsi = gains(2);
Kv = gains(3);


% The DISTURBANCE GAIN matrix
G = [0 0; 0 0; 1 0; 0 1];
% The Process Noise Covariance matrix
% represents lack of knowledge in teh model's ability to estimate states
Q = zeros(2);
% The Measurement Noise Covariance matrix
% quantifies uncertainty in the emasuremens used to update the states
R = zeros(4);
R(1, 1) = sigmaPositionSensor * sigmaPositionSensor;
R(2, 2) = sigmaPositionSensor * sigmaPositionSensor;
R(3, 3) = sigmaVelocitySensor * sigmaVelocitySensor;
R(4, 4) = sigmaHeadingSensor * sigmaHeadingSensor;
% modeling the plant in linear form to find the Kalman filter gain
% note conversion of KNOTS to FPS in elements multiplied by V
A = [0 0 0 u(2) * 6076/3600; 0 0 6076/3600 0; 0 0 -Kv 0; -Kpsibyx 0 0 -Kpsi];
B = [0 0 0; 0 0 0; 0 Kv 0; Kpsibyx 0 Kpsi];
C = eye(4);
D = [0 0 0; 0 0 0; 0 0 0; 0 0 0];
% Construct state-space model using MATLAB's built-in functionality
% We "need" this to get a REALLY good estimate for Pinitial and the
% Kalman Filter Gain "K" -- which we will promptly ignore and calculate on
% our own anyways!
plant = ss(A, [B G], C, [D zeros(4, 2)], 'inputname', {'Vc' 'Psic' 'Xc' 'WonV' 'WonPsi'}, ...
'outputname', {'MeasureX' 'MeasureY' 'MeasureV' 'MeasurePsi'});
[kalmf, K, Pinitial] = kalman(plant, Q, R);
% This is our ESTIMATE of the state vector. We take it to be the same at
% the initial conditions and will then use the simulation to deviate from
% the navigation solution due to noise and sensor measurement
xest(:, i) = [0 0 60 0]';
% actual error squared! Someday it'll be interesting to know how good the
% Kalman Filter is doing, and we could plot this matrix through time
% (using the 3rd "i" index) and compare that to our estimate for P via the
% Kalman Filter update rules
% Pactual(:, :, i) = [(x(1, i) - xest(1, i)) * (x(1, i) - xest(1, i)) 0 0 0; ...
%     0 (x(2, i) - xest(2, i)) * (x(2, i) - xest(2, i)) 0 0; ...
%     0 0 (x(3, i) - xest(3, i)) * (x(3, i) - xest(3, i)) 0; ...
%     0 0 0 (x(4, i) - xest(4, i)) * (x(4, i) - xest(4, i))];
% Here we use the best guess for P from MATLAB's internal routines
% We "could" just estimate this as the identiy matrix and hope that the KF
% resolves it for us eventually
P(:, :, i) = Pinitial;


model = 'data'; % use 'simulation' or 'data' to assess the dynamics


% enter the simulation loop
while time < ENDTIME
    % while the current simulation time is less than the final end time
    Kfactor = zeros(4,1);
    Pfactor = eye(4);
    %% Should we take a measurement with our sensors?
    % IF A MEASUREMENT COMES IN, SIMULATING THE MEASUREMENT PROCESS
    % AND USING MEASUREMENT TO CORRECT KALMAN FILTER
    if mod(time, mInterval) < DT
        % IF it has been `mInterval` (i.e., our sensor's sampling rate),
        % we enter the "get a sensor measurement" part of our code

        %% SINGLE SENSING
        % make a measurement with any sensor
        % sensors are described with normal distributions with zero bias
        % and given standard deviations.
        % MATLAB provides the randn function to sample a zero-mean and
        % standard deviation of 1
        y(1, mCntr) = x(1, i) + randn * sigmaPositionSensor;
        y(2, mCntr) = x(2, i) + randn * sigmaPositionSensor;
        y(3, mCntr) = x(3, i) + randn * sigmaVelocitySensor;
        y(4, mCntr) = x(4, i) + randn * sigmaHeadingSensor;

        

        % Here is where we correct the KALMAN FILTER, would be running onboard
        % Make our estimate of what the measures would be, from xest
        yest = C * xest(:, i);
        % calculate the Kalman Filter Gain, K
        % I am using "right-hand division" so that MATLAB figures out the best
        % algorithm for inverting the matrix
        K = P(:, :, i) * C' / (C * P(:, :, i) * C' + R);
        ydiff = (y(:, mCntr) - yest);
        Kfactor = K * ydiff;
        Pfactor = (eye(size(A)) - K * C);
        
        % store the MEASUREMENT TIME interval
        m_time_arr(mCntr) = time;
        
        % INCREMENT our measurement counter index so we can store all of
        % our sensed information into nice MATLAB matrices for plotting
        mCntr = mCntr + 1;
    end

    %% Perform Runge-Kutta 4th Order integration estimates
    % Evaluate the DFUNC at the current state, forcing, and time

    % in this case, x varies throughout the interval, but the elements of
    % the control vector 'u' are constant
    % AND, THIS IS NEW: we pass it y to use instead of raw measures
    xdot1 = HW12derivative(x(:, i), xest(:, i), u, time, gains, model);
    xdot2 = HW12derivative(x(:, i) + xdot1 * DT / 2, xest(:, i), u, time + DT / 2, gains, model);
    xdot3 = HW12derivative(x(:, i) + xdot2 * DT / 2, xest(:, i), u, time + DT / 2, gains, model);
    xdot4 = HW12derivative(x(:, i) + xdot3 * DT, xest(:, i), u, time + DT, gains, model);

    % "Average" the derivative estimates from RK4 to get a best approx
    % Store this in a variable that will contain the derivative of the
    % state variable "x", for a "time" indicator "i"
    totalxdot(:, i) = (xdot1 + 2 * xdot2 + 2 * xdot3 + xdot4) / 6;
%     totalxdot(:, i) = totalxdot(:, i) + Kfactor;
    % update state by doing the actual integration
    % and store that in the next time step "i+1"
    x(:, i + 1) = x(:, i) + totalxdot(:, i) * DT;
    

    % VARIANT UISNG THE NONLINEAR MODEL
    xestdot1 = HW12derivative(xest(:, i), y(:, mCntr - 1), u, time, gains, model);
    xestdot2 = HW12derivative(xest(:, i) + xestdot1 * DT / 2, y(:, mCntr - 1), u, time + DT / 2, gains, model);
    xestdot3 = HW12derivative(xest(:, i) + xestdot2 * DT / 2, y(:, mCntr - 1), u, time + DT / 2, gains, model);
    xestdot4 = HW12derivative(xest(:, i) + xestdot3 * DT, y(:, mCntr - 1), u, time + DT, gains, model);
    totalxestdot(:, i) = (xestdot1 + 2 * xestdot2 + 2 * xestdot3 + xestdot4) / 6;
    totalxestdot(:, i) = totalxdot(:, i) + Kfactor;
    xest(:, i + 1) = xest(:, i) + totalxestdot(:, i) * DT;



    Pdot1 = A * (P(:, :, i)) * A' + G * Q * G';
    Pdot2 = A * (P(:, :, i) + Pdot1 * DT / 2) * A' + G * Q * G';
    Pdot3 = A * (P(:, :, i) + Pdot2 * DT / 2) * A' + G * Q * G';
    Pdot4 = A * (P(:, :, i) + Pdot3 * DT) * A' + G * Q * G';
    totalPdot = (Pdot1 + 2 * Pdot2 + 2 * Pdot3 + Pdot4) / 6;
    P(:, :, i + 1) = P(:, :, i) + totalPdot * DT;
    P(:, :, i + 1) = Pfactor * P(:, :, i + 1);


    i = i + 1;        % increment our iteration counter
    time = time + DT; % increment time
    time_arr(i) = time;
end

%% Plot the RESULTS
figure(100); clf;
subplot(4, 1, [1, 2])
plot(x(1, :), x(2, :), y(1, :), y(2, :))
legend('Actual', 'Measured')
title('$y$ by $x$ position, feet', 'interpreter', 'latex');

subplot(4, 1, 3)
plot(time_arr, x(3, :), m_time_arr, y(3, :))
legend('Actual', 'Measured')
title('Timeline of $V$, knots', 'interpreter', 'latex');

subplot(4, 1, 4)
plot(time_arr, x(4, :), m_time_arr, y(4, :))
legend('Actual', 'Measured')
title('Timeline of $\psi$, degrees', 'interpreter', 'latex');
