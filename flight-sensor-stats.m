%% Aerospace Autonomy 497/597
% Homework 7, Question 1
clear all; close all hidden; clc;

%% SETTING UP OUR SIMULATION OF REALITY USING NON-LINEAR MODEL
dt = 0.01;
time = 0;
endtime = 150;
i = 1;

% parameterize the statistics of our measures
sigma_SensorGood = 10; % ft
sigma_SensorCheap = 100; % ft


%% THEORETICAL ESTIMATES
% theoretical calculations of combinations' statistics
% CHEAP SENSORS
% - Using  2 of them, what is the MEAN & STD?
% - Using  5 of them, what is the MEAN & STD?
% - Using 10 of them, what is the MEAN & STD?

% GOOD SENSORS
% - Using  2 of them, what is the MEAN & STD?
% - Using  5 of them, what is the MEAN & STD?
% - Using 10 of them, what is the MEAN & STD?

% COMBO SENSORS
% - Using 10 CHEAP and 1 GOOD, what is the MEAN & STD?
% - Using 10 CHEAP and 2 GOOD, what is the MEAN & STD?

%% SIMULATION
% set the time between measurements
mInterval = 1;
mCntr = 1;

% initialize our time history of state, with state vector containing:
% [PosX(East)_feet; PosY(North)_feet; V_knots; heading_degrees]
x(:, i) = [0 0 60 0]';

% initialize our controls
u = [750 60 0]';

% initialize a TIME array
time_arr(1) = 0;

% initialize weight with user input
ng = input("Choose number of good sensors: valid inputs are 1, 2, 5, or 10 \n : ");
if not(any(ismember([1,2,5,10], ng)))
    disp("please enter valid input for the number of good sensors: the program will now be restarted")
    return
end

nc = input("Choose number of cheap sensors: valid inputs are 1, 2, 5, or 10 \n : ");
if not(any(ismember([1,2,5,10], nc)))
    disp("please enter valid input for the number of cheap sensors: the program will now be restarted")
    return
end

dispPlots = input("Would you like to display the plots associated with the sensors: valid inputs are 0 for no and 1 for yes \n (0 ~ no) / (1 ~ yes): ");
if not(any(ismember([0,1], dispPlots))) 
    disp("please enter valid input for displaying plots: the program will now be restarted")
    return
end

%% Utilizes weight function which takes input from user on number of good and cheap sensors and then calulates the weigth for each based on equation in notes using standard deviation 
wg = weight(false, nc, ng,sigma_SensorCheap, sigma_SensorGood );
wc = weight(true, nc, ng,sigma_SensorCheap, sigma_SensorGood );

lastXMeasuresC = input("Choose last values to measure average over for cheap sensor: valid inputs are 2, 5, or 10 \n : ");
if not(any(ismember([2,5,10], lastXMeasuresC)))
    disp("please enter valid input for the last X measures of cheap sensor to average over: the program will now be restarted")
    return
end

lastXMeasuresG = input("Choose last values to measure average over for good sensor: valid inputs are 2, 5, or 10 \n : ");
if not(any(ismember([2,5,10], lastXMeasuresG)))
    disp("please enter valid input for the last X measures of good sensor to average over: the program will now be restarted")
    return
end

while time < endtime

    % calculate the 4 estimates of derivative over the time interval from t to t+dt

    % in this case, x varies throughout the interval, but the elements of
    % the control vector 'u' are constant
    xdot1 = HW7derivative(x(:, i), u, time);
    xdot2 = HW7derivative(x(:, i) + xdot1 * dt / 2, u, time + dt / 2);
    xdot3 = HW7derivative(x(:, i) + xdot2 * dt / 2, u, time + dt / 2);
    xdot4 = HW7derivative(x(:, i) + xdot3 * dt, u, time + dt);

    % calculate our 4th order Runge Kutta estimate of derivative
    totalxdot(:, i) = (xdot1 + 2 * xdot2 + 2 * xdot3 + xdot4) / 6;

    %update state
    x(:, i + 1) = x(:, i) + totalxdot(:, i) * dt;

    if mod(time, mInterval) < dt
        %% Two Good 
        mNoiseSensorGood(:, mCntr) = randn(20, 1) * sigma_SensorGood;
        mNoiseSensorCheap(:, mCntr) = randn(20, 1) * sigma_SensorCheap;
        if ng == 1 
            ySensorGood(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorGood(1, mCntr);
            ySensorGood(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorGood(2, mCntr);
        end
        if ng == 2
            ySensorGood(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorGood(1, mCntr);
            ySensorGood(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorGood(2, mCntr);
            ySensorGood(1, mCntr, 2) = x(1, i + 1) + mNoiseSensorGood(3, mCntr);
            ySensorGood(2, mCntr, 2) = x(2, i + 1) + mNoiseSensorGood(4, mCntr);
        end
        %% Five Good
        if ng == 5 
            ySensorGood(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorGood(1, mCntr);
            ySensorGood(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorGood(2, mCntr);
            ySensorGood(1, mCntr, 2) = x(1, i + 1) + mNoiseSensorGood(3, mCntr);
            ySensorGood(2, mCntr, 2) = x(2, i + 1) + mNoiseSensorGood(4, mCntr);
            ySensorGood(1, mCntr, 3) = x(1, i + 1) + mNoiseSensorGood(5, mCntr);
            ySensorGood(2, mCntr, 3) = x(2, i + 1) + mNoiseSensorGood(6, mCntr);
            ySensorGood(1, mCntr, 4) = x(1, i + 1) + mNoiseSensorGood(7, mCntr);
            ySensorGood(2, mCntr, 4) = x(2, i + 1) + mNoiseSensorGood(8, mCntr);
            ySensorGood(1, mCntr, 5) = x(1, i + 1) + mNoiseSensorGood(9, mCntr);
            ySensorGood(2, mCntr, 5) = x(2, i + 1) + mNoiseSensorGood(10, mCntr);
        end
        
        %% Ten Good
        if ng == 10
            ySensorGood(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorGood(1, mCntr);
            ySensorGood(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorGood(2, mCntr);
            ySensorGood(1, mCntr, 2) = x(1, i + 1) + mNoiseSensorGood(3, mCntr);
            ySensorGood(2, mCntr, 2) = x(2, i + 1) + mNoiseSensorGood(4, mCntr);
            ySensorGood(1, mCntr, 3) = x(1, i + 1) + mNoiseSensorGood(5, mCntr);
            ySensorGood(2, mCntr, 3) = x(2, i + 1) + mNoiseSensorGood(6, mCntr);
            ySensorGood(1, mCntr, 4) = x(1, i + 1) + mNoiseSensorGood(7, mCntr);
            ySensorGood(2, mCntr, 4) = x(2, i + 1) + mNoiseSensorGood(8, mCntr);
            ySensorGood(1, mCntr, 5) = x(1, i + 1) + mNoiseSensorGood(9, mCntr);
            ySensorGood(2, mCntr, 5) = x(2, i + 1) + mNoiseSensorGood(10, mCntr);
            ySensorGood(1, mCntr, 6) = x(1, i + 1) + mNoiseSensorGood(11, mCntr);
            ySensorGood(2, mCntr, 6) = x(2, i + 1) + mNoiseSensorGood(12, mCntr);
            ySensorGood(1, mCntr, 7) = x(1, i + 1) + mNoiseSensorGood(13, mCntr);
            ySensorGood(2, mCntr, 7) = x(2, i + 1) + mNoiseSensorGood(14, mCntr);
            ySensorGood(1, mCntr, 8) = x(1, i + 1) + mNoiseSensorGood(15, mCntr);
            ySensorGood(2, mCntr, 8) = x(2, i + 1) + mNoiseSensorGood(16, mCntr);
            ySensorGood(1, mCntr, 9) = x(1, i + 1) + mNoiseSensorGood(17, mCntr);
            ySensorGood(2, mCntr, 9) = x(2, i + 1) + mNoiseSensorGood(18, mCntr);
            ySensorGood(1, mCntr, 10) = x(1, i + 1) + mNoiseSensorGood(19, mCntr);
            ySensorGood(2, mCntr, 10) = x(2, i + 1) + mNoiseSensorGood(20, mCntr);
        end
        
        if nc == 1 
            ySensorCheap(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorCheap(1, mCntr);
            ySensorCheap(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorCheap(2, mCntr);
        end

        %% Two cheap
        if nc == 2
            ySensorCheap(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorCheap(1, mCntr);
            ySensorCheap(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorCheap(2, mCntr);
            ySensorCheap(1, mCntr, 2) = x(1, i + 1) + mNoiseSensorCheap(3, mCntr);
            ySensorCheap(2, mCntr, 2) = x(2, i + 1) + mNoiseSensorCheap(4, mCntr);
        end
        %% Five cheap
        if nc == 5 
            ySensorCheap(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorCheap(1, mCntr);
            ySensorCheap(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorCheap(2, mCntr);
            ySensorCheap(1, mCntr, 2) = x(1, i + 1) + mNoiseSensorCheap(3, mCntr);
            ySensorCheap(2, mCntr, 2) = x(2, i + 1) + mNoiseSensorCheap(4, mCntr);
            ySensorCheap(1, mCntr, 3) = x(1, i + 1) + mNoiseSensorCheap(5, mCntr);
            ySensorCheap(2, mCntr, 3) = x(2, i + 1) + mNoiseSensorCheap(6, mCntr);
            ySensorCheap(1, mCntr, 4) = x(1, i + 1) + mNoiseSensorCheap(7, mCntr);
            ySensorCheap(2, mCntr, 4) = x(2, i + 1) + mNoiseSensorCheap(8, mCntr);
            ySensorCheap(1, mCntr, 5) = x(1, i + 1) + mNoiseSensorCheap(9, mCntr);
            ySensorCheap(2, mCntr, 5) = x(2, i + 1) + mNoiseSensorCheap(10, mCntr);
        end
        
        %% Ten Cheap
        if nc == 10
            ySensorCheap(1, mCntr, 1) = x(1, i + 1) + mNoiseSensorCheap(1, mCntr);
            ySensorCheap(2, mCntr, 1) = x(2, i + 1) + mNoiseSensorCheap(2, mCntr);
            ySensorCheap(1, mCntr, 2) = x(1, i + 1) + mNoiseSensorCheap(3, mCntr);
            ySensorCheap(2, mCntr, 2) = x(2, i + 1) + mNoiseSensorCheap(4, mCntr);
            ySensorCheap(1, mCntr, 3) = x(1, i + 1) + mNoiseSensorCheap(5, mCntr);
            ySensorCheap(2, mCntr, 3) = x(2, i + 1) + mNoiseSensorCheap(6, mCntr);
            ySensorCheap(1, mCntr, 4) = x(1, i + 1) + mNoiseSensorCheap(7, mCntr);
            ySensorCheap(2, mCntr, 4) = x(2, i + 1) + mNoiseSensorCheap(8, mCntr);
            ySensorCheap(1, mCntr, 5) = x(1, i + 1) + mNoiseSensorCheap(9, mCntr);
            ySensorCheap(2, mCntr, 5) = x(2, i + 1) + mNoiseSensorCheap(10, mCntr);
            ySensorCheap(1, mCntr, 6) = x(1, i + 1) + mNoiseSensorCheap(11, mCntr);
            ySensorCheap(2, mCntr, 6) = x(2, i + 1) + mNoiseSensorCheap(12, mCntr);
            ySensorCheap(1, mCntr, 7) = x(1, i + 1) + mNoiseSensorCheap(13, mCntr);
            ySensorCheap(2, mCntr, 7) = x(2, i + 1) + mNoiseSensorCheap(14, mCntr);
            ySensorCheap(1, mCntr, 8) = x(1, i + 1) + mNoiseSensorCheap(15, mCntr);
            ySensorCheap(2, mCntr, 8) = x(2, i + 1) + mNoiseSensorCheap(16, mCntr);
            ySensorCheap(1, mCntr, 9) = x(1, i + 1) + mNoiseSensorCheap(17, mCntr);
            ySensorCheap(2, mCntr, 9) = x(2, i + 1) + mNoiseSensorCheap(18, mCntr);
            ySensorCheap(1, mCntr, 10) = x(1, i + 1) + mNoiseSensorCheap(19, mCntr);
            ySensorCheap(2, mCntr, 10) = x(2, i + 1) + mNoiseSensorCheap(20, mCntr);
        end
        % ^^ This pattern for each sensor up to ??
        % 3D indices are ( {x/y}, {measurement id}, {sensor id} )
          This is where we can make an EVEN BETTER estimate via
        % various filtering methods, by combining together past and/or
        % current measure from either/both sensors and any number of
        % sensors, how exciting!

        %% COMBINING SENSING
        % Let's combine our sensors to make a variety of COMBINED estimates
        % CHEAP
        estSensorCheap(:, mCntr) = mean(ySensorCheap(:, mCntr, array(nc)), 3);  
        errEstSensorCheap(:, mCntr) = estSensorCheap(:, mCntr) - x(1:2, i + 1);
        
        % GOOD
        estSensorGood(:, mCntr) = mean(ySensorGood(:, mCntr, array(ng)), 3);  
        errEstSensorGood(:, mCntr) = estSensorGood(:, mCntr) - x(1:2, i + 1); %XXXX;

        % COMBO
        %% Calculates combo based on sum of sensors weight and std using Combo equation specifies in notes
        estCombo(:, mCntr) = (estSensorGood(:, mCntr) * (wg * ng) + estSensorCheap(:, mCntr) * (wc * nc)) / ( ng * wg + nc * wc );
        errEstCombo(:, mCntr) = estCombo(:, mCntr) - x(1:2, i + 1);
        
        %% TIME AVERAGING
        % Note that the real trick here is to not "look back" before the
        % "start of time", i.e., before we started measuring,
        % as indicated by an index in the second indice less than 1
        % So, at the start, we just average together all the measurements we
        % have until enough time has passed to have enough measurements

        % A moving - N POINT MOVING AVG
        j = 0;
        k = 0;
        estSensorGoodLast(:, mCntr) = [0; 0];
        estSensorCheapLast(:, mCntr) = [0; 0];
        NCAvg = lastXMeasuresC;  
        NGAvg = lastXMeasuresG;
        while j < min(NCAvg, mCntr)
            estSensorCheapLast(:, mCntr) = estSensorCheapLast(:, mCntr) + ySensorCheap(:, mCntr - j, 1);
            j = j + 1;
        end
        while k < min(NGAvg, mCntr)
            estSensorGoodLast(:, mCntr) = estSensorGoodLast(:, mCntr) + ySensorGood(:, mCntr - k, 1);
            k = k + 1;
        end

        estSensorCheapLast(:, mCntr) = estSensorCheapLast(:, mCntr) / j;
        errSensorCheapLast(:, mCntr) = estSensorCheapLast(:, mCntr) - x(1:2, i + 1);
        
        estSensorGoodLast(:, mCntr) = estSensorGoodLast(:, mCntr) / k;
        errSensorGoodLast(:, mCntr) = estSensorGoodLast(:, mCntr) - x(1:2, i + 1);

        mCntr = mCntr + 1;
    end

    i = i + 1; % increment our iteration counter
    time = time + dt; % increment time
    time_arr(i) = time;
end

%% find the mean and standard deviation of each of our measures
mu_SensorGood_PosX = mean(mNoiseSensorGood(1, :));
sigma_SensorGood_PosX = std(mNoiseSensorGood(1, :));
fprintf('Sensor, Good - PosX: \t\x03bc = %.3f, \x03C3 = %.3f\n', mu_SensorGood_PosX, sigma_SensorGood_PosX)

mu_SensorGood_PosY = mean(mNoiseSensorGood(2, :));
sigma_SensorGood_PosY = std(mNoiseSensorGood(2, :));
fprintf('Sensor, Good - PosY: \t\x03bc = %.3f, \x03C3 = %.3f\n', mu_SensorGood_PosY, sigma_SensorGood_PosY)

mu_SensorCheap_PosX = mean(mNoiseSensorCheap(1, :));
sigma_SensorCheap_PosX = std(mNoiseSensorCheap(1, :));
fprintf('Sensor, Cheap- PosX: \t\x03bc = %.3f, \x03C3 = %.3f\n', mu_SensorCheap_PosX, sigma_SensorCheap_PosX)

mu_SensorCheap_PosY = mean(mNoiseSensorCheap(2, :));
sigma_SensorCheap_PosY = std(mNoiseSensorCheap(2, :));
fprintf('Sensor, Cheap- PosY: \t\x03bc = %.3f, \x03C3 = %.3f\n', mu_SensorCheap_PosY, sigma_SensorCheap_PosY)

%% Find the mean and standard deviation in error in PositionX for each of the combined estimates
fprintf('\n--------------------\nFUSING MEASUREMENT\n')
% CHEAP
fprintf('* CHEAP -----\n')
mu_errEstSensorCheap = mean(errEstSensorCheap, 2);
sigma_errEstSensorCheap = std(errEstSensorCheap, 0, 2);
fprintf('Sensor, Cheap  %d-PosX: \t\x03bc = % 8.3f\t\x03C3 = %8.3f\n', nc, mu_errEstSensorCheap(1), sigma_errEstSensorCheap(1))

% GOOD
fprintf('* GOOD -----\n')
mu_errEstSensorGood = mean(errEstSensorGood, 2);
sigma_errEstSensorGood = std(errEstSensorGood, 0, 2);
fprintf('Sensor, Good   %d-PosX: \t\x03bc = % 8.3f\t\x03C3 = %8.3f\n', ng, mu_errEstSensorGood(1), sigma_errEstSensorGood(1))

% COMBO
fprintf('* COMBO -----\n')
mu_errEstCombo = mean(errEstCombo, 2);
sigma_errEstCombo = std(errEstCombo, 0, 2);
fprintf('Sensor, Combo: Using %d CHEAP & %d GOOD-PosX: \t\x03bc = % 8.3f\t\x03C3 = %8.3f\n', nc, ng, mu_errEstCombo(1), sigma_errEstCombo(1))

% AVERAGING
fprintf('* AVERAGING CHEAP -----\n')
mu_errSensorCheapLast = mean(errSensorCheapLast, 2);
sigma_errSensorCheapLast = std(errSensorCheapLast, 0, 2);
fprintf('Sensor, Cheap  %dL-PosX: \t\x03bc = % 8.3f\t\x03C3 = %8.3f\n', lastXMeasuresC, mu_errSensorCheapLast(1), sigma_errSensorCheapLast(1))

fprintf('* AVERAGING GOOD -----\n')
mu_errSensorGoodLast = mean(errSensorGoodLast, 2);
sigma_errSensorGoodLast = std(errSensorGoodLast, 0, 2);
fprintf('Sensor, Good  %dL-PosX: \t\x03bc = % 8.3f\t\x03C3 = %8.3f\n', lastXMeasuresG,  mu_errSensorGoodLast(1), sigma_errSensorGoodLast(1))

fprintf('--------------------\n\n')

%% Plot the RESULTS
if (dispPlots)
    subplot(5, 1, [1, 2])
    plot(x(1, :), x(2, :), ySensorGood(1, :, 1), ySensorGood(2, :, 1), ySensorCheap(1, :, 1), ySensorCheap(2, :, 1))
    legend('Real Position', 'Good Sensor', 'Cheap Sensor')
    title('y by x position, feet')
    sgtitle('Single Sensor')
    
    subplot(5, 1, 3)
    plot(time_arr, x(3, :))
    title('Timeline of V, knots')
    
    subplot(5, 1, 4)
    plot(time_arr, x(4, :))
    title('Timeline of Psi, degrees')
    
    subplot(5, 1, 5)
    plot(time_arr(2:end), totalxdot(4, :))
    title('Timeline of Psi_dot, degrees/sec')
    %%
    figure(100); clf;
    plot(x(1, :), x(2, :), 'DisplayName', 'Real Position')
    hold on
    
    for i = 1:size(ySensorGood, 3)
    plot(ySensorGood(1, :, i), ySensorGood(2, :, i), 'DisplayName', sprintf('Good Sensor %d', i))
    end
    
    for i = 1:size(ySensorCheap, 3)
    plot(ySensorCheap(1, :, i), ySensorCheap(2, :, i), 'DisplayName', sprintf('Cheap Sensor %d', i))
    end
    
    legend
    title('y by x position, feet')
    
    figure(101); clf;
    plot(x(1, :), x(2, :), 'DisplayName', 'Real Position')
    hold on
    plot(estSensorCheap(1, :), estSensorCheap(2, :), 'DisplayName', 'Cheap 1')
    legend
    title('y by x position, feet')
    
    figure(102); clf;
    plot(x(1, :), x(2, :), 'DisplayName', 'Real Position')
    hold on
    plot(estSensorGood(1, :), estSensorGood(2, :), 'DisplayName', 'Good 1')
    legend
    title('y by x position, feet')
    
    figure(103); clf;
    plot(x(1, :), x(2, :), 'DisplayName', 'Real Position')
    hold on
    plot(estCombo(1, :), estCombo(2, :), 'DisplayName', 'Combo')
    legend
    title('y by x position, feet')
end
function w = weight(isCheap, nC, nG, sigmaC, sigmaG )
    
    if (isCheap)
        w = (1/sigmaC^2) * ((1/(sigmaC^2)*nC) + (1/(sigmaG^2)*nG))^(-1);
    else 
        w = (1/sigmaG^2) * ((1/(sigmaC^2)*nC) + (1/(sigmaG^2)*nG))^(-1);
    end 
end

function arr = array(num)
    temp = [];
    for index=1:num
        temp(index) = index;
    end
    arr = temp;
end
