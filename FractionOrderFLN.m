%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Nonlinear System Identification using Fractional Order Functional Link Network (FRFT-FLN) Based Adaptive Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%       Let's Start!!!  %%%%%%%%%%%%%%%%%%%%%%
%% Default commands! 
clc;
clear all;

%% Defining the number of independent trials
% Number of trials to ensure ergodicity in the results
no_of_independent_trials = 10;

% Loop over independent trials
for itr = 1:no_of_independent_trials
    
    %% Displaying the number of the current independent trial
    clc;
    disp(['Independent Trial No: ', num2str(itr)])
    
    %% Initialization of Parameters
    % Number of input samples
    no_of_inputs = 6e6; 
    
    % Generate random input signal uniformly distributed in the range [-0.5, 0.5]
    input = rand(1, no_of_inputs) - 0.5;
    
    % Input buffer length
    N = 10; 
    
    % FRFT order (determines functional expansion complexity)
    FRFT_order = 2; 
    
    % Initialize input buffer with zeros
    x_buffer = zeros(1, N); 
    
    % Length of expanded input after functional expansion
    M = (2 * FRFT_order + 1) * N + 1; 
    
    % Initialize FRFT weights
    FRFT_weights = zeros(1, M); 
    
    % Step size for weight updates
    mu_weight = 0.01; 
    
    % Step size for fractional order updates
    mu_alpha = 0.01; 
    
    % Additive white Gaussian noise with a 50 dB noise floor
    noise = awgn(input, 50) - input; 
    
    % Initial fractional order (alpha) for FRFT
    alpha(1) = 0.5; 
    
    %% Begin FRFT-FLN Processing
    for i = 1:length(input)
        
        %% Update input buffer with current input sample
        x_buffer = [input(i) x_buffer(1:end-1)];
        
        %% Desired Output Generation
        % Nonlinear system output with sigmoid nonlinearity and added noise
        q = 1.5 * input(i) - 0.3 * input(i)^2; 
        if q > 0
            rho = 4; % Nonlinearity factor for positive values
        else
            rho = 0.5; % Nonlinearity factor for negative values
        end
        desired_output(i) = 2 * ((1 / (1 + exp(-rho * q))) - 0.5) + noise(i); 
        
        %% Functional Expansion Block (FEB) Generation
        FRFT_FEB = [];
        theta = alpha(i) * pi / 2; 
        a = 6;
        
        for k = 1:N
            for l = 1:FRFT_order
                fractional_mod_sin = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                fractional_mod_cos = cos(pi * l * x_buffer(k) * (1 + cos(theta)));
                
                FRFT_FEB = [FRFT_FEB, ...
                            exp(-a * alpha(i) * abs(x_buffer(k))) * fractional_mod_sin, ...
                            exp(-a * alpha(i) * abs(x_buffer(k))) * fractional_mod_cos];
            end
        end
        
        % Final FEB contents with bias and input buffer
        FRFT_FEB_final = [1, x_buffer, FRFT_FEB];
        
        %% Compute FRFT-FLN Output
        FRFT_output(i) = FRFT_weights * FRFT_FEB_final'; 
        
        % Error computation
        error(i) = desired_output(i) - FRFT_output(i);
        
        %% Update Rule for Fractional Order (Alpha)
        z = [];
        for k = 1:N
            for l = 1:FRFT_order
                fractional_mod_sin = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                fractional_mod_cos = cos(pi * l * x_buffer(k) * (1 + cos(theta)));
                
                d_fractional_mod_sin = (pi / 2) * fractional_mod_cos * (pi * l * x_buffer(k) * sin(theta));
                d_fractional_mod_cos = -(pi / 2) * fractional_mod_sin * (pi * l * x_buffer(k) * sin(theta));
                
                exp_term = exp(-a * alpha(i) * abs(x_buffer(k)));
                d_exp_term = -a * abs(x_buffer(k)) * exp_term;
                
                z_sin = d_exp_term * fractional_mod_sin + exp_term * d_fractional_mod_sin;
                z_cos = d_exp_term * fractional_mod_cos + exp_term * d_fractional_mod_cos;
                
                z = [z, z_sin, z_cos];
            end
        end
        z_final = [0, zeros(1, N), z];
        
        alpha(i + 1) = alpha(i) + mu_alpha * error(i) * z_final * FRFT_weights';
        
        %% Weight Update Rule
        FRFT_weights = FRFT_weights + mu_weight * error(i) * FRFT_FEB_final;
    end
    
    % Store squared error for the current trial
    err(itr, :) = error.^2;
end

%% Smoothing the Error Using Moving Average Filter
disp('Please Wait! Smoothing Operation is Going On...')
length_of_smoothing_filter = 200; % Length of smoothing filter
smoothing_filter_coeff = (1 / length_of_smoothing_filter) * ones(1, length_of_smoothing_filter); % Filter coefficients

% Apply smoothing filter to the error
for i = 1:no_of_independent_trials
    err_smooth(i, :) = filter(smoothing_filter_coeff, 1, err(i, :));
end

%% Plotting the Learning Curve
plot(10 * log10(mean(err_smooth)), 'b', 'DisplayName', 'FRFT-FLN');
xlabel('Iterations');
ylabel('MSE (dB)');
grid on;
legend show;

%% Compute and Display Final MSE
FRFT_mse = 10 * log10(mean(mean(err(:, end-1000:end))));
fprintf('Average MSE Value over the last 1000 iterations is %f\n', FRFT_mse);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% That's it! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
