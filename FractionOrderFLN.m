clc; clear; close all;
rand('seed',0);

for itr = 1:10
    disp(['Independent Trial No: ', num2str(itr)])
    
    no_of_inputs = 18e6;
    input = rand(1, no_of_inputs) - 0.5;
    N = 10;
    FRFT_order = 2;
    x_buffer = zeros(1, N);
    M = (2 * FRFT_order + 1) * N + 1;

    FRFT_weights = zeros(1, M);
    %------learning rate---------------%
    mu_weight = 0.01;
    mu_theta = 0.01;
    mu_a = 0.1;
    
    noise = awgn(input, 30) - input;
    
    alpha = 0.5*ones(1,length(input)); 
    amp = 6*ones(1,length(input));
    
    for i = 1:length(input)
        x_buffer = [input(i) x_buffer(1:end-1)];
        q = 1.5 * input(i) - 0.3 * input(i)^2; 
        if q > 0
            rho = 4;
        else
            rho = 0.5;
        end
        desired_output(i) = 2 * ((1 / (1 + exp(-rho * q))) - 0.5) + noise(i); 
        
        %------------FRFT-FLN--------------------------------------%
        FRFT_FEB = [];
        theta = alpha(i); 
        
        for k = 1:N
            for l = 1:FRFT_order
                fractional_mod_sin = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                fractional_mod_cos = cos(pi * l * x_buffer(k) * (1 + cos(theta)));
                
                FRFT_FEB = [FRFT_FEB, ...
                            exp(-amp(i) * abs(x_buffer(k))) * fractional_mod_sin, ...
                            exp(-amp(i) * abs(x_buffer(k))) * fractional_mod_cos];
            end
        end
        
        %-----------------output & error calculation--------------------%
        FRFT_FEB_final = [1, x_buffer, FRFT_FEB];
        FRFT_output(i) = FRFT_weights * FRFT_FEB_final';
        error_FTFTFLN(i) = desired_output(i) - FRFT_output(i);

        z = [];
        for k = 1:N
            for l = 1:FRFT_order
                fractional_mod_sin = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                fractional_mod_cos = cos(pi * l * x_buffer(k) * (1 + cos(theta)));
                
                e_exp_term = - abs(x_buffer(k)) * exp(-amp(i) * abs(x_buffer(k)));
                
                z_sin = e_exp_term * fractional_mod_sin;
                z_cos = e_exp_term * fractional_mod_cos;
                
                z = [z, z_sin, z_cos];
            end
        end
        z_final = [0, zeros(1, N), z];
        
    %% Update Rule for Fractional Order (alpha/theta)
        v = [];
        for k = 1:N
            for l = 1:FRFT_order
                v_fractional_mod_sin = cos(pi * l * x_buffer(k) * (1 + cos(theta)));
                v_fractional_mod_cos = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                
                v_d_fractional_mod_sin =  v_fractional_mod_cos * (pi * l * x_buffer(k) * sin(theta));
                v_d_fractional_mod_cos = - v_fractional_mod_sin * (pi * l * x_buffer(k) * cos(theta));
                
                v_exp_term = exp(-amp(i) * abs(x_buffer(k)));
                v_d_exp_term = - abs(x_buffer(k)) * v_exp_term;
                
                v_sin = v_exp_term * v_d_fractional_mod_sin ;
                v_cos = v_exp_term * v_d_fractional_mod_cos ;
                
                v = [v, v_sin, v_cos];
            end
        end
        v_final = [0, zeros(1, N), v];
         

        alpha(i+1) = alpha(i) + mu_theta * error_FTFTFLN(i) * v_final * FRFT_weights';     
        amp(i + 1) = amp(i) + mu_a * error_FTFTFLN(i) * z_final * FRFT_weights';
        
        %-----------------Weight update--------------------%
        FRFT_weights = FRFT_weights + mu_weight * error_FTFTFLN(i) * FRFT_FEB_final;
    end
    Err_FRFTFLN(itr, :) = error_FTFTFLN .^ 2;
end

Error_FRFTFLN = mean(Err_FRFTFLN);

N_smooth = 500;
%==============FsLMS===========================%
Smooth_FRFTFLN = smooth(Error_FRFTFLN, N_smooth,'moving');

%----------Steady State MSE Values----------------%
MSE_FRFTFLN = 10 * log10(mean(Smooth_FRFTFLN(end - 1000:end)));

figure;
plot(10 * log10(Smooth_FRFTFLN), 'k');
xlabel('Iterations');
ylabel('MSE (dB)');
legend('FRFT-FLN');
title('MSE of FRFT-FLN Filter');
