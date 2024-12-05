% Comparision of Non-linear System Identification using
% AEFLN, TFLN and FRFT-FLN Filters

clear all;
clc;
rand('seed',0);

for itr = 1:2
    disp(['Independent Trial No: ', num2str(itr)])
    
    no_of_inputs = 18e6;
    input = rand(1, no_of_inputs) - 0.5;
    N = 10;
    FRFT_order = 2;
    x_buffer = zeros(1, N);
    M = (2 * FRFT_order + 1) * N + 1;
    
    P1 = 2;
    fln_expansion = (2 * P1 + 1) * N + 1;
    TFLN_weights = zeros(1, fln_expansion);

    P2 = 2;
    AEFLN_order = 2;
    aefln_expansion = (2 * P2 + 1) * N + 1;
    AEFLN_weights = zeros(1, aefln_expansion);


    FRFT_weights = zeros(1, M);
    %------learning rate---------------%
    mu_weight = 0.01;
    mu_aefln = 0.01;
    mu_fln = 0.01;
    mu_alpha = 0.01;
    noise = awgn(input, 30) - input;

    alpha(1) = 0.1;
    a(1) = 0.5;
    a_integral = 0;
    
    %-------------vector initialization-------------------------%
    g = zeros(1, no_of_inputs);
    g1 = zeros(1, no_of_inputs);
    g2 = zeros(1, no_of_inputs);
    g3 = zeros(1, no_of_inputs);

    for i = 3:no_of_inputs
        g1(i) = exp(0.5 * input(i)) * (sin(pi * input(i)) + 0.3 * sin(3 * pi * input(i - 2)) + 0.1 * sin(5 * pi * input(i)));
    end
    g1 = awgn(g1, 50); 

    for i = 1:no_of_inputs
        q(i) = (3/2) * input(i) - (3/10) * input(i)^2;
        rho = (q(i) > 0) * 4 + (q(i) <= 0) * 0.5;
        g2(i) = 2 * ((1 / (1 + exp(-rho * q(i)))) - 0.5);
    end
    g2 = awgn(g2, 50); 

    chi = 0.1;
    for i = 1:no_of_inputs
        if abs(input(i)) >= 0 && abs(input(i)) < chi
            g3(i) = (2 / (3 * chi)) * input(i);
        elseif abs(input(i)) >= chi && abs(input(i)) < (2 * chi)
            g3(i) = sign(input(i)) * (3 - (2 - abs(input(i)) / chi)^2) / 3;
        else
            g3(i) = sign(input(i));
        end
    end
    g3 = awgn(g3, 50);

    for i = 1:length(input)
        x_buffer = [input(i) x_buffer(1:end-1)];
         
        %--------------TFLN-----------------------------------------%
        TFLN_FEB = [1, x_buffer, sin(pi * x_buffer), cos(pi * x_buffer), sin(2 * pi * x_buffer), cos(2 * pi * x_buffer)];
        
        %---------------AEFLN---------------------------------------%
        AEFLN_FEB = [1, x_buffer, exp(-a(i) * abs(x_buffer)) .* sin(pi * x_buffer), exp(-a(i) * abs(x_buffer)) .* cos(pi * x_buffer), exp(-a(i)*abs(x_buffer)).*sin(2*pi*x_buffer),exp(-a(i)*abs(x_buffer)).*cos(2*pi*x_buffer)];
        
        %------------FRFT-FLN--------------------------------------%
        FRFT_FEB = [];
        theta = alpha(i) * pi / 2;
        amp = 6;
        
        for k = 1:N
            for l = 1:FRFT_order
                
                fractional_mod_sin = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                fractional_mod_cos = cos(pi * l * x_buffer(k) * (1 + cos(theta)));
                
                FRFT_FEB = [FRFT_FEB, exp(-amp*alpha(i) * abs(x_buffer(k))) * fractional_mod_sin, ...
                            exp(-amp*alpha(i) * abs(x_buffer(k))) * fractional_mod_cos];
            end
        end
        
        %--------------------System changing-----------------%
        if i < no_of_inputs / 3
            g(i) = g2(i);
        elseif i >= no_of_inputs / 3 && i < 2 * no_of_inputs / 3
            g(i) = g3(i);
        else
            g(i) = g1(i);
        end

        %-----------------output & error calculation--------------------%
        FRFT_FEB_final = [1, x_buffer, FRFT_FEB];
        FRFT_output(i) = FRFT_weights * FRFT_FEB_final';
        error_FTFTFLN(i) = g(i) - FRFT_output(i);
        
        TFLN_output(i) = TFLN_weights * TFLN_FEB';
        AEFLN_output(i) = AEFLN_weights * AEFLN_FEB';
        
        error_TFLN(i) = g(i) - TFLN_output(i);
        error_AEFLN(i) = g(i) - AEFLN_output(i);


        z = [];
        for k = 1:N
            for l = 1:FRFT_order
                fractional_mod_sin = sin(pi * l * x_buffer(k) * (1 + cos(theta)));
                fractional_mod_cos = cos(pi * l * x_buffer(k) * (1 + cos(theta)));

                d_fractional_mod_sin = (pi/2) * fractional_mod_cos * (pi*l*x_buffer(k) * sin(theta));
                d_fractional_mod_cos = -(pi/2) * fractional_mod_sin * (pi*l*x_buffer(k) * sin(theta));


                exp_term = exp(-amp*alpha(i) * abs(x_buffer(k)));
                d_exp_term = -amp*abs(x_buffer(k)) * exp_term;

                z_sin = d_exp_term * fractional_mod_sin ...
                        + exp_term * d_fractional_mod_sin;

                z_cos = d_exp_term * fractional_mod_cos ...
                        + exp_term * d_fractional_mod_cos;

                z = [z, z_sin, z_cos];
            end
        end
        z_final = [0, zeros(1, N), z];
         

        alpha(i+1) = alpha(i) + mu_alpha* error_FTFTFLN(i) * z_final * FRFT_weights';
        temp_vect=[zeros(1,N+1) ((-abs(x_buffer).*exp(-a(i)*abs(x_buffer))).*sin(pi*x_buffer))  ((-abs(x_buffer).*exp(-a(i)*abs(x_buffer))).*cos(pi*x_buffer)) ((-abs(x_buffer).*exp(-a(i)*abs(x_buffer))).*sin(2*pi*x_buffer))  ((-abs(x_buffer).*exp(-a(i)*abs(x_buffer))).*cos(2*pi*x_buffer))];
        a(i+1)=a(i)+0.1*error_AEFLN(i)*AEFLN_weights*temp_vect';
        
        %-----------------Weight update--------------------%
        TFLN_weights = TFLN_weights + mu_fln * error_TFLN(i) * TFLN_FEB;
        AEFLN_weights = AEFLN_weights + mu_aefln * error_AEFLN(i) * AEFLN_FEB;
        FRFT_weights = FRFT_weights + mu_weight * error_FTFTFLN(i) * FRFT_FEB_final;
    end

    Err_FLN(itr, :) = error_TFLN .^ 2;
    Err_AEFLN(itr, :) = error_AEFLN .^ 2;
    Err_FRFTFLN(itr, :) = error_FTFTFLN .^ 2;
end

Error_AEFLN = mean(Err_AEFLN);
Error_TFLN = mean(Err_FLN);
Error_FRFTFLN = mean(Err_FRFTFLN);

N_smooth = 1000;
%==============FsLMS===========================%
Smooth_AEFLN = moving_average(Error_AEFLN, N_smooth);
Smooth_TFLN = moving_average(Error_TFLN, N_smooth);
Smooth_FRFTFLN = moving_average(Error_FRFTFLN, N_smooth);

%----------Steady State MSE Values----------------%
MSE_AEFLN = 10 * log10(mean(Smooth_AEFLN(end - 1000:end)));
MSE_TFLN = 10 * log10(mean(Smooth_TFLN(end - 1000:end)));
MSE_FRFTFLN = 10 * log10(mean(Smooth_FRFTFLN(end - 1000:end)));

figure;
plot(10 * log10(Smooth_FRFTFLN), 'k');
hold on;
plot(10 * log10(Smooth_TFLN), 'b');
plot(10 * log10(Smooth_AEFLN), 'r');
xlabel('Iterations');
ylabel('MSE (dB)');
legend('FRFT-FLN', 'TFLN', 'AEFLN');
title('MSE Comparison of AEFLN, TFLN, and FRFT-FLN Filters');
hold off;


function smoothed_data = moving_average(data, window_size)
    smoothed_data = filter(ones(1, window_size) / window_size, 1, data);
end