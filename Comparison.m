%% Comparision of Non-linear System Identification using
% AEFLN, TFLN and FRFT-FLN Filters
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
    
    
    mu_theta = 0.01;
    mu_a = 0.1; 
    
    noise = awgn(input, 30) - input;

    a(1) = 0;
    alpha = 0.5*ones(1,length(input)); 
    amp = 6*ones(1,length(input));
    
    %-------------vector initialization-------------------------%
    g = zeros(1, no_of_inputs);
    g1 = zeros(1, no_of_inputs);
    g2 = zeros(1, no_of_inputs);
    g3 = zeros(1, no_of_inputs);
    g4 = zeros(1,no_of_inputs);

    for i = 3:no_of_inputs
        g1(i) = exp(0.5 * input(i)) * (sin(pi * input(i)) + 0.3 * sin(3 * pi * input(i - 2)) + 0.1 * sin(5 * pi * input(i)));
    end
    g1 = awgn(g1, 30);
    
    for i = 5:no_of_inputs
        g4(i) = 0.8*sin(pi*input(i))^3 - (2.5/(2.5 + input(i-1)^2)) - 0.15*cos(3*pi*input(i-4)) + ...
               0.1*sin(5*pi*input(i-2)) + 0.05*exp(-0.1*abs(input(i-3))) + 1.2;
%           g4(i) = sin(pi*input(i))*(input(i)^2);
    end
    g4 = awgn(g4, 30);


    for i = 1:no_of_inputs
        q(i) = (3/2) * input(i) - (3/10) * input(i)^2;
        rho = (q(i) > 0) * 4 + (q(i) <= 0) * 0.5;
        %3.5
        g2(i) = 2 * ((cos(q(i)) / (1 + exp(-3.5*rho * q(i)))) - 0.5);
    end
    g2 = awgn(g2, 30); 

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
    g3 = awgn(g3, 30);
    noise = awgn(input, 30) - input; 
    
    for i = 1:length(input)
        x_buffer = [input(i) x_buffer(1:end-1)];
        q = 1.5 * input(i) - 0.3 * input(i)^2; 
        if q > 0
            rho = 4; 
        else
            rho = 0.5; 
        end
        desired_output(i) = 2 * ((1 / (1 + exp(-rho * q))) - 0.5) + noise(i); 
        %--------------TFLN-----------------------------------------%
        TFLN_FEB = [1, x_buffer, sin(pi * x_buffer), cos(pi * x_buffer), sin(2 * pi * x_buffer), cos(2 * pi * x_buffer)];
        
        %---------------AEFLN---------------------------------------%
        AEFLN_FEB = [1, x_buffer, exp(-a(i) * abs(x_buffer)) .* sin(pi * x_buffer), exp(-a(i) * abs(x_buffer)) .* cos(pi * x_buffer), exp(-a(i)*abs(x_buffer)).*sin(2*pi*x_buffer),exp(-a(i)*abs(x_buffer)).*cos(2*pi*x_buffer)];
        
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
        
        %--------------------System changing-----------------%
%         if i < no_of_inputs / 3
%             g(i) = g2(i);
%         elseif i >= no_of_inputs / 3 && i < 2 * no_of_inputs / 3
%             g(i) = g3(i);
%         else
%             g(i) = g1(i);
%         end

        %Till now for G2 it's working great comperatively
        g(i) = g2(i);
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
Smooth_AEFLN = smooth(Error_AEFLN, N_smooth,'moving');
Smooth_TFLN = smooth(Error_TFLN, N_smooth,'moving');
Smooth_FRFTFLN = smooth(Error_FRFTFLN, N_smooth,'moving');

%----------Steady State MSE Values----------------%
MSE_AEFLN = 10 * log10(mean(Smooth_AEFLN(end - 1000:end)));
MSE_TFLN = 10 * log10(mean(Smooth_TFLN(end - 1000:end)));
MSE_FRFTFLN = 10 * log10(mean(Smooth_FRFTFLN(end - 1000:end)));

figure;
plot(10 * log10(Smooth_TFLN), 'b');
hold on;
plot(10 * log10(Smooth_AEFLN), 'r');
plot(10 * log10(Smooth_FRFTFLN), 'k');
xlabel('Iterations');
ylabel('MSE (dB)');
legend( 'TFLN', 'AEFLN',  'FRFT-FLN');
title('MSE Comparison of TFLN, AEFLN, and FRFT-FLN Filters');
hold off;


figure;
plot(alpha);
title('Alpha');
xlabel('alpha(i)');
ylabel('Iterations');

figure;
plot(amp);
title('Amp');
xlabel('amp(i)');
ylabel('Iterations');

figure;
plot(a);
title('Exponantial');
xlabel('a(i)');
ylabel('Iterations');

