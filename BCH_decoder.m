function decoded_msg= BCH_decoder(alpha,codeword,k,t,m,px)

% alpha_power = alpha.^(0:n-1); %see from primitive poly that a^7 = a^3+1 in decimal
% alpha_power = double(alpha_power.x); %.x so conversion from gf to double/ in the form a^(alpha_power),alpha^128=alpha
% 
% alpha_power_index = zeros(1, n); % find GF element in alpha_power array 
% for i = 2: n
%   alpha_power_index(i) = find(i == alpha_power(1:n))-1; %the index value for alpha_power
% end
% alpha_power_index(1) = n; %setting the first value to 1 because alpha^127 = 1 

[alpha_power_index] = gen_gf_table(px);
%the relationship between alpha_power and alpha_power_index is very
%important for decoding purpose
% y is the decimal representation, x is the power of alpha
% alpha_power(x)=y >>>>>> alpha^(x-1) = y(in decimal) %input power x, output decimal for x-1
% alpha_power_index(y)=x-1 >>>>>> alpha^(x-1) = y(in decimal) % input decimal y, output (power+1) for y-1
% example 
% alpha_power(24)=47 >>>>> alpha^23 = 47
% alpha_power_index(47) = 23 >>>> alpha^23 = 47

% calculate syndrome
syndrome = polyval(gf(codeword,m),alpha.^(1:2*t));

% error detection and correction

if isequal(syndrome,gf(zeros(1,2*t),m))
    decoded_msg = codeword(1:k);    % no error detected

else % error occured
    % to apply Berlekamp iterative method

    % Initialisation of iteration -1
    sigma_rho = [gf(zeros(1,t),m),1]; %coefficients of X in descending order, last element is 1
    sigma_mu = [gf(zeros(1,t),m),1]; 
    
    
    d_rho = gf(1,m); % assume 1 at start
    d_mu = syndrome(1); % if zero then sigma_miu+1 = sigma_miu
    l_mu = 0; %highest order of sigma_miu
    rho = -1; %rho doesn't exist for iteration -1, set to -1 for simplicity
    %rho changes every two iteration and and it equals mu-l_mu when d is not 0
    
    % iteration 0 still setting up
    if d_mu ~= 0 % discrepany between current iteration and previous one
        sigma_mu = sigma_rho+[zeros(1,t-1) d_mu 0]; %update sigma from previous sigma
        d_rho = d_mu; % update d_rho to d_mu
        l_mu = l_mu+1; %since discrepany happens sigma order increases
        rho = 0; % at d_mu ~= 0 
    end

    for mu=1:2*t-1 % start actual iteration
        d_mu = syndrome(mu+1); %update d_mu to be the next syndrome
        for i=1:l_mu
            d_mu = d_mu+sigma_mu(end-i)*syndrome(mu+1-i); % find d_mu in position of syndrome
        end
        if d_mu ~= 0 % still error
            new_sigma = conv([d_mu*(1/d_rho) zeros(1,mu-rho)],sigma_rho); %update sigma from previous values
            temp = sigma_mu;
            sigma_mu = sigma_mu + new_sigma(end-t:end); 
            
            if 2*l_mu<=mu+1 %check for the highest value and update rho,l_mu and
                sigma_rho = temp;
                d_rho = d_mu;
                l_mu = mu+1-l_mu;
                if l_mu > t
                    decoded_msg = codeword;
                    return
                end
                rho = mu; % update rho when d_mu not equal 0 and mu-l_mu is the highest
            end
        end
    end

    % determine error-location numbers by finding the roots of sigma
    roots_error_locations = roots(sigma_mu(end,:)); %the length should be equal to highest degree of sigma
    
    err_loc_index = alpha_power_index(double(roots_error_locations.x)); %find the power of alpha thus location
    
    error_pattern = gf(zeros(size(codeword)),m);
    error_pattern(err_loc_index) = gf(1,m); % flag it 1 meaning error occured at this position
    
    % error correction
    decoded_msg = gf(codeword,m) + error_pattern; 

    if any(polyval(decoded_msg,alpha)) %still error
        decoded_msg = codeword; % exceed error correction ability 
    else
        decoded_msg = double(decoded_msg.x); %successfully decoded
    end
    
    
    
end


end