%% construct L by N matrix holding unary cost 
function [cost, U, labels] = getUnaryCost(nodes, params)
    %%%%%%%%%%%%%%%%%%%%
    % Cost definition
    %%%%%%%%%%%%%%%%%%%%
    sig = params.sig;
    mu_b = params.mu_b; mu_f1 = params.mu_f1; mu_f2 = params.mu_f2;    
    const = 1/2*log(2*pi) + log(sig);
    alpha_b = (nodes - mu_b).^2./(2*sig^2) + const;
    aa = exp(-(nodes - mu_f1).^2./(2*sig^2)) +  exp(-(nodes - mu_f2).^2./(2*sig^2));
    alpha_f = -log(aa + eps) + const + log(2);
    cost = [alpha_b, alpha_f]';
    [maxLikelihood argmax] = min(cost);
    labels = argmax - 1;
    U = sum(maxLikelihood);
end
