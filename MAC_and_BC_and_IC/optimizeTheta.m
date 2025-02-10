function [theta_opt, bn_opt, Rxxs_opt, En_opt] = optimizeTheta(H, Lxu, theta_init, w, alpha, max_iters)
    theta = theta_init;
    U = length(theta);
    for iter = 1:max_iters
        % Compute current rates and energies
        [~, bn, Rxxs, En] = minPtoneMIMO(H, Lxu, theta, w);
        
        % Compute gradient of the objective with respect to theta
        grad_theta = computeGradientTheta(H, Lxu, theta, w, bn, En);
        
        % Update theta
        theta = theta + alpha * grad_theta;
        
        % Optional: Add convergence check based on the change in objective value or gradient magnitude
        if norm(grad_theta) < tolerance
            break;
        end
    end
    theta_opt = theta;
    [~, bn_opt, Rxxs_opt, En_opt] = minPtoneMIMO(H, Lxu, theta_opt, w);
end

function grad_theta = computeGradientTheta(H, Lxu, theta, w, bn, En)
    % Placeholder: this needs to be defined based on your objective and system model
    % This could involve numerical differentiation if analytical gradients are too complex
    grad_theta = zeros(size(theta)); % This needs actual implementation
end
