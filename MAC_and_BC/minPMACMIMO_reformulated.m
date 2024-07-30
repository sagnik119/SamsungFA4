function [FEAS_FLAG, bu_a, info] = minPMACMIMO_reformulated(H, Lxu, bu_min, w, cb)
% TODO: do for cb = 2

    % start timer
    tic;

    bu_min = cb*bu_min;

    % segregate into channel for each user
    Lxu_cum = cumsum(Lxu);
    U = size(Lxu, 2);
    [Ly, ~, N] = size(H);
    H_cell = cell(U, 1);
    H_cell{1} = H(:,1:Lxu(1),:);
    for i=2:U
        H_cell{i} = H(:,Lxu_cum(i-1)+1:Lxu_cum(i), :);
    end
    for i=1:U
        H_cell{i} = reshape(H_cell{i}, size(H_cell{i}, 1), size(H_cell{i}, 2), 1, size(H_cell{i}, 3));
    end
    
    
    subsets = powerSet(1:U);
    subsets = subsets(2:end);

    cvx_begin quiet
        cvx_solver mosek
        
        max_Lxu = max(Lxu);
        variable Rxx(max_Lxu, max_Lxu, U, N) hermitian semidefinite

        variable b(U, N) nonnegative
        dual variable theta
        expressions rate_sum obj_indiv

        % set unneeded elements to 0
        for u = 1:U
            for n = 1:N
                if size(Rxx(:,:,u,n), 1) > Lxu(u)
                    Rxx(Lxu(u)+1:size(Rxx(:,:,u,n)), Lxu(u)+1:size(Rxx(:,:,u,n)), u, n) == 0;
                    Rxx(1:Lxu(u), Lxu(u)+1:size(Rxx(:,:,u,n)), u, n) == 0;
                    Rxx(Lxu(u)+1:size(Rxx(:,:,u,n)), 1:Lxu(u), u, n) == 0;
                    % disp("making 0");
                end
            end
        end

        % for u = 1:U
        %     for n = 1:N
        %         Rxx(1:Lxu(u), 1:Lxu(u), u, n) == semidefinite(Lxu(u));
        %     end
        % end

        theta : sum(b, 2)' >= bu_min;
        
    %     for u = 1:U
    %         for n = 1:N
    %             b(u, n) >= 0;
    %         end
    %     end
        
    %     for u = 1:U
    %         for n = 1:N
    %             for i = 1:Lxu(1)
    %                 Rxx(i, i, u, n) >= 0;
    %             end
    %         end
    %     end
        
         for subset = subsets
             for n = 1:N
                 temp = zeros(Ly, Ly);
                 for u = subset{:}
                     temp = temp + H_cell{u}(:, :, 1, n) * Rxx(1:Lxu(u), 1:Lxu(u), u, n) * H_cell{u}(:, :, 1, n)';
                 end
%                  for u = subset{:}
%                     temp = temp + H(:,index_start(u):index_end(u), n) * Rxx(:, :, u, n) * H(:,index_start(u):index_end(u), n)';
%                  end
    
    
                 capacity_bound = log_det(eye(Ly) + temp)/log(2);
                 rate_sum = 0;
                 for u = subset{:}
                    rate_sum = rate_sum + b(u, n);
                 end
                 rate_sum <= capacity_bound;
             end
         end
         obj = 0;
         for u = 1:U
             for n =1:N
                obj = obj + w(u)*trace(Rxx(1:Lxu(u),1:Lxu(u),u, n));
             end
         end
         
         % adding that no energy can cross 1
         % TODO: replace with maximum energy limit
         for u = 1:U
             obj_indiv = 0;
             for n =1:N
                obj_indiv = obj_indiv + w(u)*trace(Rxx(1:Lxu(u),1:Lxu(u),u, n));
             end
             obj_indiv = obj_indiv / N;
             obj_indiv <= 1;
         end

         minimize (obj)
    cvx_end
    % disp(theta)
    % disp("Rxx is");
    % disp(Rxx);
    % disp("The b from initial convex programming is");
    % disp(sum(b, 2));

    bu_min = 1/cb * bu_min;
    b = 1/cb * b;
    
    % calculate energies
    Eu_a = zeros(1, U);
    for u=1:U
        for n=1:N
            Eu_a(u) = Eu_a(u) + trace(reshape(Rxx(1:Lxu(u), 1:Lxu(u), u, n), Lxu(u), Lxu(u)));
        end
    end
    Eu_a = Eu_a/N;
    % theta = theta';

    % temporary - remove later
    % disp("Rates from initial cvx");
    % disp(sum(b, 2)');

    % Cluster users by theta values
    [clusters, unique_theta, theta] = identify_clusters(theta);
    unique_theta = unique_theta';

    % disp("clusters is");
    % disp(clusters);

    theta = theta';
    
    % Generate all possible orders considering clusters
    all_orders = generate_all_orders(clusters);
    % disp("all orders is ");
    % disp(all_orders);

    flag_time_sharing_required = 0;
    if (size(all_orders, 1) > 1)
        % only one order exists, no time sharing
        flag_time_sharing_required = 1;
    else
        FEAS_FLAG = 1;
    end

    if (flag_time_sharing_required==1)
        disp("Starting time sharing");
        % Find optimal weights for these orders
        [weights, orderings, Rxx_orders, Eun_orders, bun_orders, bu_a_orders, num_orders] = find_optimal_weights(H_cell, Rxx, all_orders, unique_theta, bu_min, w, Ly, U, N, Lxu, cb);
        disp("Time sharing finished");
        FEAS_FLAG = sum(weights > 0) > 0;  % Feasibility flag 0-infeasible, 1-feasible
    end
    
    if (FEAS_FLAG == 1 && flag_time_sharing_required == 1)
        if (sum(weights > 0)) > 1  % Feasibility flag 1-no time sharing, 2-time sharing
            FEAS_FLAG = 2;
        end
    end

    if (flag_time_sharing_required == 1)
        % calculate bu_a from bu_a_orders
        weights = weights';
        weights_reshaped = reshape(weights, [1, num_orders]);
        weights_reshaped = repmat(weights_reshaped, [U, 1]);
        bu_a = sum(bu_a_orders .* weights_reshaped, 2);
        bu_a = bu_a';
    else
        [b_achieved, bu_a, Eun] = compute_achieved_rates_order_given(H, Rxx, all_orders, Ly, U, N, cb);
    end
    
    
    
    % Store information
    info.H = H;
    info.H_cell = H_cell;
    info.bu_min = bu_min;
    info.w = w;
    info.cb = cb;
    info.unique_theta = unique_theta;
    info.theta = theta;
    info.clusters = clusters;
    if (flag_time_sharing_required == 1)
        info.orderings = orderings;
        info.time_sharing_weights = weights;
        info.Rxx_orders = Rxx_orders;
        info.Eun_orders = Eun_orders;
        info.bun_orders = bun_orders;
        info.bu_a_orders = bu_a_orders;
        info.bu_a = bu_a;
    else
        info.Rxx = Rxx;
        info.Eun = Eun;
        info.bun = b_achieved;
        info.bu_a = bu_a;
    end
    info.FEAS_FLAG = FEAS_FLAG;
    info.Rxx = Rxx;
    info.Eu_a = Eu_a;

    elapsedTime = toc;
    disp(['Elapsed time for reformulated minPMAC: ', num2str(elapsedTime), ' seconds']);

end
% 
function [weights, orderings, Rxx_orders, Eun_orders, bun_orders, bu_a_orders, num_orders] = find_optimal_weights(H_cell, Rxx, all_orders, theta, bu_min, w, Ly, U, N, Lxu, cb)
    num_orders = size(all_orders, 1);
    % disp(num_orders);
    orderings = cell(num_orders, 1);
    Rxx_orders = cell(num_orders, 1);
    Eun_orders = zeros(U, N, num_orders);
    bu_a_orders = zeros(U, num_orders);
    b_achieved_orders = zeros(U, N, num_orders);

    % Compute the achieved rates for each ordering
    for i = 1:num_orders
        orderings{i} = all_orders(i, :);
        [b_achieved_orders(:, :, i), bu_a_orders(:,i), Eun_orders(:,:,i)] = compute_achieved_rates_order_given(H_cell, Rxx, orderings{i}, Ly, U, N, Lxu, cb);
    end
    
    b_achieved_orders = real(b_achieved_orders);

    % disp("b_achieved_orders is");
    % disp(size(b_achieved_orders));
    % disp(b_achieved_orders);

    % Solve for the optimal weights using mixed-integer linear programming
    cvx_begin quiet
        cvx_solver mosek
        variable weights(num_orders) nonnegative
        variable z(num_orders) binary
        minimize(sum(z))  % Minimize the number of non-zero weights
        subject to
            weight_sum = sum(weights);
            weight_sum == 1;
            % disp("bu_min is");
            % disp((bu_min));
            % disp(size(sum(sum(b_achieved_orders .* reshape(weights, [U  num_orders]), 3), 2)'));
            sum(sum(b_achieved_orders .* repmat(reshape(weights, [1 1 num_orders]), [U, N, 1]), 3), 2)' >= bu_min;
            weights <= z;  % Ensure weights are zero where z is zero
    cvx_end

    % if weights is below tolerance, make it 0
    weights(weights <= 1e-7) = 0;

    % disp("weights is");
    % disp(weights);

    % Store the results
    for i = 1:num_orders
        if weights(i) > 0
            Rxx_orders{i} = Rxx;
        end
    end
    % correct the above - obviously wrong
    
    bun_orders = b_achieved_orders;

end

function [b_achieved, bu_a, Eun] = compute_achieved_rates_order_given(H_cell, Rxx, ...
        order, Ly, U, N, Lxu, cb)
    % order given
    pi_inv = order;
    % disp("pi_inv is");
    % disp(pi_inv);

    % Initialize b_achieved
    b_achieved = zeros(U, N);
    
    % Calculate b_achieved
    for n = 1:N
        for u = 1:U
            mat1 = 0;
            for i = u:U
                % disp(length(H_cell));
                % disp(pi_inv(i));
                % disp(size(H_cell{pi_inv(i)}));
                % disp(size(Rxx(1:Lxu(pi_inv(i)), 1:Lxu(pi_inv(i)))));
                mat1 = mat1 + H_cell{pi_inv(i)}(:, :, 1, n) * Rxx(1:Lxu(pi_inv(i)), 1:Lxu(pi_inv(i)), pi_inv(i), n) * H_cell{pi_inv(i)}(:, :, 1, n)';
            end
            mat2 = 0;
            for i = u+1:U
                mat2 = mat2 + H_cell{pi_inv(i)}(:, :, 1, n) * Rxx(1:Lxu(pi_inv(i)), 1:Lxu(pi_inv(i)), pi_inv(i), n) * H_cell{pi_inv(i)}(:, :, 1, n)';
            end
            % disp("mats are");
            % disp((eye(Ly) + mat1));
            % disp((eye(Ly) + mat2))
            % disp("dets are");
            % disp(det(eye(Ly) + mat1));
            % disp(det(eye(Ly) + mat2));
            % if ((det(eye(Ly) + mat1)) == (det(eye(Ly) + mat2)))
            %     disp("Eureka");
            %     disp(u);
            % end
            b_achieved(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
         end
    end

    % convert to real number
    b_achieved = real(b_achieved);

    b_achieved = b_achieved/cb;
    
    % Compute bu_a
    bu_a = sum(b_achieved, 2);
    
    % Compute Eun
    Eun = zeros(U, N);
    for u = 1:U
        for n = 1:N
            Eun(u, n) = trace(Rxx(1:Lxu(u), 1:Lxu(u), u, n));
        end
    end

end

function pSet = powerSet(S)
    n = length(S);            
    numSubsets = 2^n;         
    pSet = cell(1, numSubsets); 

    for i = 0:numSubsets-1
        binaryIndex = bitget(i, n:-1:1);  
        subset = S(logical(binaryIndex)); 
        pSet{i+1} = subset;               
    end
end

function all_orders = generate_all_orders(clusters)
    % Generate all permutations of users considering clusters
    all_orders = perms(clusters{1});
    for i = 2:length(clusters)
        new_orders = perms(clusters{i});
        all_orders = combvec(all_orders, new_orders);
        % disp("inside loop all_orders");
        % disp(all_orders);
    end
    if length(clusters) > 1
        all_orders = all_orders';
    end
end

function [clusters, unique_theta, theta] = identify_clusters(theta)
        % Find unique thetas
        % Define the tolerance
        tolerance = 1e-5;
        
        % Sort the theta values
        sorted_theta = sort(theta);
        
        % Initialize the unique_theta array with the first value
        unique_theta = sorted_theta(1);
        
        % Loop through the sorted theta values and check the difference
        for i = 2:length(sorted_theta)
            if abs(sorted_theta(i) - unique_theta(end)) > tolerance
                unique_theta = [unique_theta; sorted_theta(i)];
            end
        end
        % disp("Size of unique theta is");
        % disp(size(unique_theta));
        
        % Initialize clusters as a cell array
        clusters = cell(length(unique_theta), 1);
        
        % Assign users to clusters based on their theta values
        for i = 1:length(unique_theta)
            clusters{i} = find(abs(theta - unique_theta(i)) < tolerance);
        end
        
        % Display the clusters and their corresponding theta values
        % disp('Clusters and their corresponding theta values:');
        % for i = 1:length(clusters)
        %     fprintf('Theta value: %f', unique_theta(i));
        %     % disp(clusters{i});
        % end
    end

function [b_achieved, bu_a, Eun] = compute_achieved_rates(H, Rxx, ...
        theta, Ly, U, N)
    % Sort theta and compute pi and pi_inv
    [~, pi_inv] = sort(theta);
    [~, pi] = sort(pi_inv);
    
    % Initialize b_achieved
    b_achieved = zeros(U, N);
    
    % Calculate b_achieved
    for n = 1:N
        for u = 1:U
            mat1 = 0;
            for i = u:U
                mat1 = mat1 + H(:, pi_inv(i), n) * Rxx(:, :, pi_inv(i), n) * H(:, pi_inv(i), n)';
            end
            mat2 = 0;
            for i = u+1:U
                mat2 = mat2 + H(:, pi_inv(i), n) * Rxx(:, :, pi_inv(i), n) * H(:, pi_inv(i), n)';
            end
            
            b_achieved(pi_inv(u), n) = log2((det(eye(Ly) + mat1))) - log2((det(eye(Ly) + mat2)));
        end
    end
    
    % Transform theta
    theta = theta';
    
    % Compute bu_a
    bu_a = sum(b_achieved, 2);
    
    % Compute Eun
    Eun = zeros(U, N);
    for u = 1:U
        for n = 1:N
            Eun(u, n) = trace(Rxx(:, :, u, n));
        end
    end
end