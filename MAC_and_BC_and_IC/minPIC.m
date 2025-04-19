function [best_sol, best_powers, best_data_rates] = minPIC(U, H, min_rate)
    % Generate all permutations of users for decoding order
    permutations = perms(1:U);
    
    % -- 1) Bracketing
    lowFactor = 1.0;
    highFactor = lowFactor;
    maxFactor = 1e6;
    found_tight = false;
    
    best_sol = Inf;
    best_powers = [];
    best_data_rates = [];

    while highFactor <= maxFactor
        [test_sol, test_pw, test_dr, c_tight] = solve_for_imp_factor(U, H, min_rate, permutations, highFactor);
        if c_tight
            found_tight = true;
            best_sol = test_sol;
            best_powers = test_pw;
            best_data_rates = test_dr;
            break
        else
            lowFactor = highFactor;
            highFactor = highFactor * 2;
        end
    end

    if ~found_tight
        fprintf('No tight solution found up to imp_factor = %g\n', highFactor);
        best_sol = Inf;
        best_powers = [];
        best_data_rates = [];
        return
    end

    % -- 2) Binary search
    best_sol_global = best_sol;
    best_powers_global = best_powers;
    best_data_rates_global = best_data_rates;

    for iter = 1:20
        if (highFactor - lowFactor) < 1e-5
            break
        end
        midFactor = 0.5*(lowFactor + highFactor);
        [test_sol_mid, test_pw_mid, test_dr_mid, c_tight_mid] = solve_for_imp_factor(U, H, min_rate, permutations, midFactor);

        if c_tight_mid
            highFactor = midFactor;
            if test_sol_mid < best_sol_global
                best_sol_global = test_sol_mid;
                best_powers_global = test_pw_mid;
                best_data_rates_global = test_dr_mid;
            end
        else
            lowFactor = midFactor;
        end
    end

    best_sol = best_sol_global;
    best_powers = best_powers_global;
    best_data_rates = best_data_rates_global;

    fprintf('Binary search done. imp_factor in [%.5f, %.5f].\n', lowFactor, highFactor);
    fprintf('Best objective value: %.4f\n', best_sol);
end

% Separate function for solving with a specific importance factor
function [best_sol_loc, best_pw_loc, best_dr_loc, all_c_tight] = solve_for_imp_factor(U, H, min_rate, permutations, imp_factor)
    best_sol_loc = Inf;
    best_pw_loc = [];
    best_dr_loc = [];
    all_c_tight = false;

    % Loop over all permutations
    for pIdx = 1:size(permutations,1)
        pi_ = permutations(pIdx, :);

        cvx_solver mosek
        cvx_begin quiet
            variable Rxx(U,U) nonnegative
            variable B_rates(U,U)  % Changed from Bvar to B_rates
            variable Cvar(U,1)  % was C

            constraints = [];

            for i = 1:U
                % S1 = {(i, j) for j=1..U}
                S1 = [repmat(i,1,U); 1:U]';
                % S2 = {(j, i) for j != i}
                S2 = [];
                for jj = 1:U
                    if jj ~= i
                        S2 = [S2; jj, i];
                    end
                end

                S1S2 = [S1; S2];
                % S3 = all pairs (j,k) not in S1S2
                S3 = [];
                for jj = 1:U
                    for kk = 1:U
                        in_S1S2 = false;
                        for row_ = 1:size(S1S2,1)
                            if (S1S2(row_,1)==jj) && (S1S2(row_,2)==kk)
                                in_S1S2 = true;
                                break
                            end
                        end
                        if ~in_S1S2
                            S3 = [S3; jj, kk];
                        end
                    end
                end

                % Sort S1+S2 by user decoding order pi_
                S1S2_list = [S1; S2];
                decoding_priority = @(x) find(pi_ == x);
                decodeOrderRank = zeros(size(S1S2_list,1),2);
                for row_ = 1:size(S1S2_list,1)
                    jVal = S1S2_list(row_,1);
                    kVal = S1S2_list(row_,2);
                    decodeOrderRank(row_,1) = decoding_priority(jVal);
                    decodeOrderRank(row_,2) = kVal;
                end
                [~, orderIdx] = sortrows(decodeOrderRank,[1,2]);
                decoding_order = S1S2_list(orderIdx,:);

                % Cvar(i) >= 0
                constraints = [constraints, Cvar(i) >= 0];

                % Cvar(i) <= 0.5 * (1/log(2)) * log_det( I + ... )
                if ~isempty(S3)
                    H_unimportant = H(i, S3(:,1));
                    diagRxx = [];
                    for rr = 1:size(S3,1)
                        diagRxx = [diagRxx; Rxx(S3(rr,1), S3(rr,2))];
                    end
                    constraints = [ constraints, ...
                        Cvar(i) <= 0.5 * (1/log(2)) * ...
                                  log_det( eye(1) + H_unimportant * diag(diagRxx) * H_unimportant' ) ];
                else
                    % If S3 empty => forces Cvar(i)=0
                    constraints = [constraints, Cvar(i) <= 0];
                end

                % Decoding-order constraints
                cumulative_sum = Cvar(i);
                for idx_ = 1:size(decoding_order,1)
                    j_ = decoding_order(idx_,1);
                    k_ = decoding_order(idx_,2);

                    usedPairs = [S3; decoding_order(1:idx_, :)];
                    H_eff = H(i, usedPairs(:,1));
                    diagRxx_eff = [];
                    for rr = 1:size(usedPairs,1)
                        diagRxx_eff = [diagRxx_eff; Rxx(usedPairs(rr,1), usedPairs(rr,2))];
                    end
                    
                    temp = B_rates(j_,k_);  % Changed from Bvar to B_rates
                    constraints = [constraints, ...
                        cumulative_sum + temp <= ...
                        0.5 * (1/log(2)) * log_det( eye(1) + H_eff * diag(diagRxx_eff) * H_eff' ) ];
                    
                    cumulative_sum = cumulative_sum + B_rates(j_,k_);  % Changed from Bvar to B_rates
                end

                % sum of B_rates(i,:) == min_rate(i)
                row_sum = 0;
                for k = 1:U
                    row_sum = row_sum + B_rates(i,k);  % Changed from Bvar to B_rates
                end
                constraints = [constraints, row_sum == min_rate(i)];
                
                % B_rates(i,:) >= 0
                for k = 1:U
                    constraints = [constraints, B_rates(i,k) >= 0];  % Changed from Bvar to B_rates
                end
            end

            % Objective
            objective = sum(sum(Rxx)) - imp_factor*sum(Cvar);
            minimize(objective)
            subject to
                constraints
        cvx_end

        % If solved, see if all Cvar(i) constraints were tight
        if contains(cvx_status,'Solved','IgnoreCase',true)
            if cvx_optval < best_sol_loc
                all_tight_check = true;
                for i = 1:U
                    S1 = [repmat(i,1,U); 1:U]';
                    S2 = [];
                    for jj = 1:U
                        if jj ~= i
                            S2 = [S2; jj, i];
                        end
                    end
                    S1S2 = [S1; S2];
                    S3 = [];
                    for jj = 1:U
                        for kk = 1:U
                            in_S1S2 = false;
                            for row_ = 1:size(S1S2,1)
                                if (S1S2(row_,1)==jj) && (S1S2(row_,2)==kk)
                                    in_S1S2 = true;
                                    break
                                end
                            end
                            if ~in_S1S2
                                S3 = [S3; jj, kk];
                            end
                        end
                    end

                    if ~isempty(S3)
                        H_unimportant = H(i, S3(:,1));
                        diagRxx_vals = zeros(size(S3,1),1);
                        for rr = 1:size(S3,1)
                            diagRxx_vals(rr) = Rxx(S3(rr,1), S3(rr,2));
                        end
                        rhs_val = 0.5 * (1/log(2)) * log_det( eye(1) + ...
                                              H_unimportant*diag(diagRxx_vals)*H_unimportant' );
                    else
                        rhs_val = 0;
                    end
                    lhs_val = Cvar(i);
                    if abs(lhs_val - rhs_val) > 1e-4
                        all_tight_check = false;
                        break
                    end
                end
                if all_tight_check
                    best_sol_loc = cvx_optval;
                    best_pw_loc = Rxx;
                    % Calculate sum_B_rates element by element
                    sum_B_rates = zeros(U,1);
                    for i = 1:U
                        for j = 1:U
                            sum_B_rates(i) = sum_B_rates(i) + B_rates(i,j);  % Changed from Bvar to B_rates
                        end
                    end
                    best_dr_loc = sum_B_rates.';
                    all_c_tight = true;
                end
            end
        end
    end
end