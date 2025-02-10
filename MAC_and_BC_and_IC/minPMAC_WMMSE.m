function [FEAS_FLAG, bu_a, info] = minPMAC_WMMSE(H, bu_min, w, cb)

    % Start timer
    tic;

    % Initialize dimensions
    [Ly, U, N] = size(H);
    Lx = ones(1, U);  % Assuming single transmit antenna per user
    max_iter = 100;
    tol = 1e-4;
    sigma2 = 1e-3;  % Noise variance (adjust as needed)

    % Initialize transmit covariance matrices Rxx (U x N cells)
    Rxx = cell(U, N);
    for u = 1:U
        for n = 1:N
            Rxx{u, n} = (1/Lx(u)) * eye(Lx(u));  % Initial power allocation
        end
    end

    % Initialize receiver filters and weights
    Gu = cell(U, N);
    Wu = cell(U, N);

    % Precompute channel matrices
    H_un = cell(U, N);
    for u = 1:U
        for n = 1:N
            H_un{u, n} = H(:, u, n);
        end
    end

    % Start main iteration loop
    for iter = 1:max_iter
        % Store previous Rxx for convergence check
        Rxx_prev = Rxx;

        % Step 1: Update receiver filters Gu
        for u = 1:U
            for n = 1:N
                % Compute interference plus noise covariance
                Sigma = sigma2 * eye(Ly);
                for i = 1:U
                    Sigma = Sigma + H_un{i, n} * Rxx{i, n} * H_un{i, n}';
                end
                % Compute MMSE receiver
                Gu{u, n} = Rxx{u, n} * H_un{u, n}' / Sigma;
            end
        end

        % Step 2: Update weights Wu
        for u = 1:U
            for n = 1:N
                % Compute MSE
                E_un = eye(Lx(u)) - Gu{u, n}' * H_un{u, n} * Rxx{u, n};
                % Update weight
                Wu{u, n} = inv(E_un);
            end
        end

        % Step 3: Update transmit covariance matrices Rxx
        cvx_begin quiet
            % cvx_solver mosek
            variables Rxx_cell{U, N}
            expressions total_power
            total_power = 0;
            for u = 1:U
                for n = 1:N
                    Rxx_cell{u, n} = semidefinite(Lx(u));
                    total_power = total_power + w(u) * trace(Rxx_cell{u, n});
                end
            end
            minimize(total_power)
            subject to
                for u = 1:U
                    for n = 1:N
                        % Rate constraint converted to MSE constraint
                        % Compute effective channel
                        Heff = H_un{u, n}' * Gu{u, n};
                        % Compute SINR
                        SINR = trace(Heff' * Wu{u, n} * Heff * Rxx_cell{u, n});
                        % Convert minimum rate to SINR threshold
                        gamma_min = 2^(cb * bu_min(u) / N) - 1;
                        SINR >= gamma_min;
                    end
                end
        cvx_end

        % Assign optimized Rxx
        for u = 1:U
            for n = 1:N
                Rxx{u, n} = Rxx_cell{u, n};
            end
        end

        % Check convergence
        diff_norm = 0;
        total_norm = 0;
        for u = 1:U
            for n = 1:N
                diff_norm = diff_norm + norm(Rxx{u, n} - Rxx_prev{u, n}, 'fro')^2;
                total_norm = total_norm + norm(Rxx_prev{u, n}, 'fro')^2;
            end
        end
        if sqrt(diff_norm / total_norm) < tol
            disp(['Converged at iteration ', num2str(iter)]);
            break;
        end
    end

    % Compute achieved rates and energies
    b = zeros(U, N);
    for u = 1:U
        for n = 1:N
            % Compute interference plus noise covariance
            Sigma = sigma2 * eye(Ly);
            for i = 1:U
                if i ~= u
                    Sigma = Sigma + H_un{i, n} * Rxx{i, n} * H_un{i, n}';
                end
            end
            % Compute SINR
            Numerator = H_un{u, n}' * Rxx{u, n} * H_un{u, n};
            Denominator = Sigma;
            SINR = real(trace(Numerator)) / real(trace(Denominator));
            % Compute rate
            b(u, n) = (1 / cb) * log2(1 + SINR);
        end
    end

    % Compute average rates and energies
    bu_a = sum(b, 2)';
    Eu_a = zeros(1, U);
    for u = 1:U
        for n = 1:N
            Eu_a(u) = Eu_a(u) + trace(Rxx{u, n});
        end
    end
    Eu_a = Eu_a / N;

    % Check feasibility
    FEAS_FLAG = all(bu_a >= bu_min);

    % Store information
    info.H = H;
    info.bu_min = bu_min;
    info.w = w;
    info.cb = cb;
    info.Rxx = Rxx;
    info.Gu = Gu;
    info.Wu = Wu;
    info.b = b;
    info.bu_a = bu_a;
    info.Eu_a = Eu_a;
    info.FEAS_FLAG = FEAS_FLAG;

    elapsedTime = toc;
    disp(['Elapsed time for WMMSE minPMAC: ', num2str(elapsedTime), ' seconds']);
end