clc

addpath('channel_files', 'FA4', 'MAC_and_BC', 'single_user', 'tmp', 'utils', ...
    'FA4/Wi-Fi_ChannelModels_B_D', 'MAC_and_BC/Subroutines', 'FA4/ML', ...
    'FA4/ML/outputs');

noise_var =  1.0/9800000;
U = 3;
N = 64;
Ly = 4;
Lx = [2, 2, 2];
% bmin = [253.6993  654.7796  500.3479];
bmin = bu_a_lin;
subsets = powerSet(1:U);
subsets = subsets(2:end);

H = reshape(H, [4, 2, 3, 64]);
H = H/sqrt(noise_var);

Rxx = cell(U, N);
tic;
cvx_begin quiet
    cvx_solver mosek
    variable Rxx(Lx(1), Lx(1), U, N) hermitian semidefinite
    variable b(U, N) nonnegative
    dual variable theta
    expressions rate_sum


    theta : sum(b, 2)' >= bmin;
    
     for subset = subsets
         for n = 1:N
             temp = zeros(Ly, Ly);
             for u = subset{:}
                temp = temp + H(:, :, u, n) * Rxx(:, :, u, n) * H(:, :,  u, n)';
             end
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
            obj = obj + trace(Rxx(:,:, u, n));
         end
     end
     minimize (obj)
cvx_end
Rxx_old = Rxx;
elapsedTime = toc;
fprintf('Elapsed time for convex optimization in reformulated routine: %.2f seconds\n', elapsedTime);


% Finding actually achieved rates at vertices of capacity region

pi_inv = [2 1 3];
b_achieved_vertex1 = zeros(U, N);


for n = 1:N
    for u = 1:U
        mat1 = 0;
        for i = u:U
            mat1 = mat1 + H(:, :,  pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:,:, pi_inv(i),n)';
        end
        mat2 = 0;
        for i = u+1:U
            mat2 = mat2 + H(:, :,  pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:,:, pi_inv(i),n)';
        end

        b_achieved_vertex1(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
    end
end


% Other vertex

pi_inv = [2, 3, 1];

[~, pi] = sort(pi_inv);
b_achieved_vertex2 = zeros(U, N);


for n = 1:N
    for u = 1:U
        mat1 = 0;
        for i = u:U
            mat1 = mat1 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end
        mat2 = 0;
        for i = u+1:U
            mat2 = mat2 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end

        b_achieved_vertex2(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
    end
end


pi_inv = [1 2 3];

[~, pi] = sort(pi_inv);
b_achieved_vertex3 = zeros(U, N);


for n = 1:N
    for u = 1:U
        mat1 = 0;
        for i = u:U
            mat1 = mat1 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end
        mat2 = 0;
        for i = u+1:U
            mat2 = mat2 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end

        b_achieved_vertex3(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
    end
end


pi_inv = [1 3 2];

[~, pi] = sort(pi_inv);
b_achieved_vertex4 = zeros(U, N);


for n = 1:N
    for u = 1:U
        mat1 = 0;
        for i = u:U
            mat1 = mat1 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end
        mat2 = 0;
        for i = u+1:U
            mat2 = mat2 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end

        b_achieved_vertex4(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
    end
end


pi_inv = [3 1 2];

[~, pi] = sort(pi_inv);
b_achieved_vertex5 = zeros(U, N);


for n = 1:N
    for u = 1:U
        mat1 = 0;
        for i = u:U
            mat1 = mat1 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end
        mat2 = 0;
        for i = u+1:U
            mat2 = mat2 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end

        b_achieved_vertex5(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
    end
end


pi_inv = [3 2 1];

[~, pi] = sort(pi_inv);
b_achieved_vertex6 = zeros(U, N);


for n = 1:N
    for u = 1:U
        mat1 = 0;
        for i = u:U
            mat1 = mat1 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end
        mat2 = 0;
        for i = u+1:U
            mat2 = mat2 + H(:, :, pi_inv(i), n)*Rxx(:,:,pi_inv(i), n)*H(:, :,pi_inv(i),n)';
        end

        b_achieved_vertex6(pi_inv(u), n) = log2(real(det(eye(Ly) + mat1))) - log2(real(det(eye(Ly) + mat2)));
    end
end


% Functions

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