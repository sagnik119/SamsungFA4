% the channel matrix is K*N by m, m is the number of transmitter's antennas,
% K is the number of users and N is the number od receiver's antennas for
% each users
% Rx is the m by m input covariance matrix 
% Rz is the K*N by K*N worst case noise minimizing           
% 0.5 * log(det(H*Rx*H'+Rz)/det(Rz)), Rz has its N by N diagonal sub-blocks  
% equal to identity matrix
% dual_gap is the duality gap, set to 1e-6 in wcnoise
% nerr is Newton's method acceptable error, set to 1e-4 in wcnoise
% sumRate is the rate-sum/real-dimension

function [Rz, sumRate] = wcnoise(Rx, H, N, dual_gap, nerr)

switch nargin
    case 3
        dual_gap = 1e-6;
        nerr = 1e-4;
    case 4
        nerr = 1e-4;
end

[n,m] = size(H);
K = n / N;

A = zeros(n,n,n*(n+1)/2);
for i = 1:n
   A(i,i,i) = 1;
end

count = n+1;

for i = 2:n
   for j = 1:i - 1
      A(i,j,count) = 1;
      A(j,i,count) = 1;
      count = count+1;
   end
end

map = zeros(n,n);
for i = 1:K
   map((i-1) * N + 1:i * N,(i-1) * N + 1:i * N) = ones(N,N);
end




NT_max_it = 1000;                              % Maximum number of Newton's 
                                               % method iterations
%dual_gap = 1e-6;
mu = 10;                                       % step size for t
alpha = 0.001;                                 % back tracking line search parameters
beta = 0.5;

count = 1;
%n%nerr = 1e-4;                                   % acceptable error for inner loop 
                                               % Newton's method

v_0 = zeros(n*(n+1)/2,1);                      % Strictly feasible point;
v_0(1:n) = 0.5 *  ones(n,1);
v = v_0;
t =1;                                         
l_v = 1;                                       % lambda(v) for newton's method termination



while (1+n)/t > dual_gap
   t =  t * mu;
   l_v = 1;
   count = 1;
   while l_v > nerr & count < NT_max_it 
      
      f_val = 0;                               % calculating function value
      Rz = zeros(n,n);
      Rzprime = zeros(n,n);
      
      for i = 1:n*(n+1)/2                      % computing Rz
         Rz = Rz + v(i) * A(:,:,i);
      end
      
      for i = 1:K
          Rzprime((i-1) * N + 1:i * N,(i-1) * N + 1:i * N) = Rz((i-1) * N + 1:i * N,(i-1) * N + 1:i * N);
      end
      
      f_val = t * log(det(H * Rx * H' + Rz)) - (t + 1) * log(det(Rz)) - log(det(eye(n) - Rzprime));
      
      S = inv(H * Rx * H' + Rz);
      Q = inv(eye(n) - Rzprime);
      Rz_inv = inv(Rz);
      g = zeros(n*(n+1)/2,1);
      h = zeros(n*(n+1)/2,n*(n+1)/2);
      
      for i = 1:n*(n+1)/2
         g(i) = t * trace(A(:,:,i) * S) - (t + 1) * trace(A(:,:,i) * Rz_inv)...
              + (sum(sum(A(:,:,i) .* map)) ~= 0) * trace(A(:,:,i) * Q) ;  % gradient
      end
      
      
      for i = 1:n*(n+1)/2
         for j = 1:n*(n+1)/2
            h(i,j) = -t * trace(A(:,:,i) * S * A(:,:,j) * S) + (t + 1) * trace(A(:,:,i) * Rz_inv * A(:,:,j) * Rz_inv)...
                     +(sum(sum(A(:,:,i) .* map)) ~= 0) * (sum(sum(A(:,:,j) .* map)) ~= 0) * trace(A(:,:,i) * Q * A(:,:,j) * Q); %hessian
         end
      end
      
         
      dv = -h\g;                               % search direction
      
      
      s = 1;                                   % checking v = v+s*dx feasible 
                                               % and also back tracking algorithm
   
      v_new = v + s * dv;
      
      f_new = 0;
      
      Rz_new = zeros(n,n);
      
      for i = 1:n*(n+1)/2
         Rz_new = Rz_new + v_new(i) * A(:,:,i);
      end
      
      for i = 1:K
          Rzprime((i-1) * N + 1:i * N,(i-1) * N + 1:i * N) = Rz_new((i-1) * N + 1:i * N,(i-1) * N + 1:i * N);
      end
      
      f_new = t * log(det(H * Rx * H' + Rz_new)) - (t + 1) * log(det(Rz_new)) - log(det(eye(n) - Rzprime));
      

      
      feas_check = 1;
      if real(eig(Rz_new)) > zeros(n,1)
          feas_check = 1;
      else
          feas_check = 0;
      end
     
      if real(eig(eye(n) - Rzprime)) > zeros(n,1)
          feas_check = 1;
      else
          feas_check = 0;
      end
     
      feas_check = feas_check * (f_new < f_val + alpha * s * g' * dv);

   
      while feas_check ~= 1                
         s = s * beta ;
         
         v_new = v + s * dv;
         f_new = 0;
      
         Rz_new = zeros(n,n);
      
         for i = 1:n*(n+1)/2
            Rz_new = Rz_new + v_new(i) * A(:,:,i);
         end
      
         for i = 1:K
             Rzprime((i-1) * N + 1:i * N,(i-1) * N + 1:i * N) = Rz_new((i-1) * N + 1:i * N,(i-1) * N + 1:i * N);
         end
      
         f_new = t * real(log(det(H * Rx * H' + Rz_new))) - (t + 1) * real(log(det(Rz_new)) - log(det(eye(n) - Rzprime)));
      

      
         feas_check = 1;
         if real(eig(Rz_new)) > zeros(n,1)
             feas_check = 1;
         else
             feas_check = 0;
         end
         if real(eig(eye(n) - Rzprime)) > zeros(n,1)
             feas_check = 1;
         else
             feas_check = 0;
         end
     
         
         feas_check = feas_check * (f_new < f_val + alpha * s * g' * dv);



      end   
      
      v = v + s * dv;                         % update v
      l_v = -g'*dv ;                           % lambda(v)^2 for Newton's method
      count = count + 1;                      % number of Newtons method iterations
   end
  
end

Rz = zeros(n,n);
      
for i = 1:n*(n+1)/2
   Rz = Rz + v(i) * A(:,:,i);
end

for i = 1:K
    Rz((i-1) * N + 1:i * N,(i-1) * N + 1:i * N) = eye(N);
end
sumRate = 0.5 * log2(det(H * Rx * H' + Rz)/det(Rz));
