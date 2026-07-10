# Lanczos-bidiagonalization-method-for-the-large-scale-total-least-squares-problem
Lanczos bidiagonalization method for the large-scale total least squares problem

% [fx,x,err] = PSLB_TLS(A,b,kmax,tol): PSLB algorithm for solving TLS problems (without restarted)

% [fx1,x,err1] = Restart_PSLB_TLS(A,b,int_it,out_it,tol,flag): a restarted block Lanczos bidiagonalization method for TLS problem

% [V, D, flag] = tls_eigs(varargin): This algorithm is built on the eigs.m, with only the stopping criterion adjusted

% [S, V, flag] = tls_svds(A,varargin): This algorithm is built on the svds.m, with only the stopping criterion adjusted

% [x, f_hist] = TLS_RQI(A, b, max_iter, tol): Rayleigh quotient iteration (RQI) iteration for Total Least Squares



%test:
% (i) Test the generic case
      
      A = sprand(m,n,0.005); 
      b = ones(size(A,1),1); b = b/norm(b);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      t = tic; [fx1,x,err1] = Restart_PSLB_TLS(A,b,50,200,1e-13,1); t = toc(t)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      opts.maxit = 200; opts.p = 50; opts.tol = 1e-13;
      t = tic; Ab = [A,b]'*[A,b]; Ab = 0.5*(Ab'+Ab); [V, fx, flag] = tls_eigs(Ab,1,'sr',opts); t = toc(t)
               x = -V(1:end-1)/V(end);
      t = tic; [S, V, flag] = tls_svds([A,b],1,'smallest',opts);tsvds = toc(t) 
                fx = S*S; x = -V(1:end-1)/V(end);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      t = tic; [x, f_hist] = TLS_RQI(A, b, 200, 1e-13); t = toc(t)


 % (ii) Test the generic case

n = 3000;
m = 10000;
Q1 = randn(n);
Q1(n, 1:2) = zeros(1, 2); 
[Q, ~] = qr(Q1);


D = (1:n)/1000;    
D(end) = 3;        
S = eye(n) + diag(D);
P = [S; zeros(m - n, n)] * Q';
A = P(:, 1:2999);
b = P(:, 3000);

[fx1,x,err1] = Restart_PSLB_TLS(A,b,50,200,1e-8,1);

[fx,x,err] = PSLB_TLS(A,b,n,1e-8);
