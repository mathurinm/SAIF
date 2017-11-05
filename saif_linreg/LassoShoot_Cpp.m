function [beta, time_shooting, m, primV, dualV, gap2, theta] = LassoShoot_Cpp(X,y,lambda, optTol)

sttime = tic;

maxIter = 100000000;
[n, p] = size(X);

m = 0;
XX2 = X'*X;
Xy2 = X'*y;
XX = diag(XX2);
mex_P = 1;
beta = (XX2 + lambda*eye(p))\(Xy2);
Xb = X*beta;
STEP = 40*log(p);
%fprintf('\n');
while m < maxIter
    
           beta_old = beta;
           if mex_P ==1  
              beta= mexStepLassoShoot(XX2', Xy2, beta, p, lambda, STEP);
           else
              beta = mexStepLassoShootN(X, XX, Xy2, beta, Xb, lambda, STEP); 
           end
           m = m + STEP;
            
            primV= objective(X, y, beta, lambda);
            theta = (y - X*beta)/lambda;
            Xth = X'*theta;
           % tt = sum(abs(D'),2);
            ub = 1/max(abs(Xth));
            lb = -ub;
            alpha = min(max(lb,theta'*y/(lambda*theta'*theta)), ub);

            theta = alpha*theta;
            dualV = dual_objective(y, theta, lambda);
            gap2 = primV - dualV;
            %fprintf('%f ', gap2);
            if gap2 < optTol
                   break; 
            end
end
cpp_gap = gap2;
time_shooting = toc(sttime);

