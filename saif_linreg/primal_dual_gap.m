function g = primal_dual_gap(X,y,lambda, beta)
    
theta = (y - X*beta)/lambda;
aXth = X'*theta;
ub = 1/max(abs(aXth));
lb = -ub;
alpha = min(max(lb,theta'*y/(lambda*theta'*theta)), ub);
theta = alpha*theta;

primV = objective(X, y, beta, lambda);
dualV = dual_objective(y, theta, lambda);
g = primV - dualV;
