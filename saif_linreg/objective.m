
function obj = objective(X, Y,  b, lambda)
    obj = 0.5*sum((Y - X*b).^2) + lambda*norm(b,1);


