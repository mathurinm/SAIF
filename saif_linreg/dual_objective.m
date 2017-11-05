function dobj = dual_objective(Y, theta, lambda)

thY = theta - Y/lambda;
dobj = Y'*Y/2 - thY'*thY*lambda^2/2;