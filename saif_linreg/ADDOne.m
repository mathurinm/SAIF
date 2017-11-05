function  [aX, aIdx, beta, rX, rIdx, rbeta] = ADDOne(im, aX, aIdx, beta, rX, rIdx, rbeta, y, lambda)


thisX = rX(:,im);
th = y - aX*beta;
a =thisX'*thisX;
c = thisX'*th;

ca= c/a;
bti = sign(ca)*max([abs(ca)-(lambda/a), 0]);
beta = [beta; bti];
%beta = [beta; rbeta(im)];
%beta = [beta; 0]; % zeros(length(im),1)];
aX = [aX, thisX];
aIdx = [aIdx; rIdx(im)'];

rX(:,im) = [];
rIdx(im) = [];
rbeta(im) = [];

