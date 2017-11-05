function [Tbeta, time_saif, ite] = saif(X, y, lambda, optTol, verbose)

maxIter = 100000000;
mexC_P = 1;

[en, ep] = size(X);

tstart= tic;
BK = 3;
Xy = X'*y;
absXYl = abs(Xy/lambda);


[~, im] = sort(absXYl, 'descend');
im = im(1:BK);
aX = X(:,im);
aIdx = [im];

rIdx = [1:ep]';
rIdx(im) = [];
rX = X(:, rIdx);
rN = length(rIdx);
rbeta = zeros(rN,1);

beta = LassoShoot_Cpp(aX,y,lambda, 100*optTol);


rmN = 0;
adN = BK;
rmLab = [];

mmXY = max(absXYl);
mdXY = median(absXYl);

ite = 0;
XX = sum(X.*X, 1)';
Xnorm = sqrt(XX);

AConst = 5/8;  
hPerameter = log(mmXY + mdXY)^1.4*log(ep)*AConst;
MaxNX = ceil(hPerameter); 
iNAdd = ceil(hPerameter*6/5);  
HiNAdd = ceil(iNAdd/2);
gfact = 1/(2*max(absXYl));


if gfact > 0.01
    gfact = 1; 
end

gincf = 10;
gfactor = gfact;


SConst =  6; 
STEP = ceil(SConst*log(mmXY + mdXY)^1.6*log(ep));  
gSTEP= ceil(STEP*3/2);  

IsEnhanced = 1;
IsAdd = 1;

yOlam = y/lambda;
mextime = 0;
XXProdTime = 0;
oldgap2 = 1000000;
NoProgress = 0;

while  1 
            theta = (y - aX*beta)/lambda;
            aXth = aX'*theta;
            ub = 1/max(abs(aXth));
            lb = -ub;
            alpha = min(max(lb,theta'*y/(lambda*theta'*theta)), ub);
            theta = alpha*theta;
 
            primV = objective(aX, y, beta, lambda);
            dualV = dual_objective(y, theta, lambda);
            gap2 = primV - dualV;
            if gap2 <0
                gap2=0;
            end
            gapo = gap2;
            
            if verbose == 1
                gapo
            end
            
            gap2 = gap2*gfact;
            rr = sqrt(2*gap2)/lambda;

%%%%%%%%%%%%%%%%%%%%%% DEL %%%%%%%%%%%%%%%%%%%%%%%%           
            Xc = aX'*theta;
            xs = Xnorm(aIdx);
            L = Xc - xs*rr;
            H = Xc + xs*rr;
            
            rmlab = (H<1)&(L>-1);
            rmN = rmN + sum(rmlab);
            stlab = ~rmlab;
 
            
            if sum(stlab) ~=0
                rmLab = [ rmLab; aIdx(rmlab)];
                rX = [rX, aX(:, rmlab)];
                rIdx = [rIdx; aIdx(rmlab)];
                rbeta = [rbeta; beta(rmlab)];

                beta = beta(stlab);
                aIdx = aIdx(stlab);
                aX = aX(:, stlab);
            end
     
           rXth = rX'*theta;
           rXnorm = Xnorm(rIdx);
           
           if ~IsAdd
                if gap2 < optTol
                   break; 
                else
                    if gap2 > oldgap2
                       break; 
                    end
                end

                if mexC_P == 1
                                mstart = tic;
                                XX2 = aX'*aX;
                                Xy2 = aX'*y;
                                p = size(aX,2);
                                XXProdTime = XXProdTime + toc(mstart);
                                mstart = tic;
                                beta= mexStepLassoShoot(XX2, Xy2, beta, p, lambda, gSTEP);
                                mextime = mextime + toc(mstart);
                                ite = ite + gSTEP;
                else
                                mstart = tic;
                                xx = XX(aIdx);
                                aXy = Xy(aIdx);
                                aXb = aX*beta;
                                XXProdTime = XXProdTime + toc(mstart);
                                mstart = tic;
                                beta = mexStepLassoShootN(aX, xx, aXy, beta, aXb, lambda, gSTEP);
                                mextime = mextime + toc(mstart); 
                end
                oldgap2 = gap2;
                continue;
           end
        
 %%%%%%%%%%%%%%%%%%%%%%% STOP ADD? %%%%%%%%%%%%%%%%%%%     
         rXr = rXnorm*rr;
         mxth = max(abs(rXth) + rXr);
         if mxth < 1  
            if gfact < 1
                gfact = min(gfact*gincf, 1); 
                if gfact >0.55
                    gfact = 1;
                end
            else
                IsAdd = 0;
                %disp('no add');
                continue;
            end
         end   
 
 %%%%%%%%%%%%%%%%%%%%%  ADD  %%%%%%%%%%%%%%%%%%%%%%%%%%%          
          arXth = abs(rXth);
          maxlist = arXth + rXr;
   
          for iad = 1:iNAdd  
               [marXth, imx] = max(arXth);
               iminx = abs(marXth - rXr(imx));
               tt = iminx - maxlist;
               if (sum(tt < 0) < MaxNX)
                    im = imx;
                    [aX, aIdx, beta, rX, rIdx, rbeta] = ADDOne(im, aX, aIdx, beta, ...
                                                            rX, rIdx, rbeta,y, lambda);
                    arXth(im) = []; 
                    maxlist(im) = [];
                    rXr(im) = [];
                    adN = adN + 1;
               else
                  break; 
               end
          end
           
           if (iad == 1)&(NoProgress==1)
                [~, imx] = sort(maxlist, 'descend');
                im = imx(1:HiNAdd);
                aX = [aX, rX(:, im)];
                aIdx = [aIdx; rIdx(im)];
                beta = [beta; zeros(HiNAdd,1)];
                rIdx(im) = [];
                rX(:,im) = [];
                adN = adN + HiNAdd;
           end

                mstart = tic;
                XX2 = aX'*aX;
                Xy2 = aX'*y;
                p = size(aX,2);
                XXProdTime = XXProdTime + toc(mstart);
                mstart = tic;
                oldbeta = beta;
                beta= mexStepLassoShoot(XX2, Xy2, beta, p, lambda, STEP);
                if norm(oldbeta - beta) < 1.0E-8
                    NoProgress = 1;
                else
                    NoProgress = 0;
                end
                
                mextime = mextime + toc(mstart);
             
            ite = ite + STEP;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tbeta = zeros(ep,1);
Tbeta(aIdx) = beta;
time_saif = toc(tstart);

if verbose == 1
    disp('**************** infor ******************');
    rmN
    adN
    N_nonezero = sum(Tbeta~=0)
    dgap = primal_dual_gap(X,y,lambda, Tbeta)
    %mextime
    %XXProdTime
    %gfactor
    time_saif
end

    
