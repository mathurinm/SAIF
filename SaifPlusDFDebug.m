function [Tbeta, time_saif, ite] = SaifPlusDFDebug(X, y, lambda, optTol, verbose)
%optTol = 10E-6;
%lammax = max(abs(X'*Y));
%lg_data = {};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    SAIF   algorithm                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxIter = 100000000;
mexC_P = 1;

[en, ep] = size(X);

tstart= tic;
BK = 3;
%y = Y;
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
adN = 0;
rmLab = [];

%absXYl = abs(X'*y);
mmXY = max(absXYl);
mdXY = median(absXYl);

ite = 0;
XX = sum(X.*X, 1)';
Xnorm = sqrt(XX);

AConst = 5/8;   %fypical 0.2   %% typical 5/8;
hPerameter = log(mmXY + mdXY)^1.4*log(ep)*AConst %; %%bcancer 1;
MaxNX = ceil(hPerameter); 
iNAdd = ceil(hPerameter*6/5);  %%typical 6/5; 0.2

gfact = 1/(2*max(absXYl))
%gfact = min(exp(-(mmXY + mdXY)/76), 1)
%gfact = min(exp(-(mmXY + mdXY)/140), 1)%1E-4%

if gfact > 0.01
    gfact = 1; 
end
% elseif gfact < 5E-4
%     gfact = 5E-4;
% end
gincf = 10
gfactor = gfact
%pause;

SConst =  6; %1%;  ftypical 3 %% typical 6;
STEP = ceil(SConst*log(mmXY + mdXY)^1.6*log(ep));  %%bcancer 1.6;
gSTEP= ceil(STEP*3/2);  %% bcancer 3/2;

IsEnhanced = 1;
IsAdd = 1;

%q2 = (1/lambda - 1/lamMax)*y;
yOlam = y/lambda;
mextime = 0;
XXProdTime = 0;
%fact = 0.00001
%gfact = fact;
oldgap2 = 1000000;
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
            gap2 = gap2*gfact;
          
            %zct2 = theta;
            rr = sqrt(2*gap2)/lambda;
            %rro = sqrt(2*gapo)/lambda;
            
            
%%%%%%%%%%%%%%%%%%%%%% DEL %%%%%%%%%%%%%%%%%%%%%%%%           
            Xc = aX'*theta;
            xs = Xnorm(aIdx);%sqrt(sum(aX.*aX,1));
            L = Xc - xs*rr;
            H = Xc + xs*rr;
            
            rmlab = (H<1)&(L>-1);
            rmN = rmN + sum(rmlab);
            stlab = ~rmlab;
            
            %rmN
            %adN
            %rmlab
            
            if sum(stlab) ~=0
                rmLab = [ rmLab; aIdx(rmlab)];
                rX = [rX, aX(:, rmlab)];
                %rmft = aIdx(rmlab)
                rIdx = [rIdx; aIdx(rmlab)];
                rbeta = [rbeta; beta(rmlab)];

                beta = beta(stlab);
                aIdx = aIdx(stlab);
                aX = aX(:, stlab);
            end
     
           rXth = rX'*theta;
           rXnorm = Xnorm(rIdx);%sqrt(sum(rX.*rX,1))'; 
           
           if ~IsAdd
                %%%% how to stop  to determine the DEL opts?
                if gap2 < optTol
                    %disp('gap2');
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
        
 %%%%%%%%%%%%%%%%%%%%%%% STOP?? %%%%%%%%%%%%%%%%%%%     
         rXr = rXnorm*rr;
         mxth = max(abs(rXth) + rXr);
         if mxth < 1  %%%?????????
            if gfact < 1
                gfact = min(gfact*gincf, 1); %g_inf
                if gfact >0.55
                    gfact = 1;
                end
            else
                IsAdd = 0;
                disp('no add');
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
           

             if mexC_P == 1
                mstart = tic;
                XX2 = aX'*aX;
                Xy2 = aX'*y;
                p = size(aX,2);
                XXProdTime = XXProdTime + toc(mstart);
                mstart = tic;
                beta= mexStepLassoShoot(XX2, Xy2, beta, p, lambda, STEP);
                mextime = mextime + toc(mstart);
             else
                mstart = tic;
                xx = XX(aIdx);
                aXy = Xy(aIdx);
                aXb = aX*beta;
                XXProdTime = XXProdTime + toc(mstart);
                mstart = tic;
                beta = mexStepLassoShootN(aX, xx, aXy, beta, aXb, lambda, STEP);
                mextime = mextime + toc(mstart);  
             end
             
            ite = ite + STEP;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('**************** statistics ******************');
Tbeta = zeros(ep,1);
Tbeta(aIdx) = beta;
time_saif = toc(tstart);

if verbose == 1
    rmN
    adN
    N_nonezero = sum(Tbeta~=0)
    dgap = primal_dual_gap(X,y,lambda, Tbeta)
    mextime
    XXProdTime
    gfactor
    time_saif
end

    
%save('2saif_log_data.mat','lg_data');