clear all; clc; close all; rng('shuffle');
%% settings
Methods = {'MVESA'};
numMethods = length(Methods) ;
numTrials =10;

theta=[-67.7  -63.3 10.9  ]; %-50.2 -21.5 6.2 10.9 24.5 63.4 

K=length(theta);
N=20;
tmp=randperm(N-1)';
array_type='full'; %'arbitrary'

SNRs=10;
numSNRs=numel(SNRs);
%%
res_angle=zeros(5,numMethods); res_num=zeros(5,numMethods); 
for T=max(20, 1)  %0,5,20
    for M=15
        switch array_type
            case 'full'
                Mcal=[0:M-1]';
                
            case 'sparse'
                Mcal=[sort(tmp(1:M-1),'ascend'); N]-1;
                
            case 'arbitrary'
                Mcal=([1:M-1])+1*(rand(1,M-1)-0.5);
                Mcal=(M-1)*([0,Mcal]'/max(Mcal));
        end        
        
        A=exp(-1i*pi*Mcal*sind(theta));
        X=randn(K,T)+1i*randn(K,T);
        R= 1 + .2.*randn(K,T);  % magnitudes of the complex amplitudes
        X= R.*X;  % complex amplitudes
        Y0=A*X;
        
        all_Angle_errors=nan(numMethods,11);       
        t=1;
        
        for SNR=SNRs(1)
            Pn  = mean(mean(abs(Y0).^2))*10^(-SNR/10);       % noise power
            eps = sqrt(0.5*Pn).*(randn(M,T)+1i*randn(M,T));  % complex Gaussian noise

            if SNR==-10
                eps=0;
            end
            Y=Y0+eps;
            
            theta_results=cell(numMethods,numTrials);  pow_results=cell(numMethods,numTrials);      
            angle_errors=nan(numMethods,numTrials);
            f_errors=nan(numMethods,numTrials);
            num_errors1=nan(numMethods,numTrials);  num_errors2=nan(numMethods,numTrials);            
            Iters=zeros(numMethods,numTrials);   times=zeros(numMethods,numTrials);
            
            for methodIndex=1
                methodName = [Methods{methodIndex}];
                methodHandle = str2func(methodName) ;               
                X0 = sprintf('method: %s  ---------',methodName);disp(X0)              
                
                for trial=1:numTrials
                    tic;
                    [xEst,thetaEst,XEst, iter]=methodHandle(Y,Mcal,K);
                    times(methodIndex, trial)=toc;
                    Iters(methodIndex, trial)=iter;
                    theta_results{methodIndex, trial}=thetaEst;
                    pow_results{methodIndex, trial}=xEst;
                    
                  %% 计算误差                  
                    C=abs(theta-thetaEst');
                    [Matching,~]=Hungarian (C);
                    [a,b]=find(Matching==1);%代表第b个任务，分给了第a个人
                  if numel(thetaEst)<K
                      num_errors1(methodIndex, trial)=numel(thetaEst)-K;
                  else                      
                      angle_errors(methodIndex, trial)=mean(abs(thetaEst(a)-theta));
                      f_errors(methodIndex, trial)=norm(sind(thetaEst(a))-sind(theta),2)/K;
                      num_errors1(methodIndex, trial)=numel(thetaEst)-K;
                  end      
                
                end

            end
            all_Angle_errors(:, t)=mean(angle_errors,2);
            t=t+1;
        end
    end
    angle_errors
end