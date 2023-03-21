function [FV] = getME_uneqE2(E,P,Ks,Y,Omega) 
M=size(Y,1);
NP=size(E,1);
FV=zeros(NP,2);
FV(:,1)=Ks;
% sigma=1;

for j=1:NP
    thetaEst=E{j};
        
    A=exp(-1i*pi*Omega*sin(thetaEst*pi/180));
    FV(j,2)=norm(Y-A*P{j,1},'fro');    
    
%     r=Y-A*P{j};
%     fvs=1-sum(exp(-r.*conj(r)/2./sigma/sigma)/M)';
%     FV(j,2)=mean(fvs);
end
end