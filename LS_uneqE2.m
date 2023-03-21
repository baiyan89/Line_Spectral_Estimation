function  [Ks, P]=LS_uneqE2(Y,Omega,E)
%% least squares
M=size(Y,1);
NP=size(E,1);
P =cell(NP, 1);
Ks=zeros(NP,1);

for j=1:NP 
    e=E{j};
    Ks(j)=numel(e);   
    
    A=exp(-1i*pi*Omega*sind(e));
    xe=pinv(A)*Y;
%     xe=(A'*A)\(A'*Y); 
    P{j}=xe;
end
end