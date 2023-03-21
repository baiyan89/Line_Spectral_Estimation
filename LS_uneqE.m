function  [Ks, P]=LS_uneqE(Y,E)
%% least squares
M=size(Y,1);
NP=size(E,1);
P =cell(NP, 1);
Ks=zeros(NP,1);

for j=1:NP 
    e=E{j};
    Ks(j)=numel(e);   
    
    A=exp(-1i*pi*(0:M-1)'*sind(e));
    xe=pinv(A)*Y;
%     xe=(A'*A)\(A'*Y); 
    P{j}=xe;
end
end