function [power,theta_est,XEst, G]=MVESA(Y,Omega,K)   
rng('shuffle'); % rand、randi 和 randn 会在每次调用 rng 时生成不同的数字序列。

eta=10^(-6); tt=0;

[M,~]=size(Y); 
NP=30; % population size
Gmax=200;

%% Intialize maks and A
[E, ~] = initialize_uneqE(Y,NP);
[Ks,P] =LS_uneqE2(Y,Omega,E);
FV = getME_uneqE2(E,P,Ks,Y,Omega); 

[~,~,FrontNo,CrowdDis] = P_sort_new(M,NP,FV,1);

archiveE=[];
archiveP=[];
archiveFV=[];

old_s=ones(1,K); 
%% Main loop
G=1;
while G<=Gmax
    %% 打乱解的顺序
    sec =TournamentSelection(2,NP,FrontNo,-CrowdDis);%会丢失一些non-dominated解
   
    E0=E(sec,:);
    P0=P(sec);
    
    make_plot=0;
    if make_plot==1
        plot(FV(:,1),FV(:,2),'rx')
        hold on;
        xlabel('number of targets');
        ylabel('ME');
    end
    
    %% produce offspring solutions   
    Offspring_E = crossover2(E0,P0,M); 
    Offspring_E=mutation(Offspring_E, []);
    [Ks2,Q] =LS_uneqE2(Y,Omega,Offspring_E);
    FV2 = getME_uneqE2(Offspring_E,Q,Ks2,Y,Omega);
    
    
    if make_plot==1
        plot(FV2(:,1),FV2(:,2),'b+')
        hold on;
    end
      
    
    %% combination         
    P=[P; Q];    
    E=[E; Offspring_E];
    FV=[FV;FV2];
    

    %% get rid of repetitive solutions
    [FV,ia]=unique(FV,'rows');
    P=P(ia);
    E=E(ia);
    
    %% selection 
    [FV,loc]=unique(FV,'rows');
    P=P(loc);
    E=E(loc);
    [Next,non_index,FrontNo,CrowdDis] = P_sort_new(M,NP,FV,1); 
    nonFV=FV(non_index,:);
    nonE=E(non_index,:);
    nonP=P(non_index);    
      
    
    P=P(Next);    
    E=E(Next,:); 
    FV=FV(Next,:);
    
    if make_plot==1
        plot(FV(:,1),FV(:,2),'k.')
        hold off;
        legend('old','new','select')
        title(['G=',num2str(G)]);
        axis([1 8 0 30])
        pause(0.1)   
    end    
    
    %% update archive: 
    % 先保存每个K下的min f2对应的solution
    comFV=[archiveFV; nonFV];
    comP=[archiveP; nonP];
    comE=[archiveE; nonE];
    uKs=unique(comFV(:,1));%uKs=uKs(uKs<=M);
    len_u=numel(uKs);
    
    
    % reduce size of archive (<M)
    archiveE=cell(len_u,1);
    archiveP=cell(len_u,1);
    archiveFV=zeros(len_u,2);
    for u=1:len_u
        loc1=find(comFV(:,1)==uKs(u));
        [~, loc2]=min(comFV(loc1,2));
        archiveE(u,1)=comE(loc1(loc2(1)));
        archiveP(u,1)=comP(loc1(loc2(1)));
        archiveFV(u,:)=comFV(loc1(loc2(1)),:);
    end 
   
    
    % truncation: 优先去掉K+1解中power较小元素，看其fv是不是比K解要小，若小，则替换k解
    for count0= 2  : len_u  %
        [~,lock]=min(abs(archiveFV(:,1)-count0));
        e=archiveE{lock};
        power=sum(abs(archiveP{lock}),2);
        [~,sec]=sort(power,'descend');
        [~,diff_max]=max(abs(diff(power(sec))));
        if isempty(diff_max)
            del_num=1;
        else
            del_num=randi(diff_max(1));
        end
        
        e(sec(end-del_num+1:end))=[];
        [Ks0,p] =LS_uneqE(Y,{e});
        power=p{1,1};
        fv= getME_uneqE({e},p,Ks0,Y);
        
        %如果archive没有K对应的解
        locK=find(archiveFV(:,1)==numel(e));
        if isempty(locK)
            archiveE{end+1}=e;
            archiveP{end+1}=power;
            archiveFV(end+1,:)=fv;
        elseif fv(2)<archiveFV(max(lock-1,1),2)
            lock=max(lock-1,1);
            archiveE{lock}=e;
            archiveP{lock}=power;
            archiveFV(lock,:)=fv;
        end
    end 

loc_ar1=find(archiveFV(:,1)==K);
if ~isempty(loc_ar1)
    [~,loc_ar2]=min(archiveFV(loc_ar1,2));
    new_s=archiveE{loc_ar1(loc_ar2)};
    if sum(abs(new_s-old_s))<eta
        tt=tt+1;
    end
    old_s=new_s;
end

G=G+1;    
end
archiveFV=archiveFV(archiveFV(:,1)<=M,:);
archiveE=archiveE(archiveFV(:,1)<=M,:);
%% find knee
archiveFV=[0, norm(Y,'fro'); archiveFV];
knees1=findknees1(archiveFV); lock1=round(knees1);
lock=lock1-1;
theta_est=archiveE{lock};

%% assume K is known
[~,x]=LS_uneqE(Y,{theta_est});
x=x{1,1}';
XEst=x';
power=mean(abs(x),1);
% plot(archiveFV(:,1),archiveFV(:,2),'-ro','MarkerSize',8,'LineWidth',1.5)
% a=archiveE{4};[~,x]=LS_uneqE(Y,{a});
% x=x{1,1}'; [a' diag(x'*x)]
% 
% a=archiveE{5};[~,x]=LS_uneqE(Y,{a});
% x=x{1,1}'; [a' diag(x'*x)]
% 
% a=archiveE{6};[~,x]=LS_uneqE(Y,{a});
% x=x{1,1}'; [a' diag(x'*x)]
% 
% a=archiveE{7};[~,x]=LS_uneqE(Y,{a});
% x=x{1,1}'; [a' diag(x'*x)]
end  