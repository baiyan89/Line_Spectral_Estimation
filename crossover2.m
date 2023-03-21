function Offspring_E=crossover2(E,P,M)
% execute crossover
% E: DOAs of solutions;
% Ks: sparsity of solutions

Offspring_E=cell(size(E));
for i=1 : size(E,1)/2
    g1=E{i}';
    power1=sum(abs(P{i}),2);
    
    g2_ind=i+size(E,1)/2;
    g2=E{g2_ind}';        
    
    while length(g1)==1 && length(g2)==1
        g2_ind=randi([1, size(E,1)]);
        g2=E{g2_ind}';
    end    
    power2=sum(abs(P{g2_ind}),2);
    
%     g1=[-26.4474  -24.2987    5.8604   55.6261   64.5779]';
%     g2=[-88.6248  -69.9628  -65.3046  -23.9567   17.1366   28.8912   62.3002]';
%     power1=[36.3001   41.4124   14.5992   26.4648   23.9506]';
%     power2=[410.1390  784.1889  454.0848    9.2945   20.2779   13.8903   82.5250]';
    
    g=[g1;g2];
    [gsort,secg]=sort(g,'ascend'); power=[power1;power2]; power=power(secg);
    [~,secg]=sort(secg,'ascend');

    %% get distance between all metavariable pairs in g1 and g2
    d=abs(g1 - g2');
    
    % Get most similar counterpart to each metavariable
    [~,n1] = min(d, [], 2);%n1: g1离第几个g2最近
    [~,n2] = min(d, [], 1);%n2: g2离第几个g1最近    
     
    
    % find mutual links
    sec1=[1:1:length(n1)]';
    sec2=[1:1:length(n2)]';   
    s1=[sec1, n1];%看g1->g2的单向link
    s2=[n2', sec2];%看g1<-g2的单向link    
    [coms,locs1,locs2]=intersect(s1,s2,'rows');
    
    S1=gsort; S2=gsort; pS1=power; pS2=power;
    S1(secg(length(g1)+1:end))=nan;  pS1(secg(length(g1)+1:end))=nan;
    S2(secg(1:length(g1)))=nan;  pS2(secg(1:length(g1)))=nan;
    
    [~, locS2]=ismember(g2(locs2),S2);
    [~, locS1]=ismember(g1(locs1),S1);
    
    S1(locS2)=[];  pS1(locS2)=[];
    S2(locS1)=[];  pS2(locS1)=[];
    
    % crossover
%     com_num=numel(unique(n1));%mutual-link的个数
%     crosP_num=randi([1, max(length(S1)-1, 1)]);%cross points的个数
    crosP_num=randi([1, min([length(g1), length(g2), 1])]);%cross points的个数
    sec=randperm(length(S1)-1);
    crosL=sort(sec(1:crosP_num)); crosL=[0,crosL,length(S1)];
    
    
    child1 = []; child2 = [];   pChild1=[];pChild2=[];
    for j=1:crosP_num+1
        if crosL(j)+1<=crosL(j+1)
            if mod(j,2)==1
                child1=[child1; S1(crosL(j)+1:crosL(j+1))];
                child2=[child2; S2(crosL(j)+1:crosL(j+1))];
                pChild1=[pChild1; pS1(crosL(j)+1:crosL(j+1))];
                pChild2=[pChild2; pS2(crosL(j)+1:crosL(j+1))];
                
            else
                child1=[child1; S2(crosL(j)+1:crosL(j+1))];
                child2=[child2; S1(crosL(j)+1:crosL(j+1))]; 
                pChild1=[pChild1; pS2(crosL(j)+1:crosL(j+1))];
                pChild2=[pChild2; pS1(crosL(j)+1:crosL(j+1))];                
            end
        end
    end
    child1(isnan(child1))=[];
    child2(isnan(child2))=[];
    pChild1(isnan(pChild1))=[];
    pChild2(isnan(pChild2))=[];
    
    %% 满足最大source number<=M
    %去掉的element个数随机，优先去掉power较小的element
%     if length(child1)>=M
%         [~,sec]=sort(pChild1,'ascend');
%         del_num=max(randi(length(child1)-(M-1)),randi(numel(child1)) );
%         child1(sec(1:del_num))=[];
%     end
% 
%     if length(child2)>=M
%         [~,sec]=sort(pChild2,'ascend');
%         del_num=max(randi(length(child2)-(M-1)), randi(numel(child2)));
%         child2(sec(1:del_num))=[];
%     end
 
   
    
    %% 保存子代解到Offspring_E
    Offspring_E{i}=child1';
    Offspring_E{i+size(E,1)/2}=child2';
end
Offspring_E(cellfun(@isempty,Offspring_E))=[];
end

% function c=ainb(a,b)
% if isempty(a)==1&&isempty(b)==1
% ? ? c=[];
% ? ? return;
% end
% if length(a)>length(b)
% ? ? d=a;
% ? ? a=b;
% ? ? b=d;
% end
% c=nan(length(a),1);%%提前预设
% l=1;
% for i=1:length(a)
% ? ? if any(b==a(i))
% ? ? ? c(l)=a(i);
% ? ? ? l=l+1;
% ? ? ?end
% end
% c(isnan(c)==1,)=[];
% end
