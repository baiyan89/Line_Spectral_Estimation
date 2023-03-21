function Offspring_E=mutation(Offspring_E,P)
NP=size(Offspring_E,1);

for jj=1:NP
    e=Offspring_E{jj};
    D=numel(e);
%     s=rand;
%     if s<0.33 && length(e)>1
%         mutation_strategy=2;
%         
%     elseif s>0.66 && length(e)<10-1
%         mutation_strategy=3;
%         
%     else
%         mutation_strategy=1;        
%     end   
    mutation_strategy=1;
    
    switch mutation_strategy  %在无噪声的情况下，发现mutation_strategy取1要比随机取1,2,3的结果要好。
        case 1
            % mutation parameters
            ProM = 1/D;  DisM =15; 
            
            % 定义boundaries
            MaxValue=90*ones(1,D);
            MinValue=-90*ones(1,D);
            
            % polynomial mutation
            k    = rand(1,D);
            miu  = rand(1,D);
            Temp = k<=ProM & miu<0.5;
            e(Temp) = e(Temp)+(MaxValue(Temp)-MinValue(Temp)).*((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(e(Temp)-MinValue(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
            Temp = k<=ProM & miu>=0.5;
            e(Temp) = e(Temp)+(MaxValue(Temp)-MinValue(Temp)).*(1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-e(Temp))./(MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
           
            % sort
            Offspring_E{jj}=sort(e,'ascend');
    
        case 2 % randomly delete one elemeent
            loc=randi([1, length(e)]);
            Offspring_E{jj}(loc)=[];           
            
        case 3 % randomly add one element
            e=[e, 180*rand-90];
            Offspring_E{jj}=sort(e,'ascend');           
    end
end