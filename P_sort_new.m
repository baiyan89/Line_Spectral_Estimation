function [Next,non_index,FrontNo,CrowdDis] = P_sort_new(~,NP,FV,selection)
% FV: NP*2, where NP is pupulation size, 2 is the number of objectives.
% Next: 被选择的solution index，共选择NP个
% non_index: 非支配解的index
% FrontNo：population中每个解的front rank
% CrowdDis：crowding distance of all solutions

[N,M] = size(FV);
selection=1;
switch selection
    case 1  % non_dominated sorting and crowding distance
        MaxFno = 0;
        Sorted = false(1,N);
        [FV,rank] = sortrows(FV);
        FrontNo = zeros(1,N) + inf;
        while sum(Sorted)<NP
            MaxFno = MaxFno + 1;
            ThisFront = false(1,N);
            for i = 1 : N
                if ~Sorted(i)
                    x = 0;
                    for j = 1 : N
                        if ThisFront(j)
                            x = 2;
                            for j2 = 2 : M
                                if FV(i,j2) < FV(j,j2)
                                    x = 0;
                                    break;
                                end
                            end
                            if x == 2
                                break;
                            end
                        end
                    end
                    if x ~= 2
                        ThisFront(i) = true;
                        Sorted(i) = true;
                    end
                end
            end
            FrontNo(rank(ThisFront)) = MaxFno;
        end
        
        Next = FrontNo < MaxFno;
        CrowdDis = CrowdingDistance(FV,FrontNo);
        Last     = find(FrontNo==MaxFno);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:NP-sum(Next)))) = true;
%         Next=find(Next~=0);
        non_index=find(FrontNo==1);
        FrontNo    = FrontNo(Next);
        CrowdDis   = CrowdDis(Next);
        
    case 2 %按k一个个选
        Next=[];a=FV;Fros=cell(1,M-1); kk=1;
        for jj=1:M-1
            locs=find(a(:,1)==jj);
            [~,sec]=sort(a(locs,2),'ascend');
            Fros{jj}=locs(sec);
        end
        while numel(Next)<NP
            for jj=1:M-1
                Next=[Next,Fros{jj}(kk)];
                
                if numel(Next)==NP
                    break;
                end
            end
            kk=kk+1;
        end
        non_index=Next(1:M-1);
end
end