function index = TournamentSelection(K,N,varargin)
%TournamentSelection - Tournament selection.
%
%   P = TournamentSelection(K,N,fitness1,fitness2,...) returns the indices
%   of N solutions by K-tournament selection based on their fitness values.
%   In each selection, the candidate having the MINIMUM fitness1 value will
%   be selected; if more than one candidates have the same minimum value of
%   fitness1, then compare their fitness2 values, and so on.
%
%   Example:
%       P = TournamentSelection(2,100,FrontNo)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%     varargin = cellfun(@(S)reshape(S,[],1),varargin,'UniformOutput',false);
%     [~,rank] = sortrows([varargin{:}]);
%     [~,rank] = sort(rank);
%     Parents  = randi(length(varargin{1}),K,N);
%     [~,best] = min(rank(Parents),[],1);
%     index    = Parents(best+(0:N-1)*K);   
    
   %% 保证每一对父代中，至少有一个父代是非支配解
    non_index0=find(varargin{1}==1);
    if length(non_index0)<N/2
        non_index=[non_index0, non_index0(randi(length(non_index0),1,N/2-length(non_index0)))];
    else
%         eval('error');
        non_index=non_index0(1:N/2);
    end   
    mating_strategy=1;
    switch mating_strategy
        case 1 %非支配解与随机选的某个解配对
            all_index=randi(N,1,N/2);
            
        case 2 %一个非支配解和一个支配解配对
            all_index=setdiff(1:N/2,non_index0);
            all_index=all_index(randi(numel(all_index),1,N/2));
            
    end
    index=nan(1,N);
    index(1:2:N)=non_index;
    index(2:2:N)=all_index;
end