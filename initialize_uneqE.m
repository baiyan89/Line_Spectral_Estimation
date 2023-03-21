function [E, E2] = initialize_uneqE(Y,NP)
% intialize a set of DOA combinations.
[M,~]=size(Y);
E = cell(NP,1); 
E2 = cell(NP,1); 

integers=randi([1, M-1], 1, NP);
%% intialization 1: ±Èinitializaiton 2ÒªºÃ
for i=1:NP     
    E{i}=180*rand(1,integers(i))-90;
    E{i}=sort(E{i},2);
    E2{i}=-E{i};
end

% E{1}=[ -83.76 -6.69 12.35  ];
%% intialization 2
% sec=-90:1:90;
% for i=1:NP     
%     E{i}=sec(randi([1 numel(sec)],1,integers(i)));
%     E{i}=min(max(E{i}+rand(1,integers(i))-0.5, -90),90);
%     E{i}=sort(E{i},2);
%     E2{i}=-E{i};
% end
end
