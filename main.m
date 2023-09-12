clear memory
clear all;
addpath('./datasets');
addpath([pwd, '/funs']);

%% option
filename = 'bbcsport3vbigRnSp.mat';
% Datanames :
% \'bbcsport3vbigRnSp.mat' \ 3Source.mat \'ORL.mat'\ 
repNum=25;

%% data_preprocessing
load(filename);
if strcmp('ORL.mat',filename)
    Y=gt;
     best_view = 1;    
elseif strcmp('3Source.mat',filename)
        for i=1:length(X)
            X{i}=X{i}';
        end
        Y=truth;
        best_view = 3;
elseif strcmp('bbcsport3vbigRnSp.mat',filename)
    Y=truth;
    best_view = 3;
end
view_num=length(X);
t1=clock;

%% data pre processing
for v = 1:view_num
    X{v}=NormalizeData(X{v});
end
num_clust = length(unique(Y));

%% iteration
for iteration=1:repNum
    if mod(iteration-1,5)==0
        [X,~,mapping] = shuffle_data(X,Y,best_view);
    end
    [result,P] = T_UMC(X,Y,filename,best_view,iteration);
    tempACC(iteration) = result(1);
    tempNMI(iteration) = result(2);
    tempPurity(iteration) = result(3);
    tempARI(iteration) = result(4);
    tempFscore(iteration) = result(5);
    tempPrecision(iteration) = result(6);
    tempRecall(iteration) = result(7);
end

ACC = [mean(tempACC),std(tempACC)];
NMI = [mean(tempNMI),std(tempNMI)];
Purity = [mean(tempPurity),std(tempPurity)];
ARI = [mean(tempARI),std(tempARI)];
Fscore = [mean(tempFscore),std(tempFscore)];
Precision = [mean(tempPrecision),std(tempPrecision)];
Recall = [mean(tempRecall),std(tempRecall)];
t2 = clock;
Time = etime(t2,t1)/repNum;

ACC = roundn(ACC,-4)
NMI = roundn(NMI,-4)
Purity = roundn(Purity,-4)
ARI = roundn(ARI,-4)
Fscore =roundn(Fscore,-4)
Precision =roundn(Precision,-4)
Recall =roundn(Recall,-4)
Time = roundn(Time,-4)

% compute re-coupling accuracy
for v = 1: length(X)
    a1{v} = mappingsACC(P{v},mapping{best_view},1);%CPA@1
    a2{v} = mappingsACC(P{v},mapping{best_view},3);%CPA@3
    a3{v} = mappingsACC(P{v},mapping{best_view},10);%CPA@10
end
a1
a2
a3



