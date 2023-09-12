function[best_view,sill,sumD]=select_best_view(X,num_clust)
% select the best view
num_views = length(X);
sill = [];
for a = 1: num_views
     S{a} = L2_distance_1(X{a},X{a});
    [idx{a},~,sumD{a}]=kmeans(S{a},num_clust,'maxiter',1000,'replicates',20,'EmptyAction','singleton');
    sill = [sill,mean(silhouette(S{a},idx{a},'cityblock'))];
    sumD{a} = sum(sumD{a});
end
[~,best_view] = max(sill);
end