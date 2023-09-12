function [datas,Yi,mapping] = shuffle_data( raw_data,Y,best_view )
% x is dv*n
n_views = length(raw_data); 
for i =1 :n_views
   raw_data{i}= raw_data{i}'; 
end
datas = {};
for a = 1 : n_views 
    n_view_a = size(raw_data{a}, 1);
    d_view_a = size(raw_data{a}, 2);
    mapping{a} = randperm(n_view_a);  % random mapping
    if a == best_view
        mapping{a} =(1:n_view_a); 
    end
    datas{a} = zeros(size(raw_data{a}));
    for i = 1:n_view_a
       datas{a}(i,:) = raw_data{a}(mapping{a}(i),:);
    end
end
Yi = zeros(n_view_a,1);
    for i = 1:n_view_a
        Yi(i,:)  = Y(mapping{best_view}(i),:);
    end
for i =1 :n_views
   datas{i}= datas{i}'; 
end
end

