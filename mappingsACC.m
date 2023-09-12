function [ acc ] = mappingsACC(P,map,k)
[A,C] = sort(P,2,'descend');
[tops] = C(:,1:k);
cnt = 0;
for i = 1:size(P,1)   
    for j = 1:k
       if(map(1,i) == tops(i,j))
          cnt = cnt+1; 
       end
    end
end
acc = cnt/size(P,1);
end

