function[My_result, P]=T_UMC(X,gt,filename,best_view,epoch)
%% Initialization
cls_num = length(unique(gt));
Label = double(gt);
num_views = length(X); 
N = size(X{1},2); % view number; sample number
for k = 1:num_views
    pre_P{k} = zeros(N,N);
    P{k} = zeros(N,N);
    Z{k} = zeros(N,N);
    E{k} = zeros(size(X{k},1),N);
    % auxiliary parameter
    M{k} = eye(N,N);
    K{k} = zeros(N,N);
    % multiipliers
    Y{k} = zeros(N,N);
    J{k} = zeros(size(X{k},1),N);
    G{k} = zeros(N,N);
    w{k} = 1/ num_views;
end
M{best_view} = eye(N,N);
P{best_view} = eye(N,N);
I = eye(N,N);
sX = [N, N, num_views];
iter = 0;Isconverg = 0; num = 0;
epson = 1e-7; mu = 10e-5; max_mu = 10e10; pho_mu = 2;
% construct laplacian matrix
for v = 1:num_views
    fea_v = X{v};   
    Wl = L2_distance_(fea_v); 
    d = sum(Wl); 
    Dl = zeros(N);
    for i = 1:N
        Dl(i,i) = d(i);
    end
    temp_ = diag(sqrt( diag(Dl).^(-1) ));
    L{v} = eye(N)-temp_*Wl*temp_; 
end 

%% setting parameters
if strcmp('ORL.mat',filename)
    alpha =0.5;
    beta = 0.3;
    gamma = 0.5;
elseif strcmp('3Source.mat',filename) 
    alpha = 0.3;
    beta = 0.8;
    gamma = 0.3;
elseif strcmp('bbcsport3vbigRnSp.mat',filename) 
    alpha = 1;
    beta = 1;
    gamma = 0.5;
else
    error('filename not defined!')
end

%% iteration 
while(Isconverg ==0)
    iter = iter+1;
    fprintf('---------processing epoch %d ----- iter %d--------\n', epoch,iter);
    num = num + 1;
    for k = 1: num_views 
        %% Update Z^k
        tmp_A = (mu*X{k}'*X{k}) ;
        tmp_B = (2*alpha*L{k}+2*gamma*M{k}*M{k}'+mu*P{k}*P{k}');
        tmp_C = (2*gamma*Z{best_view}*M{k}' + mu*X{k}'*X{k} - mu*X{k}'*E{k} + X{k}'*J{k} + mu*K{k}*P{k}' - Y{k}*P{k}');
        Z{k} = lyap(tmp_A,tmp_B,-tmp_C);
        clear tmp_A tmp_B tmp_C
        
        %% Update M^v
        M{k} = pinv(2*gamma*Z{k}'*Z{k}+mu*I)*(2*gamma*Z{k}'*Z{best_view}+mu*P{k}+G{k});

        %% Update P^v
        P{k} = pinv(Z{k}'*Z{k}+I)*(M{k}-G{k}./mu+Z{k}'*K{k}-Z{k}'*Y{k}./mu);
        
        %% Update E^k
        for vi = 1 : size(X,2)
                if vi == 1
                F = [X{1}-X{1}*Z{1}+J{1}/mu];
                else
                F =[F; X{vi}-X{vi}*Z{vi}+J{vi}/mu;];
                end
        end
        [Econcat] = solve_l1l2(F,beta/mu);
        for ei = 1:num_views
            if ei == 1
                concatlength = 0;
            end
             E{ei} = Econcat((concatlength+1):(concatlength+size(X{ei},1)),:);  
             concatlength = concatlength + size(X{ei});
        end    
        tmp_PZ{k} = Z{k}*P{k};   
        %% Update Jv Gv
        J{k} = J{k} + mu * (X{k}-X{k}*Z{k}-E{k});
        G{k} = G{k} + mu * (P{k}-M{k});
    end
    %% Update K
    Z_tensor = cat(3, tmp_PZ{:,:});
    clear tmp_PZ
    Y_tensor = cat(3,Y{:,:});
    z = Z_tensor(:);
    y = Y_tensor(:);
    [k, ~] = wshrinkObj_weight(z + 1/mu*y,1/mu,sX,0,3);
    K_tensor = reshape(k,sX);
    %% Update Y
    y = y+mu*(z-k); 
    %% Update mu
    mu = min(mu*pho_mu, max_mu);
    %% Converge Condition
    Isconverg = 1;
    max_Z = 0;
    max_P_M = 0;
    for k = 1: num_views
        if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
            history.norm_Z= norm(X{k}-X{k}*Z{k}-E{k},inf);
            fprintf('norm_Z %7.10f\n', history.norm_Z);
            Isconverg = 0;
            max_Z=max(max_Z,history.norm_Z );
        end
        K{k}=K_tensor(:,:,k);
        Y_tensor = reshape(y,sX);
        Y{k}=Y_tensor(:,:,k);
        if (norm(P{k}-M{k},inf)>epson)
            history.norm_P_M= norm(P{k}-M{k},inf);
            fprintf('norm_P_M %7.10f\n', history.norm_P_M);
            Isconverg = 0;
            max_P_M=max(max_P_M,history.norm_P_M );
            pre_P{k}=P{k};
        end
    if (iter > 300)
        Isconverg = 1;
        error('---------NOT CONVERGENT!--------\n');
    end
    end
end 
tmp=zeros(N,N);
for v=1: num_views
        tmp = tmp + Z{v}*P{v};
end
Clus = SpectralClustering(1/(2*num_views)*(abs(tmp)+abs(tmp')), cls_num );
[result(num,:)] = ClusteringMeasure1(Label, Clus);
My_result = [result(num,:)];   


    