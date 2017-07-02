function [U_vec,V_vec,DetectLayerNum] = IntDriver(X_input,NetInfo,SimInfo,Conf)

Default_lambda_N = 0.3;
if ~isfield(Conf,'lambda_N') || Conf.lambda_N <= 0
    Conf.lambda_N = Default_lambda_N;
end
Default_lambda_S = 0.7;
if ~isfield(Conf,'lambda_S') || Conf.lambda_S <= 0
    Conf.lambda_S = Default_lambda_S;
end
Default_lambda_R = 0.01;
if ~isfield(Conf,'lambda_R') || Conf.lambda_R <= 0
    Conf.lambda_R = Default_lambda_R;
end 
DefaultMax_K = 5;
if ~isfield(Conf,'Max_K') || Max_K <= 0
    Conf.Max_K = DefaultMax_K;
end
DefaultLeastProp = 0.15;
if ~isfield(Conf,'LeastProp') || Conf.LeastProp <= 0 || Conf.LeastProp >= 1 
    Conf.LeastProp = DefaultLeastProp;
end
[Num_Sample, Num_Gene] = size(X_input);
LeastNumLayer = max(ceil(Conf.LeastProp*Num_Sample),10);

U_vec = zeros(Num_Sample,DefaultMax_K);
V_vec = zeros(Num_Gene,DefaultMax_K);

LapMat = Conf.lambda_S*SimInfo + Conf.lambda_N*NetInfo;
clear SimInfo NetInfo

selected_idx_sample = zeros(Num_Sample,1);
NumSampleRemain = Num_Sample;

DetectLayerNum = 0;
while sum(selected_idx_sample) < Num_Sample && DetectLayerNum < (Conf.Max_K-1) && NumSampleRemain >= LeastNumLayer
    
    [U_vec_1,V_vec_cur] = GetRankOneLayer(X_input,Conf,LapMat,LeastNumLayer);
    DetectLayerNum = DetectLayerNum + 1;
    cur_idx_s = (U_vec_1~=0);
    U_vec_cur = zeros(Num_Sample,1);
    U_vec_cur(selected_idx_sample==0) = U_vec_1;
    selected_idx_sample(selected_idx_sample==0) = cur_idx_s;
    U_vec(:,DetectLayerNum) = U_vec_cur;
    V_vec(:,DetectLayerNum)  = V_vec_cur;
    X_input = X_input(~cur_idx_s,:);
    
    NumSampleRemain = sum(~selected_idx_sample);
    if NumSampleRemain == 0
        return;
    end
end
% --- If remain samples exist---
if NumSampleRemain > 0
    NumSampleRemain = size(X_input,1);
    U_vec_1 = ones(NumSampleRemain,1);
    V_vec_cur = V_vec_Est(X_input,U_vec_1,Conf,LapMat);
    DetectLayerNum = DetectLayerNum + 1;
    U_vec_cur = zeros(Num_Sample,1);
    U_vec_cur(selected_idx_sample==0) = U_vec_1;
    U_vec(:,DetectLayerNum) = U_vec_cur;
    V_vec(:,DetectLayerNum)  = V_vec_cur;
end

ind_sample = (sum(abs(U_vec),1)>0);
ind_gene = (sum(abs(V_vec),1)>0);
U_vec = U_vec(:,ind_sample);
V_vec = V_vec(:,ind_gene);

end


function [U_vec,V_vec,iter] = GetRankOneLayer(X_input,Conf,LapMat,leastGroupNum)

[n_samp,p_gene] = size(X_input);
if sum(sum(abs(X_input))) == 0
    U_vec = ones(n_samp,1)/n_samp;
    V_vec = zeros(p_gene,1);
    iter = 0;
    return;
end

try
    [U_vec_0,~,~] = svds(X_input,1);
    if sum(U_vec_0<-(10^-4)) > sum(U_vec_0>(10^-4))
        % ensure most element of u0 is positive
        U_vec_0 = -U_vec_0;
    end
    U_vec_0(U_vec_0<=0)=0;
    U_vec_0 = U_vec_0/max(U_vec_0);
    V_vec_0 = X_input'*U_vec_0;
catch
    U_vec_0 = ones(n_samp,1);
    V_vec_0 = sum(X_input,1)';
end

SampleNumInclude = 0;
Relate_merr = 0.005;
eps_temp = 10^-5;
niter = 50;

Relate_U_vec_d = 10;
Relate_V_vec_d = 10;
iter = 0;
while (Relate_U_vec_d > Relate_merr || Relate_V_vec_d > Relate_merr ...
        || SampleNumInclude < leastGroupNum) && iter < niter
    iter = iter+1;
    % ---- V_vec updating ----
    V_vec_1 = V_vec_Est(X_input,U_vec_0,Conf,LapMat);
    % ---- U_vec updating ----
    U_vec_1 = U_vec_binary(X_input,V_vec_1,leastGroupNum);
    SampleNumInclude = sum(U_vec_1~=0);
    % ---- Residual ----
	Relate_U_vec_d = sum((U_vec_0-U_vec_1).^2)/(sum(U_vec_0.^2)+eps_temp);
    Relate_V_vec_d = sum((V_vec_0-V_vec_1).^2)/(sum(V_vec_0.^2)+eps_temp);
    U_vec_0 = U_vec_1;
	V_vec_0 = V_vec_1;
end
U_vec = U_vec_1;
V_vec = V_vec_1;
end

function V_vec_hat = V_vec_Est(X_input,U_vec_fix,Conf,LapMat)
GeneLen = size(X_input,2);
lambda_lap = Conf.lambda_N + Conf.lambda_S;
A_mat = (sum(U_vec_fix.^2)+Conf.lambda_R)*speye(GeneLen) ...
    + (lambda_lap*LapMat);
b_vec = X_input'*U_vec_fix;

V_vec_hat = A_mat\b_vec;
end


function U_vec_hat = U_vec_binary(X_input,V_vec_fix,leastGroupNum)
vec_x = X_input*V_vec_fix/sum(V_vec_fix.^2);
th_candi = sort(vec_x,'descend');
TH_vec = min(th_candi(leastGroupNum),0.5);
U_vec_hat = (vec_x >= TH_vec);
end

