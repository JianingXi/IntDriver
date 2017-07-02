% --- configure ---
Conf.lambda_N = 0.3;
Conf.lambda_S = 0.7;
Conf.lambda_R = 0.01;

load('./FunctionalInformation/GO_similarity_matrix.mat');
load('./FunctionalInformation/network_adj_matrix.mat');
load('./FunctionalInformation/GeneSymbol_net.mat');

Input_data = load('./Input/BRCA.mat');

[U_sample_indicator,V_gene_score] = ...
    IntDriver(Input_data.mutation_mat,network_adj_matrix,Go_sim_matrix,Conf);
[~,ind_gene] = sort(max(V_gene_score,[],2),'descend');
GeneSelected = GeneSymbol_net(ind_gene(1:200));

mkdir('./Output')
save('./Output/Result_BRCA.mat','U_sample_indicator','V_gene_score','GeneSelected');
clear network_adj_matrix Go_sim_matrix Conf Input_data ind_gene
