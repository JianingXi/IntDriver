# IntDriver
Discovering potential driver genes through an integrated model of somatic mutation profiles and gene functional information

![image](https://github.com/JianingXi/IntDriver/blob/master/image/splash.jpg)

Developer: Jianing Xi <xjn@mail.ustc.edu.cn> from Health Informatics Lab, School of Information Science and Technology, University of Science and Technology of China

## Instructions to IntDriver (version 1.0.0)

Requirement
------------------------
* 8GB memory
* MATLAB R2015a or later


Gene functional information
------------------------
1. Interactome information
The file `./FunctionalInformation/GeneSymbol_net.mat` is the list of investigated genes, which is also the genes included in interaction network [iRefIndex 9](http://irefindex.org).
The file `./FunctionalInformation/network_adj_matrix.mat` is the adjanceny matrix of interaction network [iRefIndex 9](http://irefindex.org). The entry Aij equaling 1 represents connection between the i-th gene and j-th gene related to the list above, and 0 otherwise.
2. Gene similarity information from Gene Ontology terms
The file `./FunctionalInformation/GO_similarity_matrix.mat` is the GO-based similarity matrix of the investigated genes, which is calculated by a previous published software CS2 [https://academic.oup.com/bioinformatics/article/25/9/1178/204335/GS2-an-efficiently-computable-measure-of-GO-based](https://academic.oup.com/bioinformatics/article/25/9/1178/204335/GS2-an-efficiently-computable-measure-of-GO-based).

Input data descriptions
------------------------
We provide the input data files of TCGA somatic mutation data of three cancers from [cBioPortal](http://www.cbioportal.org/data_sets.jsp) respectively.

        =================================================================================================
        | CANCER NAME                             |SAMPLE SIZE | FILE DIRECTORY                         |
        =================================================================================================
        |Breast invasive carcinoma (BRCA)         | 507        |`./Input/BRCA.mat`                      |
        -------------------------------------------------------------------------------------------------
        |Kidney renal clear cell carcinoma (KIRC) | 424        |`./Input/KIRC.mat`                      |
        -------------------------------------------------------------------------------------------------
        |Lung squamous cell carcinoma (LUSC)      | 183        |`./Input/LUSC.mat`                      |
        -------------------------------------------------------------------------------------------------

Output data descriptions
------------------------
The descriptions of output variables of mCGfinder are provided below:

        =================================================================================================
        | VARIABLE NAME        | DESCRIPTION                                                            |
        =================================================================================================
        |GeneSelected          |The top 200 genes selected by IntDriver as potential driver genes.      |
        -------------------------------------------------------------------------------------------------
        |U_sample_indicator    |The sample indicator matrix indicates the assignment of tumor samples to|
        |                      |each layer of the reconstructed matrix. The entry U_ij = 1 indicate that|
        |                      |the i-th samples are assigned to the j-th layer, and 0 otherwise.       |
        -------------------------------------------------------------------------------------------------
        |V_gene_score          |The gene score matrix for each layers of the reconstructed matrix and   |
        |                      |the investigated genes. Higher scores represent the larger possibility  |
        |                      |of the related genes to be potential driver genes.                      |
        -------------------------------------------------------------------------------------------------

run IntDriver
------------------------
To apply IntDriver on the example input datasets with the default configurations, please run the script file `./demo.m` and the result file will be automatically saved as `.mat` file in directory `./Output` when the program is finished.

Parameter configurations
------------------------
The configurations of IntDriver can be set as a Struct variable `Conf` in `./demo.m`, with their descriptions provided below:

        =================================================================================================
        | PARAMETER NAME       | DESCRIPTION                                                            |
        =================================================================================================
        |Conf.lambda_N         |The tuning parameter of interactome information for integrating network |
        |                      |information into the predicted mutation scores of genes.                |
        -------------------------------------------------------------------------------------------------
        |Conf.lambda_S         |The tuning parameter of gene similarity information, which is used to   |
        |                      |incorporate gene simliarity (like Gene Ontology based similarity) into  |
        |                      |the mutation scores of genes.                                           |
        -------------------------------------------------------------------------------------------------
        |Conf.lambda_R         |The tuning parameter of Frobenius norm based regularization to prevent  |
        |                      |over-fitting problem. The default value is 0.01.                        |
        -------------------------------------------------------------------------------------------------
        |Conf.Max_K            |Maximum rank/ number of layers of the reconstructed matrix. The default |
        |                      |value is 5.                                                             |
        -------------------------------------------------------------------------------------------------
        |Conf.LeastProp        |Least proportion of samples assigned to each layers. The default value  |
        |                      |is 15%.                                                                 |
        -------------------------------------------------------------------------------------------------


