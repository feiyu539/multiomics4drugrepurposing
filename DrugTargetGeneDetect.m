% drug reposition pipeline
% ldy modified for GitHub version
clear;

%% --- Load input data ---
driver_gene_list = csvread(['.' filesep 'input' filesep 'driver_gene_list.csv'], 1, 1);
drug_cluster_detail = csvread(['.' filesep 'input' filesep 'muscle_drug_cluster_detail.csv'], 0, 1);
gene_cluster_detail = csvread(['.' filesep 'input' filesep 'drivernetwork_detail.csv'], 0, 1);

% Placeholder for drug-target matrix, should be loaded as needed
% targetMat = [drugID, geneID] pairs
targetMat = []; % You should load or define targetMat properly

%% --- Detect cluster IDs from directory ---
files = dir(fullfile('.','bfrm'));
clusters = {files(3:end).name}; % skip . and ..
clustersMat = sort(cellfun(@(x) str2double(x(12:end)), clusters));

%% --- Process each drug cluster ---
drug_gene_Matlist = cell(1,length(clustersMat));
for cluster = 1:length(clustersMat)
    fprintf('Processing drugcluster %d ...\n', clustersMat(cluster));

    drugID = drug_cluster_detail(cluster,:);
    drugID = drugID(drugID~=0);

    % Load matrices
    clusterPath = fullfile('.', 'bfrm', ['drugcluster' num2str(clustersMat(cluster))]);
    mA = load(fullfile(clusterPath, 'mA.txt'));
    mF = load(fullfile(clusterPath, 'mF.txt'));
    mPostPib = load(fullfile(clusterPath, 'mPostPib.txt'));

    % Pre-process
    mA(:,1) = [];
    mF(1,:) = [];
    mPostPib(:,1) = [];

    % Binarize mPostPib
    c = median(mPostPib);
    binaryFilterMat = double(mPostPib > c);

    % Weighted interaction matrix
    weightMat = mA .* binaryFilterMat;

    % Calculate etaMat
    etaMat = zeros(size(mF));
    for j = 1:size(etaMat,2)
        for i = 1:size(etaMat,1)
            targetGenes = DrugTargetGeneDetect(targetMat, drugID(j));
            [~, dgi, ~] = intersect(driver_gene_list, targetGenes);
            etaMat(i,j) = sum(weightMat(dgi,j)) * mF(i,j);
        end
    end

    % Compute drug-gene interaction matrix
    drug_gene_Mat = weightMat * etaMat;
    drug_gene_Matlist{cluster} = drug_gene_Mat;
end

%% --- Concatenate drug-gene matrices ---
final_drug_geneMat = [];
for cluster = 1:length(drug_gene_Matlist)
    final_drug_geneMat = [final_drug_geneMat, drug_gene_Matlist{cluster}];
end

%% --- Create etamdMat list by gene clusters ---
etamdMatlist = cell(1, size(gene_cluster_detail,1));
for cluster = 1:size(gene_cluster_detail,1)
    geneID = gene_cluster_detail(cluster,:);
    geneID = geneID(geneID ~= 0);
    [~, dgi, ~] = intersect(driver_gene_list, geneID);
    etamdMatlist{cluster} = final_drug_geneMat(dgi,:);
end

%% --- Compute L1 norm for each cluster across drugs ---
TEMat = zeros(length(etamdMatlist), size(final_drug_geneMat,2));
for i = 1:size(TEMat,2)
    for j = 1:length(etamdMatlist)
        TEMat(j,i) = norm(etamdMatlist{j}(:,i), 1); % L1 norm
    end
end

%% --- Extract drug IDs ---
drugIDlist = [];
for cluster = 1:length(clustersMat)
    drugID = drug_cluster_detail(cluster,:);
    drugID = drugID(drugID ~= 0);
    drugIDlist = [drugIDlist, drugID];
end

%% --- Save results ---
save('processResults.mat', 'TEMat', 'drugIDlist', 'etamdMatlist', 'drug_gene_Matlist', 'driver_gene_list', 'drug_cluster_detail', 'gene_cluster_detail');
writematrix(TEMat', 'TEMat.csv');
writematrix(drugIDlist', 'druglist.csv');

%% --- Function Definition ---
function [targetGenes] = DrugTargetGeneDetect(targetMat, drugID)
    targetGenes = [];
    for i = 1:size(targetMat,1)
        if targetMat(i,1) == drugID
            targetGenes = [targetGenes, targetMat(i,2)];
        end
    end
end
