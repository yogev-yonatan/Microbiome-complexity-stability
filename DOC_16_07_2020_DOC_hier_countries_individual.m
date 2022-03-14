%% Loading the data
otu_table = readtable('spongeemp.sub5tk.txt');

%% Loading metadata
load metadata.mat;
otu_table_temp = otu_table;

% Creating struct
sponge.data = otu_table_temp{:,2:end-1};
sponge.OTU = [string(otu_table_temp{:,1}), string(otu_table_temp{:,end})];
sponge.datad = double(sponge.data>0);
sponge.sample_name = string(metadata.sample_name);

%% Filtering sponge tissue
ind = string(metadata.sample_type)=="sponge tissue";
% ind = string(metadata.sample_type)=="seawater";
sponge.data = sponge.data(:,ind);
sponge.datad = sponge.datad(:,ind);
sponge.sample_name = sponge.sample_name(ind);
metadata = metadata(ind,:);
metadata.sampling_site = string(metadata.sampling_site);

%% Finding number of clusters per country
countries = unique(string(metadata.country));
marker_size_1 = 10;
marker_size_2 = 50;
figure;
edges_1 = linspace(0, 1, 100);
edges_2 = linspace(0, 2000, 100);
min_num_samples = 20;
min_num_samples_per_clust = 20;
retry = true;
k = 0;
clust = struct;
country_clust = strings(length(countries), 2);
num_replicats = 10;
% num_clust_vec = [10 0 0 0 2 2 0 2 4 2 1 0 3 2 2 0 2 3 0 2 1 0 0 10 2 2 2 1 0 6 1];
for i = 1:length(countries)
    retry = true;
    if sum(string(metadata.country) == countries(i))<min_num_samples
        continue
    end
    while retry 
        disp(countries(i));
        % Creating data
        ind = metadata.country == countries(i);
        data = sponge.data(:, ind);
        datad = sponge.datad(:, ind);
        metadata_temp = metadata(ind, :);
        
        % Filtering out OTUs not present
        ind = sum(datad, 2) ~= 0;
        data = data(ind, :);
        datad = datad(ind, :);
        OTU = sponge.OTU(ind, :);
        %     num_species = sum(datad);
        
        % All samples
        % tSNE
        subplot(2,2,1);
        y_tsne = tsne(datad', 'Distance', 'jaccard');
        scatter(y_tsne(:, 1), y_tsne(:, 2), marker_size_1, sum(datad), 'filled');
        colorbar;
        axis square;
        title('t-SNE');
        
        % PCoA
        subplot(2,2,2);
        D = pdist(datad', 'jaccard');
        y_pca = cmdscale(D);
        scatter(y_pca(:, 1), y_pca(:, 2), marker_size_1, sum(datad), 'filled');
        colorbar;
        axis square;
        title('PCoA');
        
        % Distances
        subplot(2,2,3);
        histogram(D, edges_1);
        title('Jaccard distances');
        
        % Number of species
        subplot(2,2,4);
        histogram(sum(datad), edges_2);
        title('Number of species');
        
        % Choose number of clusters
        num_clust = input('Number of clusters: ');
%         min_clust = input('Min number of clusters: ');
%         max_clust = input('Max number of clusters: ');
%         num_clust = num_clust_vec(i);
        
        % Clustering Hierarchical
%         Z = linkage(D, 'average');
%         C = cluster(Z, 'maxclust', num_clust);
%         C = cluster(Z, 'cutoff', clust_cutoff);
%         C = cluster(Z, 'criterion', 'distance', 'cutoff', J_cutoff);
%         stemp.C = C;
%         I = inconsistent(Z);
        
%         Clustering kmediods
        C = kmedoids(datad', num_clust, 'Replicates', num_replicats, ...
            'Distance', 'jaccard');

        % Silhoutte
%         sil = zeros(max_clust - min_clust, 1);
%         C_cell = cell(length(sil));
%         mediods_cell = cell(length(sil));
%         num_clust_vec = min_clust:max_clust;
%         k_temp = 0;
%         if min_clust ~= 1
%             for j = min_clust:max_clust
%                 j
%                 k_temp = k_temp + 1;
%                 [C, mediods] = kmedoids(datad', j, 'Replicates', num_replicats, ...
%                     'Distance', 'jaccard');
%                 C_cell{k_temp} = C;
%                 mediods_cell{k_temp} = mediods;
% %                 sil(k_temp) = silhouette(datad', C, 'Jaccard');
%                 temp = evalclusters(datad', C, 'silhouette');
%                 sil(k_temp) = temp.CriterionValues;
%             end
%             plot(min_clust:max_clust ,sil);
%             pause;
%             [~, best_clust] = max(sil);
%             C = C_cell{best_clust};
%             mediods = mediods_cell{best_clust};
%             num_clust = num_clust_vec(best_clust);
%         else
%             [C, mediods] = kmedoids(datad', 1);
%         end
        
        % Plotting clusters
        
        for j = 1:num_clust
            if sum(C==j)>=min_num_samples_per_clust
                % All t-SNE
                subplot(3,2,1);
                hold off;
                gscatter(y_tsne(:, 1), y_tsne(:, 2), C, [], [], marker_size_1, false);
                hold on;
                scatter(y_tsne(C==j, 1), y_tsne(C==j, 2), marker_size_2, [0 0 0], 'X');
                axis square;
                title('All tSNE');
                
                % All PCoA
                subplot(3,2,2);
                hold off;
                gscatter(y_pca(:, 1), y_pca(:, 2), C, [], [], marker_size_1, false);
                hold on;
                scatter(y_pca(C==j, 1), y_pca(C==j, 2), marker_size_2, [0 0 0], 'X');
                axis square;
                title('All PCoA');
                
                % Cluster t-SNE
                subplot(3,2,3);
                yc_tsne = tsne(datad(:, C==j)', 'Distance', 'jaccard');
                scatter(yc_tsne(:, 1), yc_tsne(:, 2), marker_size_1, sum(datad(:, C==j)), 'filled');
                axis square;
                title('Cluster t-SNE');
                
                % Clutser PCoA
                subplot(3,2,4);
                D = pdist(datad(:, C==j)', 'jaccard');
                yc_pca = cmdscale(D);
                scatter(yc_pca(:, 1), yc_pca(:, 2), marker_size_1, sum(datad(:, C==j)), 'filled');
                axis square;
                title('Cluster PCoA');
                
                % Cluster Distances
                subplot(3,2,5);
                histogram(D, edges_1);
                title('Cluster Distances');
                
                % Cluster Number of species
                subplot(3,2,6);
                histogram(sum(datad(:, C==j)), edges_2);
                title('Cluster number of species');
%                 pause;
            end
        end
        
%         Retry?
%         retry = input('Retry: ');
%         retry = logical(retry);
        retry = false;
        
        if not(retry)
            country_clust(i, :) = [countries(i), string(num_clust)];
            for j = 1:num_clust
                if sum(C==j)>=min_num_samples_per_clust
                    k = k + 1;
                    clust(k).data = data(:, C==j);
                    clust(k).datad = datad(:, C==j);
                    clust(k).metadata = metadata_temp(C==j, :);
                    clust(k).country = countries(i);
                    clust(k).OTU = OTU;
%                     clust(k).mediod = mediods(j);
                end
            end
        end
    end
end

%% Removing outliers and DOC
outlier_const = 2;
percent_min = 0.5;
percent_max = 0.9;
figure;
for i = 1:length(clust)
    i
    
    % Remove empty species
    data = clust(i).data;
    datad = clust(i).datad;
    OTU = clust(i).OTU;
    metadata_temp = clust(i).metadata;
    
    ind = sum(datad, 2) ~= 0;
    data = data(ind, :);
    datad = datad(ind, :);
    OTU = OTU(ind, :);
    
    % Remove outliers by jaccard
    distances = pdist2(datad', datad', 'jaccard');
    distances = mean(distances, 1);
    [~, ind] = min(distances);
    datad_mediod = datad(:, ind);
    distances = pdist2(datad_mediod', datad', 'jaccard');
    ind = distances <= mean(distances).*outlier_const;
    
    data = data(:, ind);
    datad = datad(:, ind);
    metadata_temp = metadata_temp(ind, :);
    
    % Normalizing
    data = data./sum(data);
    
    % Calculating DOC
    [diss, overlap] = DOC(data);
    [~, index] = sort(overlap, 'ascend');
    overlap = overlap(index);
    diss = diss(index);
    
    % Calculating slope
    overlap_filt = overlap(round(end*percent_min):round(end*percent_max));
    diss_filt = diss(round(end*percent_min):round(end*percent_max));
    clust(i).fit = fit(overlap_filt, diss_filt, 'poly1');
    
    % Parameters
    overlap_mean = mean(overlap); % What about percent_min/max?
    N_median = median(sum(datad));
    N_mean = mean(sum(datad));
    N_var = var(sum(datad));
    N_std = std(sum(datad));
    
    temp = confint(clust(i).fit);
    conf_int = temp(:, 1);
    temp = coeffvalues(clust(i).fit);
    slope = temp(1);
    
    % Saving values
    clust(i).data = data;
    clust(i).datad = datad;
    clust(i).OTU = OTU;
    clust(i).metadata = metadata_temp;
    clust(i).N = sum(datad);
    clust(i).N_median = N_median;
    clust(i).N_mean = N_mean;
    clust(i).N_var = N_var;
    clust(i).N_std = N_std;
    clust(i).conf_int = conf_int;
    clust(i).overlap_mean = overlap_mean;
    clust(i).slope = slope;
    
    % DOC plot
%     plot(overlap, diss, '.', 'MarkerSize', 20);
%     pause;
end

%% Plotting DOC
figure; hold on;
% N_median = [clust.N_median];
N_mean = [clust.N_mean];
N_median = [clust.N_median];
N_std = [clust.N_std];
slope = [clust.slope];
overlap_mean = [clust.overlap_mean];
conf_int = zeros(length(clust), 2);
country = categorical([clust.country]);
for i = 1:length(clust)
    conf_int(i, 1) = clust(i).conf_int(1);
    conf_int(i, 2) = clust(i).conf_int(2);
end

min_overlap = 0.6;
ind = overlap_mean > min_overlap;

N_median = N_median(ind);
N_mean = N_mean(ind);
N_std = N_std(ind);
slope = slope(ind);
conf_int = conf_int(ind, :);
overlap_mean = overlap_mean(ind);
country = country(ind);


% slope = slope(randperm(length(slope)));

subplot(1,2,1);
hold on;
load d.mat;
% h2 = scatter(av_Species_n, D_squared, 100, [0 0 0], 's', 'filled');
% h = scatter(N_median, slope, 100, overlap_mean, 'filled');
h = scatter(N_median, slope, 100, [1 0 0], 'filled');
% h = gscatter(N_median, slope, country, [], [], 100);
eb(1) = errorbar(N_median, slope, conf_int(:, 1), conf_int(:, 2), 'vertical', 'LineStyle', 'none','Color','k','LineWidth', 0.1);
eb(2) = errorbar(N_median, slope, N_std,'horizontal','LineStyle','none','Color','k','LineWidth', 0.1);
colormap('jet');
xlabel('N');
ylabel('Slope');
set(gca, 'FontSize', 20);
title('Real');

subplot(1,2,2);
hold on;
slope = slope(randperm(length(slope)));
h = scatter(N_median, slope, 100, overlap_mean, 'filled');
eb(1) = errorbar(N_median, slope, conf_int(:, 1), conf_int(:, 2), 'vertical', 'LineStyle', 'none','Color','k','LineWidth', 0.1);
eb(2) = errorbar(N_median, slope, N_std,'horizontal','LineStyle','none','Color','k','LineWidth', 0.1);
colormap('jet');
xlabel('N');
ylabel('Slope');
set(gca, 'FontSize', 20);
title('Rand');
