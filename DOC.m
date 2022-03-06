function [overlap_vector,dissimilarity_vector] = DOC(X,dissimilarity_mode)
%Calculate the DOC of the input data X
%Yogev Yonatan, yogev.yn@gmail.com
% X- input data, rows of X correspond to species,
% columns of X correspond to observations.
% 
% dissimilarit_mode - The dissimilarity metric:
% 'e' for Euclidean distance
% 'rjsd' for root Jensen–Shannon divergence
% 'spearman' for Spearman distance
%
% overlap_vector and dissimilarity_vecotr are the output overlap and
% dissimilarity values between each pair of observations, respectievly.

O = @(x,y) 0.5*(sum(x)+sum(y));%Overlap function


normalized_results = X./sum(X,1);

exams = size(X,2);
overlap1 = zeros(exams);
dissimilarity = zeros(exams);

k = 0;
overlap_vector = zeros(exams*(exams-1)/2,1);
dissimilarity_vector = zeros(exams*(exams-1)/2,1);
switch dissimilarity_mode
    case 'rjsd'
        KLD = @(x,y) sum(x(x>0).*log(x(x>0)./y(x>0)));
        rJSD=@(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2));
        
        for i = 1:exams
            for j = 1+i:exams
                k = k +1;
                Shared_species1 = find(normalized_results(:,i)>0 & normalized_results(:,j)>0);
                overlap1(i,j) = O(normalized_results(Shared_species1,i),normalized_results(Shared_species1,j));
                xhat = normalized_results(Shared_species1,i)./sum(normalized_results(Shared_species1,i));
                yhat =normalized_results(Shared_species1,j)./sum(normalized_results(Shared_species1,j));
                dissimilarity(i,j) = rJSD(xhat,yhat);
                overlap_vector(k) = overlap1(i,j);
                dissimilarity_vector(k) = dissimilarity(i,j);
            end
        end
    case 'e'
        
        for i = 1:exams
            for j = 1+i:exams
                k = k +1;
                Shared_species1 = find(normalized_results(:,i)>0 & normalized_results(:,j)>0);
                overlap1(i,j) = O(normalized_results(Shared_species1,i),normalized_results(Shared_species1,j));
                xhat = normalized_results(Shared_species1,i)./sum(normalized_results(Shared_species1,i));
                yhat =normalized_results(Shared_species1,j)./sum(normalized_results(Shared_species1,j));
                dissimilarity(i,j) = pdist2(xhat',yhat','euclidean');
                overlap_vector(k) = overlap1(i,j);
                dissimilarity_vector(k) = dissimilarity(i,j);
            end
        end
    case 'spearman'
        for i = 1:exams
            for j = 1+i:exams
                k = k +1;
                Shared_species1 = find(normalized_results(:,i)>0 & normalized_results(:,j)>0);
                overlap1(i,j) = O(normalized_results(Shared_species1,i),normalized_results(Shared_species1,j));
                xhat = normalized_results(Shared_species1,i)./sum(normalized_results(Shared_species1,i));
                yhat =normalized_results(Shared_species1,j)./sum(normalized_results(Shared_species1,j));
                switch isempty(xhat)
                    case true
                        dissimilarity(i,j)=0;
                    case false
                        dissimilarity(i,j) = 1-corr(xhat,yhat,'type','Spearman');                        
                end
                overlap_vector(k) = overlap1(i,j);
                dissimilarity_vector(k) = dissimilarity(i,j);
            end
        end
        
    otherwise
        err('incorrect dissimilarity mode ')
end

[overlap_vector,idxs]= sort(overlap_vector);
dissimilarity_vector=dissimilarity_vector(idxs);

