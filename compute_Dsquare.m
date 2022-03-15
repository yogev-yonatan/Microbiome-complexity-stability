function D2 = compute_Dsquare(overlap,dissimilarity,p)

% This function computes D^2 from DOC
% Yogev Yonatan, yogev.yn@gmail.com
% 
% Inputs:
% overlap - vector of overlap values
% dissimilarity - vector of dissimilarity values
% p - the fraction of data points for D^2 computation (based on the highers
% overlap values' data)

% Output:
% D2 - the D^2 value

[overlap,i]=sort(overlap,'ascend');
dissimilarity=dissimilarity(i);


l = length(overlap);
idx = round(l*(1-p)):l;

f = fit(overlap(idx),dissimilarity(idx),'poly1');
D2 = (f.p1)^2;

