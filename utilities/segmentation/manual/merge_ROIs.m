function [Amerge, fluo_merged, merged_ROIs, delta_merged] = merge_ROIs(A, FOVframes_1d, overlap_thr)

nROIs_all = size(A,2);
id_ROIs = 1:1:nROIs_all;
Amerge = zeros(size(A));
merged_ROIs = zeros(size(A));
n_ROIs_merged = 0;
while ~isempty(id_ROIs)
    n_ROIs_merged = n_ROIs_merged + 1;
    id_ROI_temp = id_ROIs(1);
    merged = id_ROI_temp; %array with iD overlapping ROIs
    for j = 1:length(id_ROIs)
        if id_ROI_temp~=id_ROIs(j)
           print1 = find(A(:,id_ROI_temp)>0);
           print2 = find(A(:,id_ROIs(j))>0);
           overlap = intersect(print1,print2);
           if length(overlap) >= min(overlap_thr*length(print1),overlap_thr*length(print2)) %compute overlap
               Amerge(print1,n_ROIs_merged) = max([A(print1,id_ROI_temp) A(print1,id_ROIs(j))],[],2); %update spatial matrix
               Amerge(print2,n_ROIs_merged) = max([A(print2,id_ROI_temp) A(print2,id_ROIs(j))],[],2); %update spatial matrix
               merged = [merged id_ROIs(j)]; %update overlap vector
               A(overlap,id_ROI_temp)=max([A(overlap,id_ROI_temp) A(overlap,id_ROIs(j))],[],2); %update print of reference ROI
           end
        elseif id_ROI_temp==id_ROIs(j)
            print1 = find(A(:,id_ROI_temp)>0);
            Amerge(print1,n_ROIs_merged) = A(print1,id_ROI_temp);
        end
    end
    merged_ROIs(merged,n_ROIs_merged) = 1;
    %remove merged ROIs
    for j = 1:length(merged)
        id_ROIs(id_ROIs==merged(j)) = [];
    end
end
merged_ROIs(:,sum(merged_ROIs,1)==0) = [];
Amerge(:,sum(Amerge,1)==0) = [];
[~,t_FOV] = size(FOVframes_1d);
fluo_merged = zeros(n_ROIs_merged,t_FOV);
for i = 1:n_ROIs_merged
    fluo_merged(i,:) = mean(FOVframes_1d(find(Amerge(:,i)),:),1);
end
delta_merged = size(A,2)-size(Amerge,2);