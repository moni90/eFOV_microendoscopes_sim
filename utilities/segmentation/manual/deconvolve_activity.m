function [C,C_Df] = deconvolve_activity(caiman_path, A, fluo, FOVframes_1d)

addpath(genpath(caiman_path));

C = zeros(size(fluo));
for i = 1:size(A,2)
    [C(i,:),~,~,~,~,~] = constrained_foopsi(fluo(i,:));
end
[C_Df] = extract_DF_F(FOVframes_1d,A,C);