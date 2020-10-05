function [C,C_Df] = deconvolve_activity(A, fluo, FOVframes_1d)

addpath(genpath('/home/calcium/Monica/endoscopes_project/code/simulations_PSF_fit/Matlab code_jb_copy/ETIC/'));

C = zeros(size(fluo));
for i = 1:size(A,2)
    [C(i,:),~,~,~,~,~] = constrained_foopsi(fluo(i,:));
end
[C_Df] = extract_DF_F(FOVframes_1d,A,C);