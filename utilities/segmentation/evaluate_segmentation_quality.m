function [Recall, Precision, F1] = evaluate_segmentation_quality(A_groundtruth, A_segment)

bin_thr = 1e-3;

A_gt = 1*(A_groundtruth > bin_thr);
A_auto = 1*(A_segment > bin_thr);

Thresh=0.25;
[Recall, Precision, F1] = GetPerformance_Jaccard(A_gt,A_auto,Thresh);
