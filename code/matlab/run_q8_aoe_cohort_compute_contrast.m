function [summary,errors]=run_q8_aoe_cohort_compute_contrast(subject_numbers)
% AOE wrapper for the shared Q8 two-regressor contrast implementation.
if nargin<1,subject_numbers=1;end
[summary,errors]=run_q8_dom_cohort_compute_contrast(subject_numbers,"aoe");
end
