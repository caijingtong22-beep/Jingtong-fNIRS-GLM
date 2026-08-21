function [audit,status] = run_q9_aoe_cohort_entry(subject_numbers)
addpath(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1, [audit,status] = run_q9_dom_cohort("aoe"); else, [audit,status] = run_q9_dom_cohort("aoe",subject_numbers); end
end
