function [audit,status] = run_q9_procedure_cohort_entry(subject_numbers)
addpath(fileparts(fileparts(mfilename('fullpath'))));
if nargin<1, [audit,status] = run_q9_procedure_cohort(); else, [audit,status] = run_q9_procedure_cohort(subject_numbers); end
end
