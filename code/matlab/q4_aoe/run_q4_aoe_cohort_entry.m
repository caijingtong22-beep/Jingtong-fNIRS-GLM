function out = run_q4_aoe_cohort_entry()
addpath(fileparts(fileparts(mfilename('fullpath')))); out = run_q1_procedure_cohort_estimation("aoe","q4");
end
