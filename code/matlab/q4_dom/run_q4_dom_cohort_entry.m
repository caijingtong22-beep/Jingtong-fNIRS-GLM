function out = run_q4_dom_cohort_entry()
addpath(fileparts(fileparts(mfilename('fullpath')))); out = run_q1_procedure_cohort_estimation("dom","q4");
end
