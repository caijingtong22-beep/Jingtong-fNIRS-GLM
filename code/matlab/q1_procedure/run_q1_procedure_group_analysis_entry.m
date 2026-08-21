function out = run_q1_procedure_group_analysis_entry()
addpath(fileparts(fileparts(mfilename('fullpath')))); out = summarize_q1_procedure_explicit_baseline_group("procedure","q1");
end
