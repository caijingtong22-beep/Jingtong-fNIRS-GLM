function run_q2_procedure_dom_cohort()
% Run the shared Q1--Q4 baseline GLM for Q2 Procedure and DOM.
% The shared implementation selects HILE from the question key.
run_q1_procedure_cohort_estimation("procedure", "q2");
run_q1_procedure_cohort_estimation("dom", "q2");
end
