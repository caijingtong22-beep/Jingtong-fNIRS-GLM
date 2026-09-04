function run_q3_procedure_dom_cohort()
% Run the shared Q1--Q4 baseline GLM for Q3 Procedure and DOM.
% The shared implementation selects LIHE from the question key.
run_q1_procedure_cohort_estimation("procedure", "q3");
run_q1_procedure_cohort_estimation("dom", "q3");
end
