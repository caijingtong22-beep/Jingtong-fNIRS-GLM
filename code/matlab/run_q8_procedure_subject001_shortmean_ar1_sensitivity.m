function [contrast_results, diagnostics] = run_q8_procedure_subject001_shortmean_ar1_sensitivity(channel_id)
% Sensitivity model for Q8 Procedure subject001:
% [LILEProcedure, HIHEProcedure, intercept, cubic drift, short mean],
% estimated channel-wise with feasible AR(1) GLS.  This is not a replacement
% for the NIRS-SPM model; it tests the impact of short-channel nuisance data.

if nargin < 1, channel_id = 1; end
subject_label = 'subject001';
this_file = mfilename('fullpath'); repo_root = fileparts(fileparts(fileparts(this_file)));
glm_root = 'C:\\aliceproject\\GLM'; first_level = fullfile(glm_root, 'first_level');
q8_root = fullfile(repo_root, 'results', 'q8_hihe_vs_lile_procedure', 'cohort_specification', subject_label);
input_file = fullfile(first_level, 'nirs_spm_inputs_v2_fulltimeline_trial1', subject_label, 'subject001_nirs_data.mat');
metadata_file = fullfile(first_level, 'nirs_spm_inputs', subject_label, 'subject001_channel_metadata.csv');
output_root = fullfile(q8_root, 'shortmean_ar1_sensitivity'); if ~isfolder(output_root), mkdir(output_root); end
assert(isfile(input_file) && isfile(metadata_file), 'Missing subject001 input or channel metadata.');

data = load(input_file, 'nirs_data'); nirs_data = data.nirs_data; clear data;
metadata = readtable(metadata_file, 'TextType', 'string');
assert(height(metadata) == nirs_data.nch, 'Channel metadata mismatch.');
short_idx = find(string(metadata.channel_type) == "Short");
assert(~isempty(short_idx), 'No short channels were found.');
t = double(nirs_data.time(:)); n = numel(t); assert(channel_id >= 1 && channel_id <= nirs_data.nch, 'Invalid channel_id.');
u = linspace(-1, 1, n)';

contrast_results = table(); diagnostics = table();
for hb_label = ["HbO", "HbR"]
    spm_file = fullfile(q8_root, char(hb_label), sprintf('SPM_indiv_%s.mat', hb_label));
    S = load(spm_file, 'SPM_nirs'); SPM_nirs = S.SPM_nirs;
    assert(size(SPM_nirs.xX.X,1) == n, 'SPM/data timepoints mismatch.');
    names = string(SPM_nirs.xX.name(:)); task_idx = ~contains(lower(names), 'constant');
    assert(sum(task_idx) == 2, 'Expected exactly two task columns.');
    task_X = full(SPM_nirs.xX.X(:,task_idx));
    if hb_label == "HbO", Y = double(nirs_data.oxyData); else, Y = double(nirs_data.dxyData); end
    short_mean = mean(Y(:,short_idx), 2, 'omitnan'); short_z = zscore_safe(short_mean);
    X = [task_X, ones(n,1), u, u.^2, u.^3, short_z];
    regressor_names = ["LILEProcedure"; "HIHEProcedure"; "Intercept"; "LinearDrift"; "QuadraticDrift"; "CubicDrift"; "ShortChannelMean_z"];
    assert(rank(X) == size(X,2), '%s AR1 sensitivity design is rank deficient.', hb_label);
    writetable(array2table([t X], 'VariableNames', ['timestamp_sec'; regressor_names]), ...
        fullfile(output_root, sprintf('%s_%s_shortmean_ar1_design.csv', subject_label, hb_label)));
    c = [ -1; +1; zeros(size(X,2)-2,1) ];
    result_rows = table(); diag_rows = table(); beta_rows = table();
    for ch = 1:nirs_data.nch
        y = Y(:,ch); valid = isfinite(y) & all(isfinite(X),2); yv = y(valid); Xv = X(valid,:);
        b_ols = Xv \ yv; rho = estimate_ar1_rho(yv - Xv*b_ols, .95);
        y_ar = yv(2:end) - rho*yv(1:end-1); X_ar = Xv(2:end,:) - rho*Xv(1:end-1,:);
        b = X_ar \ y_ar; resid_ar = y_ar - X_ar*b; dof = numel(y_ar) - rank(X_ar);
        mse = sum(resid_ar.^2) / dof; cov_b = mse * pinv(X_ar'*X_ar); se = sqrt(max(diag(cov_b),0));
        cb = c'*b; cse = sqrt(max(c'*cov_b*c,0)); ct = cb/cse; cp = two_sided_t_p(ct,dof);
        fitted = Xv*b; residual = yv-fitted; r2 = 1-sum(residual.^2)/sum((yv-mean(yv)).^2);
        result_rows = [result_rows; table(hb_label,ch,string(metadata.channel_id(ch)),string(metadata.channel_type(ch)),string(metadata.roi(ch)), ...
            cb,cse,ct,cp,rho,dof,r2, 'VariableNames', {'hemoglobin','channel_index','channel_id','channel_type','roi','contrast_beta','contrast_se_ar1_gls','contrast_t_ar1_gls','contrast_p_ar1_gls','ar1_rho','residual_dof','r_squared_unwhitened'})]; %#ok<AGROW>
        raw_residual_ac1 = estimate_ar1_rho(residual,.999);
        whitened_residual_ac1 = estimate_ar1_rho(resid_ar,.999);
        diag_rows = [diag_rows; table(hb_label,ch,rho,raw_residual_ac1,whitened_residual_ac1,r2,dof, 'VariableNames', {'hemoglobin','channel_index','estimated_ar1_rho','raw_residual_ac1','whitened_residual_ac1','r_squared_unwhitened','residual_dof'})]; %#ok<AGROW>
        for k=1:numel(b)
            beta_rows = [beta_rows; table(hb_label,ch,regressor_names(k),b(k),se(k), 'VariableNames', {'hemoglobin','channel_index','regressor','beta','se_ar1_gls'})]; %#ok<AGROW>
        end
        if ch == channel_id
            fitted_full = nan(n,1); residual_full = nan(n,1); fitted_full(valid) = fitted; residual_full(valid) = residual;
            fig = figure('Visible','off','Color','w','Position',[100 100 1550 850]);
            set(fig,'DefaultAxesColor','w','DefaultAxesXColor','k','DefaultAxesYColor','k','DefaultTextColor','k'); tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
            nexttile; plot(t,y,'Color',[0 .35 .7],'LineWidth',.7); hold on; plot(t,fitted_full,'Color',[.85 .33 .1],'LineWidth',1.2); grid on; xlabel('Time (s)'); ylabel('Concentration'); title(sprintf('%s channel %d %s: short-mean + AR(1) GLS fit',subject_label,ch,hb_label)); legend({'Observed','Fitted'},'Location','eastoutside');
            nexttile; plot(t,residual_full,'Color',[0 .35 .7],'LineWidth',.7); hold on; yline(0,'--k'); grid on; xlabel('Time (s)'); ylabel('Residual'); title(sprintf('Raw residual; R^2 = %.3f, raw AC(1) = %.3f, whitened AC(1) = %.3f',r2,raw_residual_ac1,whitened_residual_ac1));
            nexttile; plot(t,[task_X short_z],'LineWidth',1); grid on; xlabel('Time (s)'); ylabel('Regressor'); title('HRF-convolved Q8 tasks and short-channel mean nuisance'); legend({'LILEProcedure','HIHEProcedure','ShortChannelMean z'},'Location','eastoutside');
            exportgraphics(fig,fullfile(output_root,sprintf('%s_channel%02d_%s_shortmean_ar1_fit.png',subject_label,ch,hb_label)),'Resolution',220); close(fig);
        end
    end
    result_rows.contrast_q_bh61 = bh_fdr(result_rows.contrast_p_ar1_gls); contrast_results = [contrast_results;result_rows]; diagnostics = [diagnostics;diag_rows]; %#ok<AGROW>
    writetable(result_rows,fullfile(output_root,sprintf('%s_%s_q8_HIHE_minus_LILE_shortmean_ar1_channel_contrasts.csv',subject_label,hb_label)));
    writetable(beta_rows,fullfile(output_root,sprintf('%s_%s_shortmean_ar1_betas.csv',subject_label,hb_label)));
end
writetable(diagnostics,fullfile(output_root,sprintf('%s_shortmean_ar1_diagnostics.csv',subject_label)));
end

function z=zscore_safe(x), z=(x-mean(x,'omitnan'))./std(x,0,'omitnan'); end
function rho=estimate_ar1_rho(r,limit), r=r(isfinite(r)); rho=sum(r(2:end).*r(1:end-1))/max(sum(r(1:end-1).^2),eps); rho=max(-limit,min(limit,rho)); end
function p=two_sided_t_p(t,df), p=betainc(df/(df+t^2),df/2,1/2); end
function q=bh_fdr(p), [s,o]=sort(p); m=numel(p); qs=s.*m./(1:m)'; qs=min(1,flipud(cummin(flipud(qs)))); q=zeros(size(p)); q(o)=qs; end
