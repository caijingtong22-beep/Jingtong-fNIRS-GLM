function [summary,errors]=run_q8_dom_cohort_estimation(subject_numbers,phase_key)
% Estimate the QA-approved Q8 DOM specifications.
if nargin<1,subject_numbers=1;end;if nargin<2,phase_key="dom";end;phase_key=lower(string(phase_key));assert(ismember(phase_key,["dom","aoe"]));subject_numbers=subject_numbers(:)';assert(all(ismember(subject_numbers,1:30)));
repo_root=fileparts(fileparts(fileparts(mfilename('fullpath'))));glm_root='C:\aliceproject\GLM';first_level=fullfile(glm_root,'first_level');root=fullfile(repo_root,'results',sprintf('q8_hihe_vs_lile_%s',phase_key),'cohort_specification');assert(isfolder(root),'Missing Q8 %s specification directory.',upper(phase_key));
addpath(first_level);addpath(fullfile(glm_root,'toolboxes','spm'));addpath(genpath(fullfile(glm_root,'toolboxes','NIRS_SPM_v4_r1')));summary=table();errors=table();
for sn=subject_numbers
    sid=sprintf('subject%03d',sn);data_file=q8_dom_data_file(sn,sid,first_level);D=load(data_file,'nirs_data');nch=D.nirs_data.nch;clear D
    for hb_label=["HbO","HbR"]
        spm_file=fullfile(root,sid,char(hb_label),sprintf('SPM_indiv_%s.mat',hb_label));timer=tic;status="success";message="";skipped=false;nr=NaN;nc=NaN;df=NaN;
        try
            assert(isfile(spm_file),'Missing specification: %s',spm_file);L=load(spm_file,'SPM_nirs');SPM_nirs=L.SPM_nirs;
            if isfield(SPM_nirs.nirs,'step')&&strcmpi(string(SPM_nirs.nirs.step),'estimation'),skipped=true;else,nirs_spm_estimation_batch_fixed(spm_file,data_file);L=load(spm_file,'SPM_nirs');SPM_nirs=L.SPM_nirs;end
            assert(strcmpi(string(SPM_nirs.nirs.step),'estimation'),'%s did not reach estimation.',hb_label);nr=size(SPM_nirs.nirs.beta,1);nc=size(SPM_nirs.nirs.beta,2);df=double(SPM_nirs.xX.erdf);assert(isequal([nr nc],[3 nch]),'%s beta dimensions are %d x %d; expected 3 x %d.',hb_label,nr,nc,nch);assert(all(isfinite(SPM_nirs.nirs.beta),'all'),'Non-finite beta.');fprintf('%s %s: beta %dx%d, df %.0f\n',sid,hb_label,nr,nc,df);
        catch ME
            status="failed";message=string(getReport(ME,'extended','hyperlinks','off'));errors=[errors;table(string(sid),sn,hb_label,string(ME.identifier),message,'VariableNames',{'subject_id','subject_num','hemoglobin','error_identifier','error_message'})]; %#ok<AGROW>
        end
        summary=[summary;table(string(sid),sn,hb_label,status,skipped,nr,nc,df,toc(timer),'VariableNames',{'subject_id','subject_num','hemoglobin','status','skipped_existing','n_betas','n_channels','effective_df','elapsed_sec'})]; %#ok<AGROW>
    end
end
tag=sprintf('%03d_to_%03d',min(subject_numbers),max(subject_numbers));writetable(summary,fullfile(root,sprintf('q8_%s_estimation_%s.csv',phase_key,tag)));writetable(errors,fullfile(root,sprintf('q8_%s_estimation_errors_%s.csv',phase_key,tag)));assert(isempty(errors),'Q8 %s estimation had %d failures.',upper(phase_key),height(errors));
end
function f=q8_dom_data_file(sn,sid,first_level)
if ismember(sn,[1 2 10]),f=fullfile(first_level,'nirs_spm_inputs_v2_fulltimeline_trial1',sid,sprintf('%s_nirs_data.mat',sid));else,f=fullfile(first_level,'nirs_spm_inputs',sid,sprintf('%s_nirs_data.mat',sid));end;assert(isfile(f),'Missing NIRS input: %s',f);
end
