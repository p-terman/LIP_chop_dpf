function LUXTrackXGridJobs(jobs, flag_restart_failed_jobs, flag_percent)
% LUXTrackXGridJobs(jobs, flag_restart_failed_jobs, flag_percent)
% 
% Uses jobs output structure from LUX01BuildRQ1s to track xgrid jobs.
% Break out of monitoring loop with liberal use of Ctrl-C.
% 
% Inputs:
%  jobs: output structure from LUX01BuildRQ1s
%  flag_restart_failed_jobs: if you want LUXTrackXGridJobs to restart
%   failed jobs, set this flag to 1; 0 otherwise [default: 1]
%   flag_percent: 1 to find percent done and estimated time remaining in job output. [default: 1] 
%
% 090519 DCM v1.0 -- Haven't fully tested action on failure yet. This needs
%  to be done
% 20090924 JJC added ability to read percent done and estimated completion time for each job
% from output of LUX01REEFrunner_v3. This is done by default, pass it a zero if you want to abstain.

if nargin < 2
    flag_restart_failed_jobs = 1;
end

if nargin<3
	flag_percent=1;
end

num_jobs_to_complete = length(jobs);
num_completed_jobs = 0;

dis('\n\n*** Now following progress. ***\n\n');

% Put in status for all jobs
for ii_job = 1:length(jobs)
    jobs(ii_job).status = 0;
end

while ~isempty(jobs)
    dis('Checking remaining jobs (%d of %d complete):', num_completed_jobs, num_jobs_to_complete);
    for ii_job = 1:length(jobs)
        [status, attributes] = system(dis('xgrid -auth Kerberos -job attributes -id %1.0f', jobs(ii_job).id));
        pause(1);
        
        if flag_percent
        [status, results] = system(dis('xgrid -auth Kerberos -job results -id %1.0f', jobs(ii_job).id));
        pause(1);
        job_is_string = 'Job is ';
        job_is_string_ind = strfind(results,job_is_string);
        percent_complete_string = '% complete';
        percent_complete_string_ind = strfind(results,percent_complete_string);
        if ~isempty(percent_complete_string_ind)
        	percent_done_start = job_is_string_ind(end)+length(job_is_string);
        	percent_done_end = percent_complete_string_ind(end)-1;
        	percent_done_string = results(percent_done_start:percent_done_end);
        	percent_done = str2num(percent_done_string);
        end
        
        est_time_string = '(estimated time remaining is ';
        est_time_string_ind = strfind(results,est_time_string);
        seconds_string = ' seconds)';
        seconds_string_ind = strfind(results,seconds_string);
        if ~isempty(seconds_string_ind)
        	est_time_start = est_time_string_ind(end) + length(est_time_string);
        	est_time_end = seconds_string_ind(end)-1;
        	est_time_remaining_seconds = str2num(results(est_time_start:est_time_end));
        	est_time_remaining_string = secs2hms(est_time_remaining_seconds);
        end
   		end     
        % parse xgrid attributes output
        job_status_indicator = 'jobStatus = ';
        k = strfind(attributes, job_status_indicator);
        job_status_str = strtok(attributes(k+length(job_status_indicator):end), ' ;');
        %keyboard;
        
        hostname_string_start = '<hostname>';
        hostname_string_end = '</hostname>';
       	hostname_start_ind = strfind(results,hostname_string_start);
        hostname_end_ind = strfind(results,hostname_string_end);
        if ~isempty(hostname_start_ind)
        	hostname_start = hostname_start_ind(1) + length(hostname_string_start);
        	hostname_end = hostname_end_ind(1)-2;
        	hostname = results(hostname_start:hostname_end);
        else
        	hostname = 'gsk-??.het.brown.edu';
        end
        
        switch job_status_str
            case 'Pending'
                dis('xgrid job %d pending on %s', jobs(ii_job).id, hostname);
            case 'Running'
                fprintf('xgrid job %d running on %s', jobs(ii_job).id, hostname);
                if flag_percent && ~isempty(percent_complete_string_ind) && ~isempty(seconds_string_ind)
                	fprintf(' - %2.1f%% complete',percent_done);
           			fprintf('  estimated time remaining is %s\n', est_time_remaining_string);
                else
                	fprintf('\n');
                end
            case 'Finished' % check to make sure all files have been created and are up-to-date
                for ii_file = 1:length(jobs(ii_job).file_list)
                    % look for file in output directory
                    file = dir(sprintf('%s/%s_f%09.0f.rq1',jobs(ii_job).output_path, ...
                        jobs(ii_job).filename_prefix, jobs(ii_job).file_list(ii_file)));
                    if isempty(file) || (jobs(ii_job).options.flag_overwrite && jobs(ii_job).timestamp > datenum(file.date))
                        dis('xgrid job %d failed on %s', jobs(ii_job).id, hostname);
                        jobs(ii_job).status = -1;
                        break;
                    end
                end
                if jobs(ii_job).status == 0
                    dis('xgrid job %d finished successfully on %s', jobs(ii_job).id, hostname);
                    jobs(ii_job).status = 1;
                end
            otherwise
                dis('job %d: invalid job identifier -- no longer monitoring', jobs(ii_job).id);
                jobs(ii_job).status = 1;
        end
    end

    ii_job = 1;
    while ii_job <= length(jobs)
        if flag_restart_failed_jobs && jobs(ii_job).status == -1 % job has failed -- resubmit
            %options.flag_overwrite=1;
            %options.flag_xgrid=1;
            jobs(ii_job).options.num_jobs = 1;
            dis('\n\nResubmitting xgrid job %d', jobs(ii_job).id);
            [status, new_jobs] = LUXBuildRQ1s(jobs(ii_job).filename_prefix, ...
                                                jobs(ii_job).file_list, ...
                                                jobs(ii_job).data_path, ...
                                                jobs(ii_job).output_path, ...
                                                jobs(ii_job).analysis_xml_file, ...
                                                jobs(ii_job).options);
            for ii_new_job = 1:length(new_jobs) % incorporate new jobs into current jobs
                %jobs(end+1) = new_jobs(ii_new_job);
                jobs(end+1).filename_prefix = new_jobs(ii_new_job).filename_prefix;
                jobs(end).file_list = new_jobs(ii_new_job).file_list;
                jobs(end).data_path = new_jobs(ii_new_job).data_path;
                jobs(end).output_path = new_jobs(ii_new_job).output_path;
                jobs(end).options = new_jobs(ii_new_job).options;
                jobs(end).status = 0;
                jobs(end).id = new_jobs(ii_new_job).id;
                jobs(end).timestamp = new_jobs(ii_new_job).timestamp;
            end
            jobs(ii_job) = []; % erase old job
        elseif jobs(ii_job).status == 1 % job has succeeded
            jobs(ii_job) = []; % erase old job
            num_completed_jobs = num_completed_jobs+1;
        else
            ii_job=ii_job+1;
        end
    end
    
    dis('\n');
end
