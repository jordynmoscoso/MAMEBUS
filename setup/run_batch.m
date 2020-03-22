function run_batch(foldername,date,ws_batch,sd_batch)
addpath ~/Desktop/MAMEBUS/setup/ % path to the setup file
%%% 
%%%
%%% run_batch
%%%
%%%
%%% creates a parameter sweep based on the input parameters

if ws_batch && sd_batch
    disp('Cannot currently run two parameter sweeps at once')
    return
end

% reference variables
ws_ref = 0.05;      % reference windstress for all sd_batch tests
sd_ref = 75;        % reference shelf depth for all ws_batch tests

% parameter sweep vectors
% notes: 
%   for wind-stress it seems most reasonalbe up until 0.0125, I haven't
%   done anything more than that, can extend the param sweeps more than
%   that, I think though.
ws_vec = round(logspace(log10(0.02),log10(0.2),10),3);
sd_vec = [200, 175, 150, 125, 100, 75, 60, 55, 50, 45];

% internal model parameters
t1day = 86400;
t1year = t1day*365;
ecosystem_model_type = 2; % npzd model
outputFreq = 10*t1day;
endTime = 25*t1year;

% cluster information: each folder has a run file, so we just need to
% create a *.sh file that uploads and downloads everything.
use_cluster = true;

% batch information
local_home_dir = '~/Desktop';
folderpath = fullfile(local_home_dir,foldername);
mkdir(folderpath)
run_batch_name = 'run_batch.sh';
run_batch_file = fopen(fullfile(folderpath,run_batch_name),'w');
upload_batch_file = fopen(fullfile(folderpath,'upload_batch.sh'),'w');
download_batch_file = fopen(fullfile(folderpath,'download_batch.sh'),'w');

% create folders for the ws_batch and an upload/download script
if ws_batch
    disp(['Creating a batch with ' num2str(length(ws_vec)) ' files'])
    wildcard = ['ws_batch_',num2str(date),'_tau_*'];
    for nWS = 1:length(ws_vec)
        % create run names
        run_name = ['ws_batch_',num2str(date),'_tau_',num2str(ws_vec(nWS))];
        setparams(folderpath,run_name,ecosystem_model_type,outputFreq,endTime,-ws_vec(nWS),sd_ref,use_cluster);
        
        % write commands to the necessary run/upload/download files
        fprintf(run_batch_file,'cd %s\n',run_name);
        fprintf(run_batch_file,'sh Build_MAMEBUS.sh\n');
        fprintf(run_batch_file,'sh run.sh\n');
        fprintf(run_batch_file,'cd ..\n');
    end
end

if sd_batch
    disp(['Creating a batch with ' num2str(length(sd_vec)) ' files'])
    wildcard=['sd_batch_',num2str(date),'_depth_*'];
    for nSD = 1:length(sd_vec)
        % create run names
        run_name = ['sd_batch_',num2str(date),'_depth_',num2str(sd_vec(nSD))];
        setparams(folderpath,run_name,ecosystem_model_type,outputFreq,endTime,-ws_ref,sd_vec(nSD),use_cluster);
        
        % write commands to the necessary run/upload/download files
        fprintf(run_batch_file,'cd %s\n',run_name);
        fprintf(run_batch_file,'sh Build_MAMEBUS.sh\n');
        fprintf(run_batch_file,'sh run.sh\n');
        fprintf(run_batch_file,'cd ..\n');
    end
end


fprintf(upload_batch_file,'rsync -av --update %s %s\n',['../' foldername],'jmoscoso@caolila.atmos.ucla.edu:/data3/jmoscoso/MAMEBUS/runs');
fprintf(download_batch_file,'rsync -av --update %s %s\n',['jmoscoso@caolila.atmos.ucla.edu:/data3/jmoscoso/MAMEBUS/runs/',foldername,'/',wildcard],'./');

end