% Pitch detection algorithms evaluation for speech
% 
% See the Bridge project website at Wireless Communication and Networking
% Group (WCNG), University of Rochester for more information:
% http://www.ece.rochester.edu/projects/wcng/project_bridge.html
%
% Written by He Ba, Na Yang, and Weiyang Cai, University of Rochester.
% Copyright (c) 2013 University of Rochester.
% Version February 2013.
%
% Permission to use, copy, modify, and distribute this software without
% fee is hereby granted FOR RESEARCH PURPOSES only, provided that this
% copyright notice appears in all copies and in all supporting
% documentation, and that the software is not redistributed for any fee.
%
% For any other uses of this software, in original or modified form,
% including but not limited to consulting, production or distribution
% in whole or in part, specific prior permission must be obtained from WCNG.
% Algorithms implemented by this software may be claimed by patents owned
% by WCNG, University of Rochester.
%
% The WCNG makes no representations about the suitability of this
% software for any purpose. It is provided "as is" without express
% or implied warranty.
%

close all; 
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Parameter settings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_per_window = 0.06; % frame length in second
f0max = 600; % higher bound of human speech pitch
f0min = 50; % lower bound of human speech pitch
timestep = 0.01; % time offset of detected pitch (second)
error_tolerance = 0.1; % error tolerance for pitch detection on speech is 10%
pthr1 = 1.1; % amplitude ratio of the first two highest Cepstrum peaks, used for voiced/unvoiced detection for Cepstrum pitch detection.
SNR_dB1 = [0 3 7 10 15 20]; % SNR levels
noise_type1 = {'babble','destroyerengine','destroyerops','factory','hfchannel','pink','volvo','white'};
groundtruth_path = 'groundtruth/speech/gt2_';
noisyfile_path = 'Generated files/Noisy speech/Noisy Arctic/';

% save workspace with timestamp
date_now = clock;
date_now = strcat(num2str(date_now(1)), num2str(date_now(2)), num2str(date_now(3)), '-', num2str(date_now(4)), num2str(date_now(5)));
saveFilename = ['Generated files/Saved workspace/Arctic-' date_now '.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Choose speech database for pitch detection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arctic database
list ={'arctic_bdl_a0003(m)',...
    'arctic_bdl_a0012(m)',...
    'arctic_clb_a0005(f)',...
    'arctic_clb_b0525(f)',...
    'arctic_clb_a0225(f)',...
    'arctic_rms_a0080(m)',...
    'arctic_rms_b0286(m)',...
    'arctic_rms_b0514(m)',...
    'arctic_slt_a0590(f)',...
    'arctic_slt_b0045(f)'};

% % Harvard Sentences database
% list = {'OSR_us_000_0015_8k_9(f)',...
% 'OSR_us_000_0017_8k_9(f)',...
% 'OSR_us_000_0017_8k_10(f)',...
% 'OSR_us_000_0057_8k_2(m)',...
% 's1_edited(f)',...
% 's2_edited(m)',...
% 's3_edited(f)',...
% 's4_edited(m)', ...
% 's5_edited(m)',...
% 's6_edited(m)'};

% % LDC database
% list ={'cc_001(m)_happy_14',...
% 'cc_001(m)_interest_10',...
% 'jg_001(f)_sadness_passive_negative_9',...
% 'mf_001(m)_boredom_6',...
% 'mf_001(m)_elation_5',...
% 'mf_001(m)_hotAnger_14',...
% 'mk_001(f)_boredom_passive_negative_22',...
% 'mk_001(f)_coldAnger_active_negative_13',...
% 'mm_001(f)_contempt_active_negative_7',...
% 'mm_001(f)_hotAnger_active_negative_16'};

% initialization
num_noise = length(noise_type1);
num_list = length(list);
num_SNR = length(SNR_dB1);

% initilize the average rate matrix for all files in a list
rate_bana_av_all1 = zeros(1, num_SNR);
rate_hps_av_all1 = zeros(1, num_SNR);
rate_yin_av_all1 = zeros(1, num_SNR);
rate_cep_av_all1 = zeros(1, num_SNR);

for l= 1:num_list
    audioname = list{l}; % one audio file name
    gtname = strcat(groundtruth_path, audioname, '.mat'); % ground truth MAT file name
    load(gtname);
    gt=y2; % load pitch groundtruth
    vm = zeros(1,length(gt)); % initialize the voice marker. 1: voiced frame, 0: unvoiced frame
    vm(gt~=0)=1; % in groundtruth, unvoiced frames are labeled as 0
    lenspeech = nnz(vm); % number of voiced frames
    
    % initialize the rate matrix for one audio file
    rate_bana = zeros(num_noise,num_SNR);
    rate_cep = zeros(num_noise,num_SNR);
    rate_hps = zeros(num_noise,num_SNR);
    rate_yin = zeros(num_noise,num_SNR);
    
    for j = 1:num_noise
        for k =1:num_SNR
            noisy_file = [audioname,'-',noise_type1{j},'_',num2str(SNR_dB1(k)),'dB.wav'];
            filename = strcat(noisyfile_path, noisy_file);
            [y,fs] = wavread(filename); % read wave file
            
            % BaNa
            f0_bana = Bana(y,fs,f0min,f0max,timestep,vm);
            frame_num = length(f0_bana);
            er = abs(f0_bana - gt); % difference between ground truth and the detected pitch
            hit = zeros(1,frame_num);
            for i = 1:frame_num
                if er(i)<= error_tolerance*gt(i) && vm(i) == 1 
                    hit(i) = 1;
                end
            end
            rate_bana(j,k) = sum(hit)/lenspeech;
            
            % Cepstrum
            [f0_cep,~] = cepstral(y, fs, f0min, f0max, timestep, pthr1);
            er = abs(f0_cep - gt); % difference between ground truth and the detected pitch
            hit = zeros(1,frame_num);
            for i = 1:frame_num
                if er(i)<= error_tolerance*gt(i) && vm(i) == 1
                    hit(i) = 1;
                end
            end
            rate_cep(j,k) = sum(hit)/lenspeech;
            
            % YIN
            % set input parameter structure
            P=[];
            P.minf0 = f0min;
%             P.maxf0 = f0max;
            P.hop = ceil(timestep*fs);
            P.range = [1 length(y)];
            P.wsize = ceil(time_per_window*fs);
            P.sr = fs;
            R=YIN(filename,P); % output R is a structure
            f0_yin = zeros(frame_num,1);
            for myp = 1: frame_num
                f0_yin(myp,1) = 440*(2^R.f0(myp)); % convert YIN's original output pitch from 440Hz octave to Hz
            end
            er = abs(f0_yin - gt); % difference between ground truth and the detected pitch
            hit = zeros(1,frame_num);
            for i = 1:frame_num
                if er(i)<= error_tolerance*gt(i) && vm(i) == 1 
                    hit(i) = 1;
                end
            end
            rate_yin(j,k) = sum(hit)/lenspeech;        
            
            % HPS
            wsize = ceil(time_per_window*fs);       % number of samples per window
            offset = ceil(timestep*fs);
            f0_hps = zeros(1,frame_num);
            for i = 1:frame_num
                temp = y(((i-1)*offset+1):((i-1)*offset+wsize));
                f0_hps(i) = hps(temp,fs);
            end
            er = abs(f0_hps - gt'); % difference between ground truth and the detected pitch
            hit = zeros(1,frame_num);
            for i = 1:frame_num
                if er(i)<= error_tolerance*gt(i) && vm(i) == 1 
                    hit(i) = 1;
                end
            end
            rate_hps(j,k) = sum(hit)/lenspeech;    
            
            % pitch detection of one wave file is finished for all
            % algorithms
            disp(['Pitch detection for file ', noisy_file, ' has finished.']);
        end
    end
    
    rate_bana_av = mean(rate_bana);
    rate_cep_av = mean(rate_cep);
    rate_hps_av = mean(rate_hps);
    rate_yin_av = mean(rate_yin);
    
    rate_bana_av_all1 = rate_bana_av_all1 + rate_bana_av;
    rate_hps_av_all1 = rate_hps_av_all1 + rate_hps_av;
    rate_yin_av_all1 = rate_yin_av_all1 + rate_yin_av;
    rate_cep_av_all1 = rate_cep_av_all1 + rate_cep_av;
end

rate_bana_av_all1 = rate_bana_av_all1/num_list;
rate_hps_av_all1 = rate_hps_av_all1/num_list;
rate_yin_av_all1 = rate_yin_av_all1/num_list;
rate_cep_av_all1 = rate_cep_av_all1/num_list;

fig1 = figure;
hold on;
algorithms = {'Bana', 'HPS', 'Yin','Cepstrum'};
linestyles = cellstr(char('--r','-k','-g','-c'));
markers=['s','d','>','x'];

plot(SNR_dB1,rate_bana_av_all1*100,[linestyles{1} markers(1)],'LineWidth',2,'MarkerSize',10);
plot(SNR_dB1,rate_hps_av_all1*100,[linestyles{2} markers(2)],'LineWidth',2,'MarkerSize',10);
plot(SNR_dB1,rate_yin_av_all1*100,[linestyles{3} markers(3)],'LineWidth',2,'MarkerSize',10);
plot(SNR_dB1,rate_cep_av_all1*100,[linestyles{4} markers(4)],'LineWidth',2,'MarkerSize',10);
legend('BaNa','HPS','YIN','Cepstrum','Location','SouthEast');
h1 = xlabel('SNR (dB)'); 
h2 = ylabel('Pitch Detection Accuracy (%)');
set(h1,'fontsize',14);
set(h2,'fontsize',14);
set(gca,'fontsize',14,'XTick',SNR_dB1);
ylim([20 100]);
grid on;

% save data
save(saveFilename);

