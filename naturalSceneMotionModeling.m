%% This file collates pieces of MATLAB code (vR2018b) that were written by James E Fitzgerald (2018) and provided as supplementary material to:
% Salazar-Gatzimas E, Agrochao M, Fitzgerald JE, Clark DA (2018) Decorrelation of parallel motion pathways explains the neuronal basis of an illusory motion percept. Current Biology.
% Please read the paper for an explanation and interpretation of these analyses.
% Please send any questions to James Fitzgerald at fitzgeraldj@janelia.hhmi.org

% Because various routines save and load matfiles, you'll need to
% edit the code before it'll run on your machine. I've indicated lines that 
% must be changed with MUSTEDIT: Strictly speaking, a subset of these lines are
% self-referential and will run without editing, but you should edit
% them anyway to avoid having a bunch of opaquely named files written on your
% computer. Note that the largest file is about 12 GB in size. 

% BLOCK 1 generates the simulated ensemble of naturalistic motions using
% the van Hateren database of natural images.

% BLOCK 2 performs analyses on the simulated motions ensemble

%% BEGIN BLOCK 1
% NOTE: BLOCK 1 of code was originally published as supplementary material to:
% Fitzgerald JE, Clark DA (2015) Nonlinear circuits for naturalistic visual motion estimation. eLife, doi: 10.7554/eLife.09123.
% Please read this paper for further explanation and interpretation of these analyses.
% Please send any questions to James Fitzgerald at fitzgeraldj@janelia.hhmi.org

%% Start over and clear the workspace.

clear all;
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

filedir = 'VH_TOP_DIR/'; % MUSTEDIT: This is the directory where the van Hateren image database is stored
cd(filedir);
filenames = dir(strcat(filedir,'vanhateren_iml/*.iml')); % These are the images

% Image parameters.
width_pix = 1536;
width_deg = width_pix/60;
height_pix = 1024;
height_deg = height_pix/60;

% Camera settings
camera_settings = dlmread('camerasettings.txt','',1,0); % This file is available through the same source as the images
luminance_conversion = camera_settings(:,5);

%% Generate 1D database of natural images. 

tic
gaussian_filter = fspecial('gaussian',[1 height_pix],5.7*60/(2*sqrt(2*log(2))));
gaussian_filter = repmat(gaussian_filter,width_pix,1);
N = length(filenames);
oneD_nat_luminance = zeros(width_pix,N);
oneD_nat_contrast = zeros(width_pix,N);
for i = 1:N;
    if eq(mod(i,100),1);
        strcat('On image =',num2str(i))
        toc
    end;
    image_number = i;
    filename = strcat('vanhateren_iml/',filenames(i).name);
    f1=fopen(filename,'rb','ieee-be');
    raw=fread(f1,[width_pix,height_pix],'uint16');
    fclose(f1);
    % Transformed into luminance.
    luminance = luminance_conversion(image_number)*raw;
    oneD_nat_luminance(:,i) = sum(gaussian_filter.*luminance,2);
    % Transformed into contrast.
    I_0 = mean(luminance(:));
    contrast = (luminance - I_0)/I_0;
    oneD_nat_contrast(:,i) = sum(gaussian_filter.*contrast,2);
end;
oneD_nat_luminance = [oneD_nat_luminance;oneD_nat_luminance(end:-1:1,:)];
oneD_nat_contrast = [oneD_nat_contrast;oneD_nat_contrast(end:-1:1,:)];

%% Save 1D images. To avoid regenerating the oneD images every time

matfile = 'IMAGES1D.mat'; % MUSTEDIT: This is the one-dimensional database of natural images
save(matfile,'oneD_nat_luminance','oneD_nat_contrast','N');

%% Load 1D images. To avoid regenerating the oneD images every time

matfile = 'IMAGES1D.mat'; % MUSTEDIT: This is the one-dimensional database of natural images
load(matfile);

%% Generate sampled photoreceptor responses with left-right symmetry enforced.

tic;

da_sampling = 0.1;
PR_spacing = 5.1/da_sampling;
pix_per_PR = ceil(5.7*5/da_sampling);
num_samples = 1000000;
num_PRs = 12;
h_initial = 0.005;
h_final = 0.005;
h_factor = h_final/h_initial;

v_sd = 90;
T = 0.8;
psi_dot = zeros(num_samples,1);
image_idx = zeros(num_samples,1);
a_final_list = zeros(num_samples,1);
VBar_t = zeros(1+T/h_final,num_PRs,num_samples);
spacetime = zeros(1+T/h_initial,pix_per_PR,num_PRs);

gaussian_filter = fspecial('gaussian',[1 pix_per_PR],5.7/(da_sampling*2*sqrt(2*log(2))));
gaussian_filter = repmat(gaussian_filter,1+T/h_initial,1);

PR_tau = 0.01/h_initial;
filter_length = 5*PR_tau;
ex_filter = exp(-(filter_length:-1:0)/(PR_tau));
ex_filter = ex_filter/sum(ex_filter);
ex_filter = [ex_filter';zeros(filter_length,1)];

if gt(da_sampling,1/60);
    downsampling_factor = da_sampling*60;
    num_downsampled_pix = floor(3072/downsampling_factor);
    oneD_nat_contrast_sim = zeros(num_downsampled_pix,4167);
    for i = 1:num_downsampled_pix;
        oneD_nat_contrast_sim(i,:) = mean(oneD_nat_contrast((1+(i-1)*downsampling_factor):i*downsampling_factor,:),1);
    end;
elseif eq(da_sampling,1/60)
    downsampling_factor = 1;
    num_downsampled_pix = 3072;
    oneD_nat_contrast_sim = oneD_nat_contrast;
else
    'Invalid da_sampling'
    return;
end;

close all;
figure,
imagesc(oneD_nat_contrast')
colormap('gray')
figure,
imagesc(oneD_nat_contrast_sim')
colormap('gray')
%
pixel_vector = (0:pix_per_PR-1);
whos
for i = 1:num_samples/2;
    if eq(mod(i,5000),1)
        i
        toc;
    end;
    % Choose a natural image row.
    idx = ceil(rand*N);
    image_idx(i) = idx;
    oneD = oneD_nat_contrast_sim(:,idx);
    % Choose a velocity.
    v = -abs(v_sd*randn);
    psi_dot(i) = v;
    % Create spacetime plot
    a_final = ceil(rand*num_downsampled_pix);
	a_final_list(i) = a_final;
    for j = 1:num_PRs;
        a_f = a_final+(j-1)*PR_spacing;
        for t = 0:round(T/h_initial);
            pixel_shift = -v*t*h_initial/da_sampling;
            relevant_pixels = round(a_f-pixel_shift+pixel_vector);
            relevant_pixels = 1+mod(relevant_pixels,num_downsampled_pix);
            spacetime(round(T/h_initial)-t+1,:,j) = oneD(relevant_pixels);
        end;
        spatially_filtered = sum(spacetime(:,:,j).*gaussian_filter,2);
        tmp = imfilter(spatially_filtered,ex_filter);
        VBar_t(:,j,i) = tmp(1:h_factor:end);
    end;
    % Create symmetric spacetime plot
    v = -v;
    psi_dot(i+num_samples/2) = v;
    image_idx(i+num_samples/2) = idx;
    a_final = mod(-a_final-pix_per_PR-(num_PRs-1)*PR_spacing-1,num_downsampled_pix)+1;
    a_final_list(i+num_samples/2) = a_final;
    for j = 1:num_PRs;
        a_f = a_final+(j-1)*PR_spacing;
        for t = 0:round(T/h_initial);
            pixel_shift = -v*t*h_initial/da_sampling;
            relevant_pixels = round(a_f-pixel_shift+pixel_vector);
            relevant_pixels = 1+mod(relevant_pixels,num_downsampled_pix);
            spacetime(round(T/h_initial)-t+1,:,j) = oneD(relevant_pixels);
        end;
        spatially_filtered = sum(spacetime(:,:,j).*gaussian_filter,2);
        tmp = imfilter(spatially_filtered,ex_filter);
        VBar_t(:,j,i+num_samples/2) = tmp(1:h_factor:end);
    end;
end;

toc

%% Save matfile.

matfile = 'ENSEMBLE_DIR/ENSEMBLE.mat'; % MUSTEDIT: This is the simulated ensemble of natural images
save(matfile,'psi_dot','VBar_t','image_idx','a_final_list','-v7.3')

%% Load matfile.

matfile = 'ENSEMBLE_DIR/ENSEMBLE.mat'; % MUSTEDIT: This is the simulated ensemble of natural images
load(matfile)

%% END BLOCK 1

%% BEGIN BLOCK 2

%% Change to base directory for Drosophila-like spatiotemporal sampling

clear all;
base_dir = 'ENSEMBLE_DIR/'; % MUSTEDIT: This is the directory where you stored your simulated ensemble of motions
cd(base_dir);

%% Load file

matfile = strcat(base_dir,'ENSEMBLE.mat'); % MUSTEDIT: This is the simulated ensemble of natural images
load(matfile);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

%% Define the low-pass and high-pass HRC filters

clear('lp_filter','hp_filter');

h_final = 0.005;
tau = 0.02/h_final;
T = 0.8;
lp_filter = (T/h_final:-1:0).*exp(-(T/h_final:-1:0)/tau);
lp_filter = lp_filter/sum(lp_filter);
hp_filter(1+T/h_final) = 0;
hp_filter(1:T/h_final) = -diff(lp_filter);

%% Clean up workspace.

clearvars -except base_dir matfile psi_dot lp_filter hp_filter VBar_t

%% Set up subsets of the data.
% Error bars are estimated by considering how statistical properties vary
% across multiple random halves of the data.

num_sims = length(psi_dot);
num_trials = 5;
heldout_indices = NaN(num_trials,num_sims/2);
testing_indices = NaN(num_trials,num_sims/2);
for i = 1:num_trials;
    random_permutation = randperm(num_sims/2);
    heldout_indices(i,:) = [random_permutation(1:floor(num_sims/4)) num_sims/2+random_permutation(1:floor(num_sims/4))];
    testing_indices(i,:) = [random_permutation((floor(num_sims/4)+1):num_sims/2) num_sims/2+random_permutation((floor(num_sims/4)+1):num_sims/2)];
end;

%% Generate full set of filtered signals.

[tmp num_PRs num_sims] = size(VBar_t);
f_a = zeros(num_sims,num_PRs);
g_a = zeros(num_sims,num_PRs);
for n = 1:num_sims;
    for i = 1:num_PRs;
        f_a(n,i) = lp_filter*VBar_t(:,i,n);
        g_a(n,i) = hp_filter*VBar_t(:,i,n);
    end;
end;
f_a_p = f_a.*gt(f_a,0);
f_a_m = f_a.*lt(f_a,0);
g_a_p = g_a.*gt(g_a,0);
g_a_m = g_a.*lt(g_a,0);

%% Compute the full suite of interesting local motion detectors. 

HRC_local_a = f_a_p(:,1:(end-1)).*g_a_p(:,2:end) - g_a_p(:,1:(end-1)).*f_a_p(:,2:end) + f_a_m(:,1:(end-1)).*g_a_m(:,2:end) - g_a_m(:,1:(end-1)).*f_a_m(:,2:end) + f_a_p(:,1:(end-1)).*g_a_m(:,2:end) - g_a_m(:,1:(end-1)).*f_a_p(:,2:end) + f_a_m(:,1:(end-1)).*g_a_p(:,2:end) - g_a_p(:,1:(end-1)).*f_a_m(:,2:end);
phi_local_a = f_a_p(:,1:(end-1)).*g_a_p(:,2:end) - g_a_p(:,1:(end-1)).*f_a_p(:,2:end) + f_a_m(:,1:(end-1)).*g_a_m(:,2:end) - g_a_m(:,1:(end-1)).*f_a_m(:,2:end);

gC_local_a = f_a_p(:,1:(end-2)).*g_a_p(:,2:(end-1)) - g_a_p(:,2:(end-1)).*f_a_p(:,3:end) + f_a_m(:,1:(end-2)).*g_a_m(:,2:(end-1)) - g_a_m(:,2:(end-1)).*f_a_m(:,3:end) + f_a_p(:,1:(end-2)).*g_a_m(:,2:(end-1)) - g_a_m(:,2:(end-1)).*f_a_p(:,3:end) + f_a_m(:,1:(end-2)).*g_a_p(:,2:(end-1)) - g_a_p(:,2:(end-1)).*f_a_m(:,3:end);
gC_phi_local_a = f_a_p(:,1:(end-2)).*g_a_p(:,2:(end-1)) - g_a_p(:,2:(end-1)).*f_a_p(:,3:end) + f_a_m(:,1:(end-2)).*g_a_m(:,2:(end-1)) - g_a_m(:,2:(end-1)).*f_a_m(:,3:end);

%% Compute the full set of correlation accuracies.

corr_HRCvC_dist = NaN(num_trials,num_PRs-2,4);
for i = 1:num_trials;
    for d = 1:(num_PRs-2);
        corr_HRCvC_dist(i,d,1) = corr(sum(HRC_local_a(testing_indices(i,:),1:d),2),psi_dot(testing_indices(i,:)));
        corr_HRCvC_dist(i,d,2) = corr(sum(phi_local_a(testing_indices(i,:),1:d),2),psi_dot(testing_indices(i,:)));
        corr_HRCvC_dist(i,d,3) = corr(sum(gC_local_a(testing_indices(i,:),1:d),2),psi_dot(testing_indices(i,:)));
        corr_HRCvC_dist(i,d,4) = corr(sum(gC_phi_local_a(testing_indices(i,:),1:d),2),psi_dot(testing_indices(i,:)));
    end;
end;

%% Plot the performance measures (Fig. 4C)

close all;

figure,
hold on;
errorbar(mean(corr_HRCvC_dist(:,:,1),1),std(corr_HRCvC_dist(:,:,1),0,1),'-or','LineWidth',3)
errorbar(mean(corr_HRCvC_dist(:,:,2),1),std(corr_HRCvC_dist(:,:,2),0,1),'--or','LineWidth',3)
errorbar(mean(corr_HRCvC_dist(:,:,3),1),std(corr_HRCvC_dist(:,:,3),0,1),'-ob','LineWidth',3)
errorbar(mean(corr_HRCvC_dist(:,:,4),1),std(corr_HRCvC_dist(:,:,4),0,1),'--ob','LineWidth',3)
hold off;
legend('HRC','HRC phi only','shared non-delay','shared non-delay phi only')
xlabel('number of detectors')
ylabel('correlation with velocity')

mean_corr_HRCvC_dist = squeeze(mean(corr_HRCvC_dist,1));
std_corr_HRCvC_dist = squeeze(std(corr_HRCvC_dist,0,1));

%% Compute the coactivation statistics of the octant signals in the shared non-delay line model.

gC_local_1_components = [f_a_p(:,1).*g_a_p(:,2), g_a_p(:,2).*f_a_p(:,3), f_a_m(:,1).*g_a_m(:,2), g_a_m(:,2).*f_a_m(:,3), f_a_p(:,1).*g_a_m(:,2), g_a_m(:,2).*f_a_p(:,3), f_a_m(:,1).*g_a_p(:,2), g_a_p(:,2).*f_a_m(:,3)];
gC_local_1_components_coactivation = gC_local_1_components'*gC_local_1_components;
for i = 1:8;
    gC_local_1_components_coactivation(i,:) = gC_local_1_components_coactivation(i,:)/sum(abs(gC_local_1_components_coactivation(i,:)));
end;

%% Plot rows of the coactivation matrix. (Fig. 4D)

close all;
figure,
subplot(4,1,1)
bar(gC_local_1_components_coactivation(1,:));
title('pp prog')
subplot(4,1,2)
bar(gC_local_1_components_coactivation(2,:));
title('pp regressive')
ylabel('coactivation level (a.u)')
subplot(4,1,3)
bar(gC_local_1_components_coactivation(3,:));
title('mm progressive')
subplot(4,1,4)
bar(gC_local_1_components_coactivation(4,:));
title('mm regressive')
xlabel('octant index (1. pp prog, 2. pp reg, 3. mm prog, 4. mm reg, 5. pm prog, 6. pm reg, 7. mp prog, 8. mp reg)')

%% Compute the direction selective components of the octant signals in the shared non-delay line model. 

pp_DS = f_a_p(:,1).*g_a_p(:,2) - g_a_p(:,2).*f_a_p(:,3);
mm_DS = f_a_m(:,1).*g_a_m(:,2) - g_a_m(:,2).*f_a_m(:,3);
pm_DS = f_a_p(:,1).*g_a_m(:,2) - g_a_m(:,2).*f_a_p(:,3);
mp_DS = f_a_m(:,1).*g_a_p(:,2) - g_a_p(:,2).*f_a_m(:,3);

%% Define models for T4 and T5 neurons. 

T4_progressive = f_a_p(:,1).*g_a_p(:,2) - g_a_p(:,2).*f_a_m(:,3);
T4_regressive = g_a_p(:,2).*f_a_p(:,3) - f_a_m(:,1).*g_a_p(:,2);
T5_progressive = f_a_m(:,1).*g_a_m(:,2) - g_a_m(:,2).*f_a_p(:,3);
T5_regressive = g_a_m(:,2).*f_a_m(:,3) - f_a_p(:,1).*g_a_m(:,2);

neuronal_components = [T4_progressive, T4_regressive, T5_progressive, T5_regressive];
neuronal_components_coactivation = neuronal_components'*neuronal_components;

%% Visualize the direction-selective components model neuronal activations of naturalistic motion.
% A few illustrative examples are shown in Fig. 5A

example_naturalistic_motion = ceil(rand*num_sims);
close all;
figure,
subplot(2,1,1);
hold on;
bar(1,pp_DS(example_naturalistic_motion),'k');
bar(2,mm_DS(example_naturalistic_motion),'g');
bar(3,pm_DS(example_naturalistic_motion),'y');
bar(4,mp_DS(example_naturalistic_motion),'c');
xlabel('DS-component index (1. pp, 2. mm, 3. pm 4. mp)')
ylabel('direction selective response (a.u.)')
hold off;
subplot(2,1,2);
hold on;
bar(1,T4_progressive(example_naturalistic_motion),'r');
bar(2,T4_regressive(example_naturalistic_motion),'m');
bar(3,T5_progressive(example_naturalistic_motion),'b');
bar(4,T5_regressive(example_naturalistic_motion),'g');
hold off;
xlabel('channel index (1. T4 prog, 2. T4 reg, 3. T5 prog 4. T5 reg)')
ylabel('synthetic component response (a.u.)')

%% Visualize the channel coactivation matrix. (Fig. 5B)

close all;

figure, 
imagesc(neuronal_components_coactivation)
colormap('gray')
xlabel('channel index (1. T4 prog, 2. T4 reg, 3. T5 prog 4. T5 reg)')
ylabel('channel index (1. T4 prog, 2. T4 reg, 3. T5 prog 4. T5 reg)')

%% END BLOCK 2
