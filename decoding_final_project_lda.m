%% Probabilistic & Statistical Modelling (II) Final Project
% Done by: Ayesha Ejaz
% Topic: Multi-model, Multivariate Analysis of Tactile Mental Imagery in Primary Somatosensory Cortex

%% The Decoding Toolbox (TDT) was used for performing LDA. 

%% First, setting the defaults and defining the analysis 

clearvars;
clc;
data_dir = 'C:\Users\Medion\Downloads\Decoding_project';
addpath(data_dir) % Adding path to all subject folders
addpath('C:\Users\Medion\Downloads\tdt_3.999E2\decoding_toolbox') % Adding path to tdt toolbox
assert(~isempty(which('decoding_defaults.m', 'function')), 'TDT not found in path, please add')
addpath('C:\spm12\spm12') % Adding path to SPM 12
assert((~isempty(which('spm.m', 'function')) || ~isempty(which('BrikInfo.m', 'function'))) , 'Neither SPM nor AFNI found in path, please add (or remove this assert if you really dont need to read brain images)')

clear cfg

for subj = 1:10 % Iterating over all subjects

decoding_defaults; % Setting the defaults

cfg.analysis = 'ROI'; % Defining the analysis method
cfg.searchlight.radius = 3; % The unit by default is voxels
  
cfg.results.dir = sprintf('Decoding_project/results_ImagPress_vs_ImagVibro/sub-%03d', subj); % Defining where results would be saved
cfg.results.overwrite = 1 % In case you're running the analysis again and don't want the previous results, otherwise set to 0.

%% Second, getting the file names, labels and run number of each brain image
%% file to use for decoding.


subjfolder = sprintf('/sub-%03d/1st_level_good_bad_Imag/', subj); % Indicating path to subject folders

beta_loc = fullfile(data_dir, subjfolder); % Defining the path to SPM.mat and all related beta files

labelname1 = 'ImagPress'; % Specifying the label names that were given to the regressors of interest. This was rewritten two more times to perform pairwise classification between ImagPress, ImagFlutt and ImagVibro.
labelname2 = 'ImagVibro';


cfg.files.mask = {'C:\Users\Medion\Downloads\Decoding_project\rois\rPSC_1_TR50_right_CUT_Stim_vs_Null.nii', 'C:\Users\Medion\Downloads\Decoding_project\rois\rPSC_2_TR50_right_CUT_Stim_vs_Null.nii', 'C:\Users\Medion\Downloads\Decoding_project\rois\rPSC_3b_TR50_right_CUT_Stim_vs_Null.nii', 'C:\Users\Medion\Downloads\Decoding_project\rois\rSII_TR50_left_CUT_Stim_vs_Null.nii', 'C:\Users\Medion\Downloads\Decoding_project\rois\rSII_TR50_right_CUT_Stim_vs_Null.nii'}
% Multiple ROI masks from regions BA1, BA2, BA3b, ipsilateral S2 and contralateral S2 to define which voxels to use in the analysis

regressor_names = design_from_spm(beta_loc); % Extracting all beta names and corresponding run numbers from the SPM.mat file


cfg = decoding_describe_data(cfg,{labelname1 labelname2},[1 -1],regressor_names,beta_loc); % Extracting the file names and run numbers of each label


%% Third, creating the design for decoding analysis


 cfg.design = make_design_cv(cfg); % Creating a leave-one-run-out cross-validation design


%% Fourth, setting the additional parameters manually


cfg.verbose = 1; % How much output you want to see on the screen while the program is running. In this case, 0 = no output.

cfg.decoding.method = 'classification'; % This is our default anyway.

cfg.decoding.train.classification.model_parameters.shrinkage = 'lw2'; % Specific for lda in TDT

cfg.results.output = {'accuracy_minus_chance'};

cfg.decoding.software = 'lda'; 

%% Fifth, plotting


cfg.plot_selected_voxels = 0; % In this case, plot nothing online

cfg.plot_design = 1; % This will call display_design(cfg);

display_design(cfg); % Allows you to look at your design after plotting

%% Sixth, Running the decoding analysis

results = decoding(cfg);
end 

%% Seventh, performing group Analysis

%Getting all the pairwise classification accuracies in one matrix for all
%subjects with each row being one ROI and each column a separate subject.

subjects = [1 2 3 4 5 6 7 8 9 10];
Group_Results_ImagPress_vs_ImagFlutt = [];
Group_Results_ImagPress_vs_ImagVibro = [];
Group_Results_ImagFlutt_vs_ImagVibro = [];



for subject=subjects

    subject = num2str(subject);

    ROI_Results1 = load([pwd '\results_ ImagPress_vs_ImagFlutt\sub-00' subject '/res_accuracy_minus_chance.mat']);
    ROI_Results2 = load([pwd '\results_ImagPress_vs_ImagVibro\sub-00' subject '/res_accuracy_minus_chance.mat']);
    ROI_Results3 = load([pwd '\results_ImagFlutt_vs_ImagVibro\sub-00' subject '/res_accuracy_minus_chance.mat']);

    Group_Results_ImagPress_vs_ImagFlutt = [Group_Results_ImagPress_vs_ImagFlutt, ROI_Results1.results.accuracy_minus_chance.output];
    Group_Results_ImagPress_vs_ImagVibro = [Group_Results_ImagPress_vs_ImagVibro, ROI_Results2.results.accuracy_minus_chance.output];
    Group_Results_ImagFlutt_vs_ImagVibro = [Group_Results_ImagFlutt_vs_ImagVibro, ROI_Results3.results.accuracy_minus_chance.output];

end

% Calculating the mean across all subjects for the five ROIs
ImagPress_vs_ImagFlutt_mean = mean(Group_Results_ImagPress_vs_ImagFlutt, 2);
ImagPress_vs_ImagVibro_mean = mean(Group_Results_ImagPress_vs_ImagVibro, 2);
ImagFlutt_vs_ImagVibro_mean = mean(Group_Results_ImagFlutt_vs_ImagVibro, 2);

% Getting the index of the ROI for each subject in which the decoding
% accuracy was highest
[~, ImagPress_vs_ImagFlutt_highest]= max(Group_Results_ImagPress_vs_ImagFlutt, [], 1);
[~,ImagPress_vs_ImagVibro_highest] = max(Group_Results_ImagPress_vs_ImagVibro, [], 1);
[~,ImagFlutt_vs_ImagVibro_highest] = max(Group_Results_ImagFlutt_vs_ImagVibro, [], 1);

% Getting the index of the ROI in which the accuracy was high for the most
% subjects
ImagPress_vs_ImagFlutt_highest_mode= mode(ImagPress_vs_ImagFlutt_highest);
ImagPress_vs_ImagVibro_highest_mode = mode(ImagPress_vs_ImagVibro_highest);
ImagFlutt_vs_ImagVibro_highest_mode = mode(ImagFlutt_vs_ImagVibro_highest);

% Performing ttest between the pairwise classifiers
[h, p, ci, stats] = ttest2(ImagPress_vs_ImagVibro_mean, ImagPress_vs_ImagFlutt_mean); % Pairwise classifiers ttest-1
disp("T-test results:");
disp("Hypothesis test result (h): " + h);
disp("p-value (p): " + p);
disp("Confidence interval (ci): " + ci);
disp("T-test statistics (stats):");
disp(stats);


[h, p, ci, stats] = ttest2(ImagPress_vs_ImagVibro_mean, ImagFlutt_vs_ImagVibro_mean); %Pairwise classifiers ttest-2
disp("T-test results:");
disp("Hypothesis test result (h): " + h);
disp("p-value (p): " + p);
disp("Confidence interval (ci): " + ci);
disp("T-test statistics (stats):");
disp(stats);

[h, p, ci, stats] = ttest2(ImagPress_vs_ImagFlutt_mean, ImagFlutt_vs_ImagVibro_mean); % Pairwsie classifiers ttest-3
disp("T-test results:");
disp("Hypothesis test result (h): " + h);
disp("p-value (p): " + p);
disp("Confidence interval (ci): " + ci);
disp("T-test statistics (stats):");
disp(stats);

% ttest for each roi across all subjects
% Creating data matrices
dataMatrix1 = Group_Results_ImagPress_vs_ImagFlutt;
dataMatrix2 = Group_Results_ImagPress_vs_ImagVibro;
dataMatrix3 = Group_Results_ImagFlutt_vs_ImagVibro;

% Creating a cell array to hold data matrices for the comparisons
dataMatrices = {dataMatrix1, dataMatrix2, dataMatrix3};

% Initializing a results structure
results2 = struct();

% Performing t-test for each row separately and each pairwise comparison
numComparisons = length(dataMatrices);
numRows = size(Group_Results_ImagPress_vs_ImagVibro, 1);

for comparisonIndex1 = 1:numComparisons
    dataMatrix1 = dataMatrices{comparisonIndex1};
    
    for comparisonIndex2 = comparisonIndex1+1:numComparisons
        dataMatrix2 = dataMatrices{comparisonIndex2};
        
        for rowIndex = 1:numRows
            value1 = dataMatrix1(rowIndex, :);
            value2 = dataMatrix2(rowIndex,:);
            
            % Performing t-test
            [h, p, ci, stats] = ttest2(value1, value2);
           
            % Storing results in the structure
            results2(end+1).comparisonIndex1 = comparisonIndex1;
            results2(end).comparisonIndex2 = comparisonIndex2;
            results2(end).rowIndex = rowIndex;
            results2(end).h = h;
            results2(end).p = p;
            results2(end).ci = ci;
            results2(end).stats = stats;
        end
    end
end
% The first row is empty for some reason, so just removing that
results2(1) = [];


%Bonferroni Correction
pValues = [results2.p]
numComparisons = numel(pValues);

% Desired family-wise error rate (usually 0.05)
desiredFWER = 0.05;

% Calculating adjusted alpha (Bonferroni-adjusted significance level)
adjustedAlpha = desiredFWER / numComparisons;

% Determining which tests are significant after Bonferroni correction
significantTests = pValues < adjustedAlpha;


% Creating an error barplot

categories = {'1 v 2', '1 v 2', '1 v 2', '1 v 2', '1 v 2', ...
              '1 v 3', '1 v 3', '1 v 3', '1 v 3', '1 v 3', ...
              '2 v 3', '2 v 3', '2 v 3', '2 v 3', '2 v 3'};

% Definining unique categories
uniqueCategories = unique(categories);
numCategories = numel(uniqueCategories);

% Definining colors for each category
categoryColors = {"#A2142F", "#77AC30", "#0072BD"};  % Adjust colors as needed

% Creating a bar graph with different colors for each category
figure;
hold on;

% Looping through unique categories and setting the color for each category
for i = 1:numel(uniqueCategories)
    categoryData = pValues(strcmp(categories, uniqueCategories{i}));
    bar(find(strcmp(categories, uniqueCategories{i})), categoryData, 'FaceColor', categoryColors{i});
end

hold off;
% Setting x-axis tick labels as well as x and y labels and the title
xticks(1:numel(data));
xticklabels(categories);
ylabel('P-value');
xlabel('Categories of Pairwise Classification');
title('P-values of Pairwise T-Tests');

% Rotating x-axis labels for better readability
xtickangle(45);

% Extracting the confidence intervals from the results structure and
% putting it in a suitable format for including an error bar
ci_results = [results2.ci];
ciLower = ci_results(:, [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29]);
ciUpper = ci_results(:, [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30]);
ciUpper = (ciUpper)';
ciLower = (ciLower)';
error = (ciUpper - ciLower)/3.92; % 3.92 for 95% CI
hold on;
errorbar(1:numel(pValues), pValues, error, 'k', 'LineStyle', 'none', 'CapSize', 0);
legend('P-value 1 v 2', 'P-value 1 v 3', 'P-value 2 v 3', 'error');
hold off;

%% The end. 
