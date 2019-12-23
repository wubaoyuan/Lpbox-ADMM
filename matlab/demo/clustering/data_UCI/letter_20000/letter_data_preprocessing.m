clear; clear all; 
chdir('D:\Latex\MLML all models\BP_with_ADMM\Applications_and_Compared_methods\clustering')
addpath(genpath(pwd))

%% read the original file, and transform the char label to numberical label
fileID = fopen('letter.txt');
C = textscan(fileID,'%s %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n','Delimiter',',','EmptyValue',-Inf);
fclose(fileID);

dataset_matrix = zeros(20000, 16); 
for i = 2:17
     dataset_matrix(:, i-1) = C{i};
end

character_list = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'};
label_gt_vec = zeros(20000, 1); 
parfor i = 1:20000
    label_gt_vec(i) = find(strncmp(character_list, C{1}(i), 3));
end

%% save 
dlmwrite('.\data_UCI\letter_20000\letter_only_data_20000_16.txt', dataset_matrix);
dlmwrite('.\data_UCI\letter_20000\letter_label_numerical_20000_1.txt', label_gt_vec);
dlmwrite('.\data_UCI\letter_20000\letter_label_char_20000_1.txt', C{1});