
%% created for importing RT-struct to Matlab for data analysis.
% ML 08112025

% available at: 

% forked from https://github.com/viswanath-lab/Read-RTSTRUCT-to-matlab 

%% input: 
% path to study folder, study folder name, subfolder with dicom images and Rtstruct, and name of the contour of interest
%% output: 
% images as a stack, masks as a stack, and the contours as a cell array 


function [images, mask, contours] = ImportRTstruct(targetpath, Study_folder,Sub_folder,ROI_name)

    
    directory = (fullfile(targetpath, Study_folder, Sub_folder));
    files_out = scanDir(directory);
    fileList = dir(fullfile(directory,'*.mat'));
    loaded = load(fullfile(directory, fileList.name), 'contours');
    contours = loaded.contours;
    [images, ~] = dicomfolder(fullfile(targetpath, Study_folder, 'Dicoms'));
    
    ROI_names = {};

    for j = 1:size(contours,2)
	    ROI_names{end+1}=  contours(j).ROIName;
    end

    
    exact_match_mask = strcmp(ROI_names, ROI_name);
    n = find(exact_match_mask);
    mask = contours(n).Segmentation;
    mask = flip(permute(mask,[2, 1, 3]),3); %for aligning mask with the image stack
    %imagestack(1000*mask + images)
    
end
