%% H3K.m
% created by Prof. Julian Sosnik for Prof. Catheirne McuCusker
% 
%This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 international
%To view a copy of this license, visit https://creativecommons.org/licenses/by-sa/4.0/
%
%This program performs Quantificatoin of H3K27Me3 based on antibody
%staining (total fluorescence) inside nuclei stained with DAPI. The DAPI
%staining is used to segment each nuclei (in 3D) and the total fluorescence
%on the second chanel (H3K27Me3) is calculated
%
%It performs analysis of 
%consecutive Z slices from a stack of images.
%It also handles multiple image folders.

clear
clc 
 
SE = strel('square',2);
goback = pwd;
lowThresh = 1000;    % minimal nuclear size
highTrhesh = 20000;  % maximal nuclear size

%% Main routine starts here
directory = uigetdir;
cd (directory);                 %locates the folder to be analyzed
dirnames = dir(directory);


j = 0;                          %generates list of sub-directories
for i = 1:numel(dirnames)
    if dirnames(i).isdir == 1 && dirnames(i).name(1) ~='.'
    j = j+1;
    dirname(j) = dirnames(i);
    end
end

output = strcat(directory, '/output.csv');
fileID = fopen(output, 'w');
fprintf(fileID, '%s, %s, %s, %s, %s, %s, %s\n', 'File', 'Nucleus#',...
    'Nuc_Fluorescence', 'Nuc_Size', 'Total_Nuclei', 'Total_Fluo', 'AveFluo');
fclose(fileID);
    
s = numel(dirname);           %generate empty results matrices
NumNuclei = zeros(s,1);
TotalFluo = zeros(s,1);
AveFluo = zeros(s,1);

%% Initiates main loop directory by directory
for i = 1:numel(dirname)     
    cd (dirname(i).name);
    
    dirdir1 = strcat('C1_', dirname(i).name);
    dirdir2 = strcat('C2_', dirname(i).name);
    
    cd (dirdir1);
    C1filenames = dir('*.tif');   %generates list of files to process
    for j = numel(C1filenames):-1:1
        if C1filenames(j).name(1) == '.'
            C1filenames(j) = [];
        end
    end
    d = size(imread(C1filenames(1).name));
    cd ..;
    cd (dirdir2);
    C2filenames = dir('*.tif');   %generates list of files to process
    for j = numel(C2filenames):-1:1
        if C2filenames(j).name(1) == '.'
            C2filenames(j) = [];
        end
    end
    cd ..;

    I1 = zeros(d(1), d(2), numel(C1filenames));
    I1_6 = I1;
    I2 = I1;
    I2_0 = I2;
    
 %% Start of the file analysis loop (file by file in directory)
 % this loop loads up the matrix to be analyzed
    for j = 1:numel(C1filenames)
        
    %% Analysis of the image begins here
        cd (dirdir1);
        I1_0 = imread(C1filenames(j).name);
        cd ..;
        cd (dirdir2);
        I2_0(:,:,j) = imread(C2filenames(j).name);
        cd ..;
        I3_0 = I1_0;
        T = multithresh(I1_0);
        
        %Thereshold application to I1
        I1_1 = imgaussfilt(I1_0);
        I1_2 = imquantize(I1_1,T);
        I1_3 = imfill(I1_2, 'holes');
        I1_4 = logical(I1_3 - 1);
        I1_5 = imerode(I1_4, SE);
        
        I1_6(:,:,j) = I1_5;   % Make a 3d image of the individual layers
        
    end
    
        
        I1 = bwareaopen(I1_6,lowThresh) & ...
            ~bwareaopen(bwareaopen(I1_6,lowThresh),highTrhesh); % Threshshold for object size
        
                
        %Extract the nucleus from the original images to a black background
        for j = 1:numel(C1filenames)
            for k = 1:d(1)
                for l = 1:d(2)
                    if I1(k,l,j)>0
                        I2(k,l,j) = I2_0(k,l,j);
                    end
                end   
            end
        end   % This is the end of the loop loading the matrices
    
    %%Actual analysis start here
    
    CC = bwconncomp(I1,26);
    fluo = zeros(1,CC.NumObjects);
    numPixels = zeros(1,CC.NumObjects);
    %% calculate and save (append) results to file
    fileID = fopen(output, 'a');
    
    mapmat = zeros(size(I1));
    fluomat = zeros(size(I1));
    for m = 1:CC.NumObjects
        fluo(m) = sum(I2(CC.PixelIdxList{m}));
        mapmat(cell2mat(CC.PixelIdxList(:,m)))=m;
        fluomat(cell2mat(CC.PixelIdxList(:,m)))=fluo(m);
        numPixels(m) = cellfun(@numel,CC.PixelIdxList(m));
        NumNuclei(i) = CC.NumObjects;
        TotalFluo(i) = sum(I2(:));
        AveFluo(i) = TotalFluo(i) / NumNuclei(i);
        fprintf(fileID, '%s, %2d, %12.0f, %d, %8.0f, %12.0f, %12.2f\n',...
            dirname(i).name, m, fluo(m), numPixels(m), NumNuclei(i), TotalFluo(i), AveFluo(i));
    end
    
fclose(fileID);
save(output(1:end-4));
cd ..;   

end

cd(goback)