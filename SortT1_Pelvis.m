% this runs all code necessary for T1 of the pelvis 

% Input: path to dicom folder and folder with the T1 B1 Corr maps in it
% output: sorted dicom folder for analysis
% Mira Liu 11/6/2025

function SortT1_Pelvis(varargin)

    targetpath = varargin{1}; 
    targetfolder = varargin{2};
    cd (fullfile(targetpath,targetfolder)) % for 13 & 14 
    
    
    fprintf('sorting %s\n' ,fullfile(targetpath,targetfolder))
    
    dat_list = dir(fullfile('*.dcm'));
    
    SliceLocation = zeros(1,length(dat_list));
    
    
    sortdir = fullfile(targetpath,strcat(targetfolder,'_sorted'));

    mkdir(sortdir)

    k = 0;
    for j = 1:length(dat_list)
        h=dicominfo(char(fullfile(targetpath,targetfolder,dat_list(j).name)));
        %if isAttribute(dicomFile(char(fullfile(targetpath,targetfolder,dat_list(j).name))),"SliceLocation")
        try
            SliceLocation(j) = h.SliceLocation;
            k=k+1;
            name = string(k); %new index
            name = char(strcat('IM_',name,'.dcm'));
            folder = char(strcat(targetfolder,'_sorted'));
            sourcefile = char(fullfile(targetpath,targetfolder,dat_list(j).name));
            destinationfile = fullfile(targetpath,folder,name);
            copyfile(sourcefile,destinationfile);
            
        catch
            %nothing

        end
    end
         
    fprintf('Sorted %s\n' ,fullfile(targetpath))

end