
function SortIVIM_Pelvis(varargin)

targetpath = varargin{1}; %'/Users/neuroimaging/Desktop/IVIM_PD/IVIM_13/Day1/
targetfolder = 'DICOM';
cd (fullfile(targetpath,targetfolder)) % for 13 & 14 


fprintf('sorting %s\n' ,fullfile(targetpath,targetfolder))

dat_list = dir(fullfile('*.dcm'));

SliceLocation = zeros(1,length(dat_list));
BVals = zeros(1,length(dat_list));

for j = 1:length(dat_list)
    h=dicominfo(char(fullfile(targetpath,targetfolder,dat_list(j).name)));
    SliceLocation(j) = h.SliceLocation;
    %BVals(j) = h.DiffusionBValue;
    bval_name = h.SequenceName;
    bval = extractBetween(bval_name,'*ep_b', 't');
    if isempty(bval)
        bval = extractAfter(bval_name,'*ep_b');
    end
    BVals(j) = string(bval);
end

Bvals = unique(BVals);
Slices = unique(SliceLocation);

num_b = length(Bvals);


sortdir = fullfile(targetpath,strcat(targetfolder,'_sorted'));
%if exist((sortdir),'dir')
    %fprintf('Done\n')
%else
    mkdir(sortdir)
    
    %size(dir(fullfile(sortdir,'/*.dcm')),1)
    
    %BVals = [111,222,333,444,556,667,778,889,1000];
    % The following is why I hate matlab. So inelegant. So sorry.


    for dcm_number = 1:length(Slices)*num_b
         h=dicominfo(char(fullfile(targetpath,targetfolder,dat_list(dcm_number).name)));
       % and now assign it a number
         [~,s] = min(abs(Slices-(h.SliceLocation))); % get slice index of this dicom's slice location
         bval_name = h.SequenceName;
         bval = extractBetween(bval_name,'*ep_b', 't');
         if isempty(bval)
             bval = {extractAfter(bval_name,'*ep_b')};
         end
         %bval
         [~,k] = min(abs(Bvals-(str2double(bval{1}))));

         %disp([h.SliceLocation,bval])
         %disp([dcm_number, s, k, (s-1)*num_b+k])
        % from this, resort and export
         name = string((s-1)*num_b + k); %new index
         name = char(strcat('IM_',name,'.dcm'));
         folder = char(strcat(targetfolder,'_sorted'));
         sourcefile = char(fullfile(targetpath,targetfolder,dat_list(dcm_number).name));
         destinationfile = fullfile(targetpath,folder,name);
         copyfile(sourcefile,destinationfile);
    end




    fprintf('Sorted %s\n' ,fullfile(targetpath))

end