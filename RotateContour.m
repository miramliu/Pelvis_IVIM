function Rotated_RTstruct = RotateContour(RTstruct, img)
%test = dicominfo('/Users/miraliu/Downloads/DICOM_IVIM_export/2025-02__Studies/STAR-01^0002_STAR-01-0002_RTst_2025-02-05_175946_RESEARCH-20.GCO-20-2369-MARSHAL/2.16.840.1.114362.1.12289667.25143718253.710389707.526.194.dcm');
%RTstruct = dicomContours(test);
%img = dicomread('/Users/miraliu/Downloads/DICOM_IVIM_export/2025-02__Studies/STAR-01^0002_STAR-01-0002_MR_2025-02-05_175946_RESEARCH-20.GCO-20-2369-MARSHALL/2.16.840.1.114362.1.12289667.26020642722.696635700.517.4811.dcm');
%plotContour(RTstruct)

% number of anatomic volumes of interest
%size(RTstruct.ROIs,1)

%% number of anatomic volumes of interest
%size(RTstruct.ROIs.ContourData,1)

%% number of slices with an ROI for this contour
%size(RTstruct.ROIs.ContourData{contour_idx})

%% number of points for the ROI for this contour on this slice
%size(RTstruct.ROIs.ContourData{contour_idx}{slice})

%create Rotated_RTstruct
Rotated_RTstruct = RTstruct;
[nx, ny] = size(img);
figure,
for contour_idx = 1:size(RTstruct.ROIs.ContourData,1) % for each region

    for slice = 1:size(RTstruct.ROIs.ContourData{contour_idx},1) % for each slice
        ROI_points = RTstruct.ROIs.ContourData{contour_idx}{slice}; % get the points
        Rotated_ROI_points = squeeze(flip(ROI_points,2)); % rotte them
        %Rotated_RTstruct.ROIs.ContourData{contour_idx}{slice} = Rotated_ROI_points; %% read only error because Matlab
        %add to new cell array 
        Rotated_contours{slice} = Rotated_ROI_points;
        geometric_types{slice} = RTstruct.ROIs.GeometricType{contour_idx}{slice};
        
        %make mask for single slice
        ROI_points_slicelevel = ROI_points;
        %shift all to be positive because MATLABBABAOUHILKEJG
        ROI_points_slicelevel(:,1) = ROI_points_slicelevel(:,1)+100;
        ROI_points_slicelevel(:,2) = ROI_points_slicelevel(:,2)+100;
        ROI_points_slicelevel(:,3) = []; %set remove z axis
        imshow(poly2mask(ROI_points_slicelevel(:,1),ROI_points_slicelevel(:,2),nx+20,ny+20));
        pause(0.1)

     name = Rotated_RTstruct.ROIs.Name{contour_idx};
     color = Rotated_RTstruct.ROIs.Color{contour_idx};
     Rotated_RTstruct = deleteContour(Rotated_RTstruct,contour_idx);
     Rotated_RTstruct = addContour(Rotated_RTstruct,contour_idx,name,Rotated_contours,geometric_types,color);
    end
end

end
