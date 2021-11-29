% Payam Farahani
% convert .txt file of cellpose outlines to tiff file of cell masks

files=uigetfile('*.txt');
nf = length(files);

fid = fopen(files);
tline = fgetl(fid);
tlines = cell(0,1);

 while ischar(tline)
        tlines{end+1,1} = tline;
        tline = fgetl(fid);
 end
 
  fclose(fid);
  coords=cell(length(tlines),1);
  areas=zeros(length(tlines),1);
  
  for j=1:length(tlines)
        temp=strsplit(tlines{j},',');
        temp2=reshape(temp,2,length(temp)/2)';
        coords{j}=cellfun(@str2num, temp2) + 1;
        areas(j)=polyarea(coords{j}(:,1),coords{j}(:,2)) * dim_xy^2;
        
        data(j).coords = coords{j};
        data(j).areas = areas(j);
  end
  
  
  [rows cols] = size(EGF);
  
  masks = zeros(rows,cols);
  
  for j = 1:length(coords)
      temp = coords{j};
      for k = 1:length(temp(:,1))
          masks(temp(k,2)+1,temp(k,1)+1) = j;
      end
  end
  
  imshow(masks)
  
% generate test file
test_file_name_mask = sprintf('mask_seg.tif');
imwrite(masks,test_file_name_mask,'WriteMode','append','compression','none')

  
  
  
