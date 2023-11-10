%% Export all the functions for the current model
%--------------------------------------------------------------------------
% Store all the terms describing the body also in a common folder ---------
%--------------------------------------------------------------------------
disp("Exporting functions...");

if ~exist('optimize_export','var')
    optimize_export = true;
end
cF = pwd;
cd(destination_folder_j);

disp("Replacing in the Joint folder the fresnel integrals");
myFiles = dir('*.m'); 
for i = 1:length(myFiles)
  baseFileName = myFiles(i).name;
  if ~contains(baseFileName, "_Joint")
        fprintf(1, 'Now reading %s\n', baseFileName);
        %Replace the fresnel integrals
        fid  = fopen(baseFileName,'r');
        f=fread(fid,'*char')';
        fclose(fid);
        f = strrep(f,'fresnelc','myfresnelc');
        f = strrep(f,'fresnels','myfresnels');
        f = strrep(f,'sqrt','mysqrt');
        fid  = fopen(baseFileName,'w');
        fprintf(fid,'%s',f);
        fclose(fid);
  end
end


cd(cF)
cd(destination_folder_b);

disp("Replacing in the Body folder the fresnel integrals");
myFiles = dir('*.m'); 
for i = 1:length(myFiles)
  baseFileName = myFiles(i).name;
  if ~contains(baseFileName, "_Body")
        fprintf(1, 'Now reading %s\n', baseFileName);
        %Replace the fresnel integrals
        fid  = fopen(baseFileName,'r');
        f=fread(fid,'*char')';
        fclose(fid);
        f = strrep(f,'fresnelc','myfresnelc');
        f = strrep(f,'fresnels','myfresnels');
        f = strrep(f,'sqrt','mysqrt');
        fid  = fopen(baseFileName,'w');
        fprintf(fid,'%s',f);
        fclose(fid);
  end
end

cd(cF)