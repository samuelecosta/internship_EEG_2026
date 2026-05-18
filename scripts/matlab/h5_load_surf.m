function surf_struct = h5_load_surf(fileDir, pattern)
% LOAD_H5_SURFACES Finds and loads surfaces from .h5 files into a structure.
%
%   surf_struct = load_h5_surfaces() searches the current directory (pwd) 
%   for all files matching the pattern '*Surface_*.h5'.
%
%   surf_struct = load_h5_surfaces(fileDir) searches the specified directory 
%   'fileDir' using the default pattern.
%
%   surf_struct = load_h5_surfaces(fileDir, pattern) uses a custom directory 
%   and search pattern.
%
%   OUTPUT:
%   surf_struct: a structure with dynamically named fields. Each field
%                contains vertices, vertex_normals, triangles, and triangle_normals.

if nargin < 1 || isempty(fileDir)
    fileDir = pwd;
end
    
if nargin < 2 || isempty(pattern)
    pattern = '*Surface_*.h5';
end

h5_files = dir(fullfile(fileDir, pattern));

surf_struct = struct();

if isempty(h5_files)
    warning('No files found in %s with this pattern %s', fileDir, pattern);
    return; 
end

for i = 1:length(h5_files)
    file_path = fullfile(h5_files(i).folder, h5_files(i).name); 
        
    [~, base_name, ~] = fileparts(h5_files(i).name);
    pure_name = extractAfter(base_name, 'Surface_');
    entity_name = matlab.lang.makeValidName(pure_name);

    try
        vert = h5read(file_path,'/vertices');
        vert_norm = h5read(file_path,'/vertex_normals');
        tr = h5read(file_path,'/triangles');
        tr_norm = h5read(file_path,'/triangle_normals');

        surf_struct.(entity_name).vertices = vert / 1000;
        surf_struct.(entity_name).vertex_normals = vert_norm;
        surf_struct.(entity_name).triangles = tr;
        surf_struct.(entity_name).triangle_normals = tr_norm;

        centroid = mean(surf_struct.(entity_name).vertices, 2);
        vecs_from_center = surf_struct.(entity_name).vertices - centroid;
        dot_products = dot(vecs_from_center, vert_norm);

        if mean(dot_products) > 0
            dir_msg = 'Outward';
        else
            dir_msg = 'Inward';
        end

        fprintf('Entity loaded correctly: %s | Normals direction: %s\n', entity_name, dir_msg);
        
    catch ME
        warning('Failed to open file: %s\n%s', entity_name,ME.message);
    end
end

end