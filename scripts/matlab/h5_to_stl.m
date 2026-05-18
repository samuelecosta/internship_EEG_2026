function h5_to_stl(surf_struct, outputDir)
% H5_TO_STL Exports 3D surface structures to .stl files.
%
%   h5_to_stl(surf_struct) exports each surface in the structure
%   to a folder named 'exported_stl' in the current working directory.
%
%   h5_to_stl(surf_struct, outputDir) exports the surfaces to the 
%   specified directory 'outputDir'.
%
%   INPUTS:
%   surf_struct: A structure where each field represents a 3D entity and 
%                contains at least 'vertices' and 'triangles'.
%   outputDir:   (Optional) Destination folder path.

if nargin < 1 || ~isstruct(surf_struct)
    error('Not a valid structure.');
end
        
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(pwd, 'exported_stl');
end
        
if ~exist(outputDir, 'dir')
    [status, msg] = mkdir(outputDir);
    if ~status
        error('Impossible to create the folder %s: %s', outputDir, msg);
    end
end
        
surf_names = fieldnames(surf_struct);
        
if isempty(surf_names)
    warning('Structure is empty.');
    return; 
end

for i = 1:length(surf_names)
    entity_name = surf_names{i};
    temp_surf = surf_struct.(entity_name);
            
    file_path = fullfile(outputDir, ['Surface_', entity_name, '.stl']); 
            
    try
        if ~isfield(temp_surf, 'vertices')
            error('Vertices missing.');
        end
        if ~isfield(temp_surf, 'triangles')
            error('Triangles missing.');
        end

        verts = double(temp_surf.vertices);

        faces = double(temp_surf.triangles) + 1; 

        if size(verts, 2) ~= 3
                verts = verts';
        end
        if size(faces, 2) ~= 3
                faces = faces';
        end

        [faces, was_flipped] = auto_fix_winding(verts, faces);
        if was_flipped
            fprintf(' -> [Auto-Fix] normals direction for: %s\n', entity_name);
        end

        TR = triangulation(faces, verts);

        stlwrite(TR, file_path);

        fprintf('Entity exported correctly: Surface_%s.stl\n', entity_name);
                
        catch ME
            warning('Failed to export entity: %s\n%s', entity_name, ME.message);
    end
end
end