function h5_to_obj(surf_struct, outputDir)
% EXPORT_SURF_TO_OBJ Exports 3D surface structures to .obj files.
%
%   h5_to_obj(surf_struct) exports each surface in the structure
%   to a folder named 'exported_obj' in the current working directory.
%
%   h5_to_obj(surf_struct, outputDir) exports the surfaces to the 
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
    outputDir = fullfile(pwd, 'exported_obj');
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
        
    file_path = fullfile(outputDir, ['Surface_', entity_name, '.obj']); 
        
    try
        fid = fopen(file_path, 'w');
        if fid == -1
            error('Could not open file .');
        end
            
        if isfield(temp_surf, 'vertices')
            verts = temp_surf.vertices;
            if size(verts, 2) ~= 3, verts = verts'; end
            fprintf(fid, 'v %f %f %f\n', verts');
        else
            error('Vertices missing.');
        end
            
        if isfield(temp_surf, 'vertex_normals')
            fprintf(fid, 'vn %f %f %f\n', temp_surf.vertex_normals);
            has_normals = true;
        else
            has_normals = false;
        end
            
        if isfield(temp_surf, 'triangles')
            F = double(temp_surf.triangles + 1);
            if size(F, 2) ~= 3, F = F'; end

            [F, was_flipped] = auto_fix_winding(verts, F);
            if was_flipped
                fprintf(' -> [Auto-Fix] Winding invertito per: %s\n', entity_name);
            end
                
            if has_normals
                face_data = [F(1,:); F(1,:); F(2,:); F(2,:); F(3,:); F(3,:)];
                fprintf(fid, 'f %d//%d %d//%d %d//%d\n', face_data);
            else
                fprintf(fid, 'f %d %d %d\n', F);
            end
        else
                 error('Triangles missing.');
        end
            fclose(fid);
            fprintf('Entity exported correctly: %s.obj\n', entity_name);
            
    catch ME
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        warning('Failed to export entity: %s\n%s', entity_name, ME.message);
    end
end
end