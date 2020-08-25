function structure = gt_sg_sub_append_to_struct(structure,new_element,ind)
%
% structure = gt_sg_sub_append_to_struct(structure,new_element,ind)
%
% GT_SG subfunction for concatenating the structure arrays when different 
% fields may be present. 
%
% Inputs:
% structure = Overall structure to which data will be added
% new_element = Data to add
% ind = Where to add the new data
% 
% Outputs:
% structure = Overall structure with concatenated new data
%
% B.Y.QUESTE Feb 2015

% Orderfields is necessary as even if the same fields are
% present, it bugs if they're not in the same order.
if ~isstruct(new_element)
    return;
end

new_element = orderfields(new_element);

try
    % Append data from new dive. Works if they contain the same
    % fieldnames. Otherwise move to catch statement.
    structure(ind) = new_element;
catch
    try
        % If fieldnames are different, collect fieldnames of both
        % preloaded structure and current dive, then initialise missing ones
        % with empty arrays in respective structures before merging
        % again. When initialising in pre-existing structure array, this
        % is done in the 'disposable' element: structure(numel(structure) + 1).
        % This element is later deleted.
        fields_pre = fieldnames(structure);
        fields_new = fieldnames(new_element);
        
        missing_new = intersect(setxor(fields_pre,fields_new),fields_pre);
        if ~isempty(missing_new)
            for fstep = 1:numel(missing_new)
               new_element.(missing_new{fstep}) = []; 
            end
            new_element = orderfields(new_element);
        end
        
        missing_pre = intersect(setxor(fields_pre,fields_new),fields_new);
        if ~isempty(missing_pre)
            % Creates one additional element in the structure array to create a
            % disposable element where fieldnames are set to [] allowing the dynamic
            % addition of fields and avoiding the 'Subscripted assignment between
            % dissimilar structures' error. Then delete it. 
            % Meh, this is ugly, but it's surprisingly fast and flexible...
            count = numel(structure);
            structure(count+1).(fields_new{1}) = [];
            for fstep = 1:numel(missing_pre)
               structure(count+1).(missing_pre{fstep}) = [];
            end
            structure = orderfields(structure(1:count));
        end
        
        structure(ind) = new_element;
    catch
        gt_sg_sub_echo({['ERROR: Could not process dive ',num2str(ind)],'Proceeding to load the rest of the data,'});
    end
end


end