function selectivity_site = assign_site_selectivity(sig_struct, site_map_recenter_SC)

% note this function only is made to be used with the recentering dataset,
% since I did not track recording site on other mice/sessions.

selectivity_site = sig_struct.selectivity_cent;

for i = 1:size(selectivity_site, 1)
    
    % get session name - column corresponding to session name varies
    % depending on analysis.
    if size(selectivity_site, 2) <= 4
        ephys_session_name = selectivity_site{i, 3};
    else
        ephys_session_name = selectivity_site{i, 5};
    end
    ephys_date = ephys_session_name(1:6);
    ephys_mouse = ephys_session_name(8:17);

    % get recording site - manual fix for SC_Ephys_4 naming here
    if strcmp(ephys_mouse, 'SC_Ephys_4') && strcmp(ephys_date, '230805')
        ephys_date = '230804';
        site_map_ind = cell2mat(cellfun(@(x) contains(x, ephys_date) & contains(x, ephys_mouse, 'IgnoreCase', true), site_map_recenter_SC(:, 1), 'UniformOutput', false));
        site = site_map_recenter_SC{site_map_ind, 3};
    else
        site_map_ind = cell2mat(cellfun(@(x) contains(x, ephys_date) & contains(x, ephys_mouse, 'IgnoreCase', true), site_map_recenter_SC(:, 1), 'UniformOutput', false));
        site = site_map_recenter_SC{site_map_ind, 3};
    end
    if size(selectivity_site, 2) <= 4
        selectivity_site{i, 4} = site;
    else
        selectivity_site{i, 6} = site;
    end
end

end