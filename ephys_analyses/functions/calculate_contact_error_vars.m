function contact_vars_struct = calculate_contact_error_vars(sp_vars_struct, fiducial_xy)

contact_centroid_filt = {sp_vars_struct.contact_centroid_filt};
tip_xf_prot_filt = {sp_vars_struct.tip_xf_prot_filt};
tip_yf_prot_filt = {sp_vars_struct.tip_yf_prot_filt};
tip_xf_cont_filt = {sp_vars_struct.tip_xf_cont_filt};
tip_yf_cont_filt = {sp_vars_struct.tip_yf_cont_filt};
centroid_xf_cont_filt = {sp_vars_struct.centroid_xf_cont_filt};
centroid_yf_cont_filt = {sp_vars_struct.centroid_yf_cont_filt};

contact_centroid_filt_mask = cell(size(sp_vars_struct));
prot_angle = cell(size(sp_vars_struct));
prot_angle_change = cell(size(sp_vars_struct));
lat_disp_change = cell(size(sp_vars_struct));
ml_error = cell(size(sp_vars_struct));
ml_error_ang = cell(size(sp_vars_struct));
contact_ang = cell(size(sp_vars_struct));
contact_ang_hc = cell(size(sp_vars_struct));
tip_contact_vect_mag = cell(size(sp_vars_struct));
tip_contact_vect_dist = cell(size(sp_vars_struct));
tip_contact_vect_ang = cell(size(sp_vars_struct));
tip_contact_vect_sign = cell(size(sp_vars_struct));
tip_contact_mag = cell(size(sp_vars_struct));
tip_contact_dist = cell(size(sp_vars_struct));
tip_contact_ang = cell(size(sp_vars_struct));

for i = 1:numel(tip_xf_cont_filt)
    
    % convert contact_centroid_filt to mask coordinates
    contact_centroid_filt_mask{i} = conv_cont_centroid_to_mask_coords(contact_centroid_filt{i}, fiducial_xy(i, :));

    contact_centroid_ml = cell2mat(cellfun(@(x) x(1), contact_centroid_filt_mask{i}, 'UniformOutput', false));
    contact_centroid_ap = cell2mat(cellfun(@(x) x(2), contact_centroid_filt_mask{i}, 'UniformOutput', false));
    
    % calculate protrusion angle
    prot_angle{i} = asind(tip_xf_prot_filt{i}./sqrt(tip_xf_prot_filt{i}.^2 + tip_yf_prot_filt{i}.^2));

    % calculate the change in angle in ML plane from the tip position at
    % contact of previous lick to current lick's protrusion
    for m = 1:(size(tip_xf_cont_filt{i}, 1) - 1)
        prot_angle_change{i}(m, :) = asind((tip_xf_prot_filt{i}(m + 1, :) - tip_xf_cont_filt{i}(m, :))./sqrt((tip_xf_prot_filt{i}(m + 1, :).^2) + (tip_yf_prot_filt{i}(m + 1, :).^2)));
        lat_disp_change{i}(m, :) = tip_xf_prot_filt{i}(m + 1, :) - tip_xf_cont_filt{i}(m, :);
    end
    
    % calculate the 'error angle' in the ML plane - defined as the angle
    % between ML tip position at contact relative to ML location of contact on
    % tongue surface
    ml_error{i} = contact_centroid_ml - tip_xf_cont_filt{i};
    ml_error_ang{i} = asind((contact_centroid_ml - tip_xf_cont_filt{i})./sqrt(contact_centroid_ml.^2 + contact_centroid_ap.^2));

    x_temp = tip_xf_cont_filt{i};
    y_temp = tip_yf_cont_filt{i};
    cont_x_temp = contact_centroid_ml;
    cont_y_temp = contact_centroid_ap;
    contact_ang{i} = atan2d(cont_x_temp.*y_temp - cont_y_temp.*x_temp, cont_x_temp.*x_temp + cont_y_temp.*y_temp);   
    contact_ang_hc{i} = asind((contact_centroid_ml)./sqrt(contact_centroid_ml.^2 + contact_centroid_ap.^2));

    % calculate the magnitude of the vector pointing from the tongue tip at
    % contact the the contact location on the tongue surface
    tip_contact_vect_mag{i} = sqrt(((tip_xf_cont_filt{i} - contact_centroid_ml).^2) + ((tip_yf_cont_filt{i} - contact_centroid_ap).^2));
    
    % find the midline of the tongue - defined as the line from tongue centroid
    % to tongue tip at contact. calculate the perpendicular distance (in the ML
    % plane) of the spout contact point on the tongue surface from this line.
    tip_contact_vect_dist{i} = abs((tip_yf_cont_filt{i} - centroid_yf_cont_filt{i}).*(centroid_xf_cont_filt{i} - contact_centroid_ml) - (centroid_yf_cont_filt{i} - contact_centroid_ap).*(tip_xf_cont_filt{i} - centroid_xf_cont_filt{i}))./sqrt(((tip_yf_cont_filt{i} - centroid_yf_cont_filt{i}).^2) + ((tip_xf_cont_filt{i} - centroid_xf_cont_filt{i}).^2));

    % calculate the ML angle from tip at contact to the contact location on the
    % tongue surface
    tip_contact_vect_ang{i} = acosd((tip_yf_cont_filt{i} - contact_centroid_ap)./tip_contact_vect_mag{i});

    % find the midline of the tongue - defined as the line from tongue centroid
    % to tongue tip at contact. if the contact happened to the left of this
    % line, then sign it with a negative - if the contact happened to the right
    % of this line, sign it with a positive. 
    tip_contact_vect_sign{i} = ((tip_yf_cont_filt{i} - centroid_yf_cont_filt{i}).*(contact_centroid_ml - centroid_xf_cont_filt{i}) - (tip_xf_cont_filt{i} - centroid_xf_cont_filt{i}).*(contact_centroid_ap - centroid_yf_cont_filt{i}));
    
    % apply this to tip_contact_vect_mag and _dist to get signed magnitudes and
    % distances
    for m = 1:size(tip_contact_vect_sign{i}, 1)
        for n = 1:size(tip_contact_vect_sign{i}, 2)
            if tip_contact_vect_sign{i}(m, n) < 0
                tip_contact_mag{i}(m, n) = -1*tip_contact_vect_mag{i}(m, n);
                tip_contact_dist{i}(m, n) = -1*tip_contact_vect_dist{i}(m, n);
                tip_contact_ang{i}(m, n) = -1*tip_contact_vect_ang{i}(m, n);
            elseif tip_contact_vect_sign{i}(m, n) == 0
                tip_contact_mag{i}(m, n) = tip_contact_vect_mag{i}(m, n);
                tip_contact_dist{i}(m, n) = 0;
                tip_contact_ang{i}(m, n) = 0;
            elseif tip_contact_vect_sign{i}(m, n) > 0
                tip_contact_mag{i}(m, n) = tip_contact_vect_mag{i}(m, n); 
                tip_contact_dist{i}(m, n) = tip_contact_vect_dist{i}(m, n);
                tip_contact_ang{i}(m, n) = tip_contact_vect_ang{i}(m, n);
            end
        end
    end
end

contact_vars_struct = struct('prot_angle', prot_angle);
[contact_vars_struct.prot_angle_change] = deal(prot_angle_change{:});
[contact_vars_struct.lat_disp_change] = deal(lat_disp_change{:});
[contact_vars_struct.ml_error_ang] = deal(ml_error_ang{:});
[contact_vars_struct.contact_ang] = deal(contact_ang{:});
[contact_vars_struct.contact_ang_hc] = deal(contact_ang_hc{:});
[contact_vars_struct.ml_error] = deal(ml_error{:});
[contact_vars_struct.tip_contact_mag] = deal(tip_contact_mag{:});
[contact_vars_struct.tip_contact_dist] = deal(tip_contact_dist{:});
[contact_vars_struct.tip_contact_ang] = deal(tip_contact_ang{:});
[contact_vars_struct.contact_centroid_filt_mask] = deal(contact_centroid_filt_mask{:});
