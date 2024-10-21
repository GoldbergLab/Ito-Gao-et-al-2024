function [selectivity_bin, selectivity_rand_bin] = calculate_percent_site_selectivity(selectivity_site, num_shuffles)

hemi = nan(1, size(selectivity_site, 1));
for i = 1:size(selectivity_site, 1)
    if contains(selectivity_site{i, 3}, 'right')
        hemi(1, i) = 3;
    elseif contains(selectivity_site{i, 3}, 'left')
        hemi(1, i) = 1;
    end
end

center_sites = [selectivity_site{[selectivity_site{:, 1}] == 2, 4}];
contra_sites = [selectivity_site{([selectivity_site{:, 1}] == 1 & hemi == 3)| ([selectivity_site{:, 1}] == 3 & hemi == 1), 4}];
ipsi_sites = [selectivity_site{([selectivity_site{:, 1}] == 1 & hemi == 1) | ([selectivity_site{:, 1}] == 3 & hemi == 3), 4}];
center_contra_sites = [selectivity_site{([selectivity_site{:, 1}] == 4 & hemi == 3) | ([selectivity_site{:, 1}] == 5 & hemi == 1), 4}];
center_ipsi_sites = [selectivity_site{([selectivity_site{:, 1}] == 4 & hemi == 1)| ([selectivity_site{:, 1}] == 5 & hemi == 3), 4}];
site_bin_size = 0.2;
site_bin = 2.9:site_bin_size:3.9;

center_bin = nan(numel(site_bin) - 1, 1);
contra_bin = nan(numel(site_bin) - 1, 1);
ipsi_bin = nan(numel(site_bin) - 1, 1);
center_contra_bin = nan(numel(site_bin) - 1, 1);
center_ipsi_bin = nan(numel(site_bin) - 1, 1);
tot = nan(numel(site_bin) - 1, 1);
for i = 1:numel(site_bin) - 1
    center_bin(i, 1) = sum(center_sites >= site_bin(i) & center_sites < site_bin(i + 1));
    contra_bin(i, 1) = sum(contra_sites >= site_bin(i) & contra_sites < site_bin(i + 1));
    ipsi_bin(i, 1) = sum(ipsi_sites >= site_bin(i) & ipsi_sites < site_bin(i + 1));
    center_contra_bin(i, 1) = sum(center_contra_sites >= site_bin(i) & center_contra_sites < site_bin(i + 1));
    center_ipsi_bin(i, 1) = sum(center_ipsi_sites >= site_bin(i) & center_ipsi_sites < site_bin(i + 1));

    tot(i, 1) = sum(center_bin(i, 1) + contra_bin(i, 1) + ipsi_bin(i, 1) + center_contra_bin(i, 1) + center_ipsi_bin(i, 1));

    center_bin(i, 1) = center_bin(i, 1)/tot(i, 1);
    contra_bin(i, 1) = contra_bin(i, 1)/tot(i, 1);
    ipsi_bin(i, 1) = ipsi_bin(i, 1)/tot(i, 1);
    center_contra_bin(i, 1) = center_contra_bin(i, 1)/tot(i, 1);
    center_ipsi_bin(i, 1) = center_ipsi_bin(i, 1)/tot(i, 1);

end

center_bin_rand = nan(numel(site_bin) - 1, num_shuffles);
contra_bin_rand = nan(numel(site_bin) - 1, num_shuffles);
ipsi_bin_rand = nan(numel(site_bin) - 1, num_shuffles);
center_contra_bin_rand = nan(numel(site_bin) - 1, num_shuffles);
center_ipsi_bin_rand = nan(numel(site_bin) - 1, num_shuffles);
for j = 1:num_shuffles
    site_rand_ind = randi(size(selectivity_site, 1), 1, size(selectivity_site, 1));
    selectivity_rand = selectivity_site(site_rand_ind, :);

    hemi_rand = nan(1, size(selectivity_rand, 1));
    for i = 1:size(selectivity_rand, 1)
        if contains(selectivity_rand{i, 3}, 'right')
            hemi_rand(1, i) = 3;
        elseif contains(selectivity_rand{i, 3}, 'left')
            hemi_rand(1, i) = 1;
        end
    end

    center_sites_rand = [selectivity_rand{[selectivity_rand{:, 1}] == 2, 4}];
    contra_sites_rand = [selectivity_rand{([selectivity_rand{:, 1}] == 1 & hemi_rand == 3)| ([selectivity_rand{:, 1}] == 3 & hemi_rand == 1), 4}];
    ipsi_sites_rand = [selectivity_rand{([selectivity_rand{:, 1}] == 1 & hemi_rand == 1) | ([selectivity_rand{:, 1}] == 3 & hemi_rand == 3), 4}];
    center_contra_sites_rand = [selectivity_rand{([selectivity_rand{:, 1}] == 4 & hemi_rand == 3) | ([selectivity_rand{:, 1}] == 5 & hemi_rand == 1), 4}];
    center_ipsi_sites_rand = [selectivity_rand{([selectivity_rand{:, 1}] == 4 & hemi_rand == 1)| ([selectivity_rand{:, 1}] == 5 & hemi_rand == 3), 4}];
    for i = 1:numel(site_bin) - 1

        center_bin_rand(i, j) = sum(center_sites_rand >= site_bin(i) & center_sites_rand < site_bin(i + 1));
        contra_bin_rand(i, j) = sum(contra_sites_rand >= site_bin(i) & contra_sites_rand < site_bin(i + 1));
        ipsi_bin_rand(i, j) = sum(ipsi_sites_rand >= site_bin(i) & ipsi_sites_rand < site_bin(i + 1));
        center_contra_bin_rand(i, j) = sum(center_contra_sites_rand >= site_bin(i) & center_contra_sites_rand < site_bin(i + 1));
        center_ipsi_bin_rand(i, j) = sum(center_ipsi_sites_rand >= site_bin(i) & center_ipsi_sites_rand < site_bin(i + 1));

        tot_rand = sum(center_bin_rand(i, j) + contra_bin_rand(i, j) + ipsi_bin_rand(i, j) + center_contra_bin_rand(i, j) + center_ipsi_bin_rand(i, j));

        center_bin_rand(i, j) = center_bin_rand(i, j)/tot_rand;
        contra_bin_rand(i, j) = contra_bin_rand(i, j)/tot_rand;
        ipsi_bin_rand(i, j) = ipsi_bin_rand(i, j)/tot_rand;
        center_contra_bin_rand(i, j) = center_contra_bin_rand(i, j)/tot_rand;
        center_ipsi_bin_rand(i, j) = center_ipsi_bin_rand(i, j)/tot_rand;

    end
end

selectivity_bin = {center_bin; contra_bin; ipsi_bin; center_contra_bin; center_ipsi_bin};
selectivity_rand_bin = {center_bin_rand; contra_bin_rand; ipsi_bin_rand; center_contra_bin_rand; center_ipsi_bin_rand};