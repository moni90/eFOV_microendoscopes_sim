function [start_point_gmm,start_point_gauss] = opt_starting_point(z_profile)
start_point_gmm = [...
        nanmin(z_profile), 1, -10, -10, 0.5, 10, 10;
        nanmin(z_profile), 1, -10, -30, 0.5, 10, 10;
        nanmin(z_profile), 1, -10, -50, 0.5, 10, 10;
        nanmin(z_profile), 1, -10, -70, 0.5, 10, 10;
        nanmin(z_profile), 1, -10, -90, 0.5, 10, 10;
        nanmin(z_profile), 1, -20, -20, 0.5, 10, 10;
        nanmin(z_profile), 1, -20, -40, 0.5, 10, 10;
        nanmin(z_profile), 1, -20, -60, 0.5, 10, 10;
        nanmin(z_profile), 1, -20, -80, 0.5, 10, 10;
        nanmin(z_profile), 1, -30, -60, 0.5, 10, 10;
        nanmin(z_profile), 1, -30, -80, 0.5, 10, 10;
        ];

    start_point_gauss = [...
        nanmin(z_profile), 1, -5, 1;
        nanmin(z_profile), 1, -10, 1;
        nanmin(z_profile), 1, -15, 1;
        nanmin(z_profile), 1, -20, 1;
        nanmin(z_profile), 1, -25, 1;
        nanmin(z_profile), 1, -30, 1;
        nanmin(z_profile), 1, -35, 1;
        nanmin(z_profile), 1, -40, 1;
        nanmin(z_profile), 1, -45, 1;
        nanmin(z_profile), 1, -50, 1;
        nanmin(z_profile), 1, -55, 1;
        nanmin(z_profile), 1, -60, 1;
        nanmin(z_profile), 1, -65, 1;
        nanmin(z_profile), 1, -70, 1;
        nanmin(z_profile), 1, -75, 1;
        nanmin(z_profile), 1, -80, 1;
        nanmin(z_profile), 1, -85, 1;
        nanmin(z_profile), 1, -90, 1;
        ];
end