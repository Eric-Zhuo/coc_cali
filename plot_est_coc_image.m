% Image Plot of Central Omnidirectional Camera

function [] = plot_est_coc_image(world_point_image, world_point, mirror_contour_point, K, l, Rotation, Translation, k, d_0, Color)

    [~, ~, Image_num] = size(world_point_image);
    
    u0 = K(1, 3);
    v0 = K(2, 3);



    for n = 1 : Image_num
        
        h = figure;

        hold on

        mirror_handle = plot(mirror_contour_point(1, :), mirror_contour_point(2, :), '-', 'Color', [0.5;0.5;0.5], 'LineWidth', 1); 

        point_principal_handle = plot(u0, v0, 'ro', 'MarkerSize', 8, 'LineWidth', 1);

        color_num = Color(mod(n-1, size(Color,1)) + 1, :);

        [est_world_point_image, ~] = generate_cos_world_point_image(world_point, 1, K, l, Rotation(:, :, n), Translation(:, n), k, []);

        chessboard_corner_handle = plot(world_point_image(1, :, n), world_point_image(2, :, n), 'o', 'Color', color_num, 'MarkerSize', 6); 

        est_image_handle = plot(est_world_point_image(1, :), est_world_point_image(2, :), 'g+');

        hold off


        all_x = mirror_contour_point(1,:); all_y = mirror_contour_point(2,:);

        if ~isempty(all_x)
            dist = sqrt((all_x - u0).^2 + (all_y - v0).^2);
            maxd = max(dist);
            margin = 0.1 * maxd;
            xlim([u0 - maxd - margin, u0 + maxd + margin]);
            ylim([v0 - maxd - margin, v0 + maxd + margin]);
        end

        t = title(['Estimated coc image ', num2str(n)]);
        set(t, 'FontName', 'Arial');

        le1 = legend([mirror_handle, point_principal_handle, chessboard_corner_handle, est_image_handle], ...
            {'Mirror contour', 'Principal point', 'Chessboard corner','Estimated image point'}, ...
            'Location', 'northeast', 'Orientation', 'horizontal');
        le1.NumColumns = 2;
        set(le1, 'FontName', 'Arial', 'FontSize', 14);

        set(gca,'xcolor','k','ycolor','k','linewidth',0.5,'FontSize',18)

        A1x = xlabel('$u$ in pixels');
        set(A1x,'Interpreter','latex','FontName','Times New Roman','FontSize',24)
        A1y = ylabel('$v$ in pixels');
        set(A1y,'Interpreter','latex','FontName','Times New Roman','FontSize',24)

        box on
        % grid on
        set(gca, 'GridAlpha', 0.05);

        set(h, 'units', 'centimeters', 'PaperUnits', 'centimeters');
        set(h, 'PaperPosition', [0, 0, 25, 25]);
        set(h, 'PaperSize', [25, 25]);
        print(h, ['Estimated coc image ', num2str(n)], '-dpng', '-r300');

    end

end
