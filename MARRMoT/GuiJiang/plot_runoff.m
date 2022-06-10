function plot_runoff(t, Q_obs, Q_sim, warmup, cal_idx, eval_idx)
    figure('color', 'w');
    box on;
    hold all;
    
    ymax = ceil(max(max(Q_obs), max(Q_sim)));
    % ymin = 0;
    
    % Flows
    l(1) = plot(t, Q_obs, 'k');
    l(2) = plot(t, Q_sim, 'r');

    % Dividing line
    l(3) = plot([t(max(cal_idx)), t(max(cal_idx))], [0, ymax], '--b', 'linewidth', 2);
    l(4) = plot([t(warmup), t(warmup)], [0, ymax], '--g', 'linewidth', 2);

    % Legend & text
    l = legend(l, 'Q_{obs}', 'Q_{sim}', 'Cal // Eval', 'warmup // Cal', 'Location', 'northwest');
    title('Model calibration and evaluation results')
    ylabel('Streamflow [mm/hr]')
    xlabel('Time [ymd_h]')
    
    of_name = 'of_KGE';
    of_cal = feval(of_name, Q_obs, Q_sim, cal_idx);
    of_eval = feval(of_name, Q_obs, Q_sim, eval_idx);
    txt_cal = sprintf('Calibration period \nKGE = %.2f ', of_cal);
    txt_eval = sprintf('Evaluation period \nKGE = %.2f', of_eval);
    text(t(round(mean(cal_idx))), ymax*0.8, txt_cal, 'fontsize', 16, 'HorizontalAlignment', 'left');
    text(t(round(mean(eval_idx))), ymax*0.8, txt_eval, 'fontsize', 16, 'HorizontalAlignment', 'left');
    set(gca, 'fontsize', 16);

    % Other settings
    datetick;
    ylim([0, ymax])
    set(gca, 'TickLength', [0.005, 0.005])
end
