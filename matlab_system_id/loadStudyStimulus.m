function stimulus = loadStudyStimulus(condition)
% Inspect the new acquisition profiles pulled from the other PC.
arguments
    condition (1,1) string {mustBeMember(condition, ...
        ["study_sweep90","study_multicosine10"])} = "study_sweep90"
end
root = fileparts(fileparts(mfilename('fullpath')));
p = jsondecode(fileread(fullfile(root,'research_sync','study_profiles',condition+".json")));
if condition == "study_sweep90"
    xva = p.samples_xva;
    t = (0:size(xva,1)-1)'*p.sample_step_s;
else
    % Read the exported design samples including the startup/shutdown envelope.
    c = readtable(fullfile(root,'2611','261124','stimulus_design','multicosine10.csv'));
    t = c.time_s;
    xva = [c.ap_displacement_m c.ap_velocity_m_s c.ap_acceleration_m_s2];
end
stimulus = struct('Profile',p,'Time',t,'Position',timeseries(xva(:,1),t), ...
    'Velocity',timeseries(xva(:,2),t),'Acceleration',timeseries(xva(:,3),t));
figure('Name',condition);
tiledlayout(3,1);
labels = {'Visual position [m]','Visual velocity [m/s]','Visual acceleration [m/s^2]'};
for k=1:3
    nexttile; plot(t,xva(:,k)); grid on; ylabel(labels{k});
end
xlabel('Time [s]');
end
