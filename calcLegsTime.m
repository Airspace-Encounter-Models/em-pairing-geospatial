function [course, d_nm,t_legs_s] = calcLegsTime(lat,lon,v_kts)
% Calculate distance between legs
[course,d_nm] = legs(lat,lon,'rh');

d_nm = abs(d_nm);

% Add functionality to account for v=0 and require altitude

% Calculate time for each leg
if numel(v_kts) == 1
    % Constant airspeed
    t_legs_hr = d_nm / abs(v_kts);
else
    t_legs_hr = d_nm ./ abs(v_kts(1:end-1));
end

% Account for when v_kts = 0
t_legs_hr(isinf(t_legs_hr)) = 0;

t_legs_s = t_legs_hr * 3600;

% Ensure positive time
% We need to do this because helicopters can fly backwards
t_legs_s = abs(t_legs_s);
