function [vel_profs] = create_snapshot_speed_profiles(vel_dir,inland_idx,inland_vel,ref_adjust,vdtmin,vdtmax,zdatetime)
%%% Extract seasonal speed profiles from transect of ITS_LIVE velocities

%INPUTS:
%vel_dir = velocity data directory
%inland_idx = index or the point along the profile where the referencing begins
%inland_vel = starting reference for velocity pts (code set-up that it grabs data for all points less than this reference... points move inland)
%ref_adjust = determines how much you are shifting the points in the seasonal profiles when converting to flow-following reference frame
%vterm_idx = annual terminus reference point
%vdtmin = minimum time separation for image pairs used to calculate velocity (days)
%vdtmax = maximum time separation for image pairs used to calculate velocity (days)
%zdatetime = berg_dates; %datetimes of DEMs

%OUTPUTS:
%vel_profs = average profiles of speed data with date ranges overlapping each DEM (m/yr)

%loop through the pts
vel_pts = dir(vel_dir);
vel_profs = NaN(max(inland_idx),length(zdatetime));
for i = 1:length(vel_pts)
    if contains(vel_pts(i).name,'velocity')
        pt_ref = str2num(vel_pts(i).name(end-5:end-4));

        if pt_ref <= inland_vel %setting this <= and making rel_ref = inland_idx - pt_ref+1 gives the inland pt on the glacier

            %read the file
            V = readtable([vel_dir,vel_pts(i).name]);

            %filter out all the velocities based on temporal resolution
            short_dts = find(V.days_dt>vdtmin & V.days_dt<vdtmax); %get rid of all velocities with coarse temporal resolution
            vel_dates = V.mid_date(short_dts); vel_dts = V.days_dt(short_dts);
            vel_dateranges = [vel_dates - 0.5*vel_dts, vel_dates + 0.5*vel_dts];
            % vmos = month(vel_dates); vyrs = year(vel_dates);
            vels = V.velocity_m_yr_(short_dts); vels(vels == 0) = NaN;

            %calculate the average of all speeds at a point that overlap each DEM date
            for p = 1:length(zdatetime)
                %find velocity image pairs that overlap the DEM date
                in = isbetween(zdatetime(p),vel_dateranges(:,1),vel_dateranges(:,2));
                if sum(in) > 0

                    %isolate the velocities
                    v_temp = vels(in);

                    %use the DEMs from that season to come
                    %up with the position wrt the terminus
                    rel_ref = inland_idx(p)-pt_ref+ref_adjust;
                    if rel_ref >=1
                        vel_profs(rel_ref,p) = nanmean(v_temp);
                    end

                    clear v_temp rel_ref;
                end
                clear in;
            end

            clear V short_dts vel_dates vel_dts vmos vyrs vels vel_dateranges;
        end
        clear pt_ref;
    end
end
vel_profs = vel_profs'; %rotate to match H_profs

%display DEM dates without speed observations
for p = 1:length(zdatetime)
    if sum(~isnan(vel_profs(p,:))) == 0
        disp(['No speeds for ',char(zdatetime(p))]);
    end
end


end