% Grasses on a Sandy Beach: Beneficial Relationship with Burial
% A Model of Obstacle Dunes from Disrupted Wind Field
% Suspension and Diffusion
% Calculating Plant Area Covered
% written by Emily Fairfax
% April 23rd, 2016


%% Introduction to the Model

% *** Motivation *** 
% The motivation for this model is to demonstrate how a plant's
% physiological response to burial by sand affects both it's growth pattern
% through time as well as the morphology of obstacle dunes which form
% downwind of the plant. It takes into account sand falling out of
% suspension behind the plant in the plant's "windshadow," as well as
% hillslope diffusion of the sand as it falls out of suspension and is
% added to the bed.

% *** Relationship to Burial in Sand ***
% In this particular version of the model, the plant has a beneficial
% relationship with burial, such as that found in Sea Oats and American
% Dunegrass. This means that the growth rate of the plant is based on the
% following equation: G = c*gB(1-B), where G is the growth rate in vertical
% height/time, g is a scaling factor of a base growth rate, and B is the
% scaled biomass of the plant - in this case B is represented by the
% plant's height as a percentage of a max possible growth height. The
% "growth factor" term, c, accounts for how much of the plant is buried in
% the sand. There is an optimal burial percentage defined in the model -
% for burials less than or equal to this percentage, the amount of benefit
% the plant receives from burial increases with increasing burial. Above
% this percentage and increasing burial begins to decrease the benefit the
% plant receives from being buried. In this model the plant is never buried
% enough to cause negative growth.

% *** Model Output ***
% This model will output a movie of the dune forming as it runs, and then
% upon completion produces figures showing the following values through
% time:total accumulated sand in the dune, plant height, percentage of
% plant buried in the sand, growth rate and aboveground exposed plant
% biomass. It also produces a plot of the growth factor term as it related
% to the buried biomass percentage. 

% *** Notes on Running the Model ***
% Since this model includes diffusion, many of the parameters are fairly
% sensitive - in particular the time step and spacial steps, dx and dt.
% Changing these will likely result in an unstable model.

% *** Enjoy!*** 
% Feel free to email me with any questions at emily.fairfax@gmail.com

clear all

%% Initialize
% Constants
    % Topography Constants
    beach_length = 15; % length of the beach in meters
    beach_slope = .003; % slope of beach (for linear topography)
    
    % Grass Constants
    max_grass_height = 1; %maximum height grass can grow
    min_growth = 0.015; %grass grows at a minimum of this amount per timestep, m
    grass_height = 0.05; % height of grass tuft in meters, initialize with some starting height
    block_factor = 0.3; % effective height of grass that acts as a wind block
    optimal_coverage_percent = 0.08; %*100 = percent covered for best growth
    
    % Suspended Load Constants
    suspended_sand_factor = .0006; % meters of sediment per timestep that could fall out of suspension if wind speed -> 0 
    
    %Diffusion Constants
    k = 10; %hillslope diffusivity constant in meter^2/timestep
    rho_sand = 1700; %density of soil in kg/m^3

% Space Array
    dx = .1; % horizontal step in meters
    xmin = 0; % min horizontal space in meters
    xmax = beach_length; % max horizontal space in meters
    x = xmin:dx:xmax; % x array
    xedge = xmin+(dx/2):dx:xmax-(dx/2); % xedge array  
    
% Time Array
    dt = .25; % time step in days
    max_days = 365; % number of minutes to run the simulation
    tmax = max_days; % tmax in days
    t = 0:dt:tmax; % time array
    imax = length(t); % loop length
    nplot = max_days; % number of plots
    tplot = tmax/nplot; % times at which plots will be made
    
% Variable Arrays
    N = length(x); % nodes
    
    % Sand Arrays
    sand_elevation = zeros(1,N); %create array for elevation profile of earth
    sand_elevation(1:N) = .1+beach_slope.*x; % linear valley profile
    original_sand_no_grass = sand_elevation; %make a copy of the sand array that will not change within the loop
    
    % Grass Arrays
    % Grass Location and Dimensions
    grass1_start = (2*xmax/10)/dx; %the grass starts at this position in the array
    grass1_end = (2*xmax/10+.5)/dx; %the grass ends at this position in the array
    grass_length = (grass1_end-grass1_start)*dx; %the length of the grass in meters
    grass1_end_effective_blocking_height = sand_elevation(grass1_start)+block_factor*grass_height; % find the elevation of the effective blocking height of the grass
    
    % Grass Elevations
    grass_elevation = sand_elevation; %create array for elevation including grass
    grass_elevation(grass1_start:grass1_end) = sand_elevation(grass1_start)+grass_height; % grass tuft is elevated above sand
    grass1_end_elevation = grass_elevation(grass1_end); % find the elevation of the grass at the end of the tuft
    
    % Wind Shadow Behind Grass
    grass1_shadow_end = (grass1_end+10/dx*grass_height); % grass wind shadow extends ten times the height of the tuft behind the grass
    grass1_shadow_length = length(x(grass1_end:grass1_shadow_end)); %this is the length of the wind shadow in the array
    grass1_shadow_end_elevation = grass_elevation(grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    % Wind Disturbance Arrays
    % Treat Wind Disturbance as a Linear Trajectory from Effective Blocking
    % Height to End of Shadow Sand Elevation, simple y = mx+b formula
    wind_disturbance_height = grass1_shadow_end_elevation-((grass1_shadow_end*dx*(grass1_shadow_end_elevation-grass1_end_effective_blocking_height)/((grass1_shadow_end*dx)-grass1_start*dx)))+((grass1_shadow_end_elevation-grass1_end_effective_blocking_height)/((grass1_shadow_end*dx)-grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    
    %Flux Arrays
    Q_diffusion = zeros(1,N); %preallocate Q array for diffusion
    dQdx = zeros(1,N-2); %preallocate dQdx array
    
%Plotting Arrays
    bottomline=zeros(N:1);%Create a bottom line array for use in filling plots, needs to be same length as N
    bottomline(1:N)=-200000; %arbitrary very negative number
    range_of_percents = zeros(1,imax-1);
    range_of_percents = 0:20/(imax-1):20;
    range_of_zeros = zeros(1,imax-1);
    range_of_zeros_for_time = zeros(1,imax);

%% Run
for i = 2:imax %this is the time loop
    
    % Wind Speed Dependent Suspended Sediments
        % Calculate How Much Sediment to Deposit Along Wind Shadow
        percent_max_windspeed = ones(1,N-1); %initialize windspeed factor so that all wind is max windspeed

        for j = 1:(grass1_shadow_length) %loop through the wind shadow (zone of deposition), this is the only place we are going to alter the wind speed

            if sand_elevation(grass1_end+j)<=wind_disturbance_height(grass1_end+j) %if the sand isn't taller than the grass in a given unit cell
                % Linear Wind Speed Profile
                linear_wind_profile = 1/(grass1_shadow_end-grass1_end); %the slope of the line that goes from the end of the tuft to the end of the wind shadow
                percent_max_windspeed(grass1_end+j) = linear_wind_profile*x(grass1_end+j)-linear_wind_profile*grass1_end; %the line with values from 1->0 over the length of the wind shadow representing how much sand is going to fall out of suspension as a fraction of the total possible amount to fall out of suspension
            end

        end
        % Net Flux of Sand from Suspended Sediments
        dQdx_suspended_sand = suspended_sand_factor*(1 - percent_max_windspeed); %net flux into each box from suspended load being deposited is dependent on the wind speed profile as affected by the grass tuft
    
    % Diffusion of Sand
        %Flux of Sand from Diffusion
        Q_diffusion(1:N-1) = k.*diff(sand_elevation(1:N))./dx; %carry soil out of each dx box according to flux eq. for hillslope diffusion
        
        % Net Flux of Sand from Diffusion
        dQdx_diffusion(1:N-2) = (1/rho_sand)*(diff(Q_diffusion(1:N-1))./dx); %figure out the net flux in each box from sand diffusion, independent of deposited sand
        dQdX_diffusion(1:N-1) = [dQdx_diffusion(1) dQdx_diffusion]; %flux sand out of left hand side of box, ultimate downslope position
    
    % Combined Net Fluxes of Sand into Each Cell: Suspended Sediment + Diffusion
    dQdx(1:N-1) = dQdx_suspended_sand(1:N-1)+dQdX_diffusion(1:N-1);%+dQdx_diffusion(1:N-1); % diff Q to get dQdx
    dQdx(N) = dQdx(N-1); % net flux out of last box is same as box before it, allows model to drain...boundary condition

    % Sand Height Changes
    dHdt = dQdx; % conservation statement: change in sand height is a function of net flux of sand into box, no source or sink terms
    sand_elevation = sand_elevation + dHdt.*dt; % update total height of sand
    Hneg = find(sand_elevation<0); % find negative sand heights, shouldn't be any but just in case
    sand_elevation(Hneg)=0; % set negative sand heights to zero... no such thing as negative sand.
    
    %Calculate Area of Plant that is Covered by Sand using Trapezoidal
    %Integration
    plant_area_total(i) = trapz(x(grass1_start:grass1_end),grass_elevation(grass1_start:grass1_end))-trapz(x(grass1_start:grass1_end),original_sand_no_grass(grass1_start:grass1_end)); % Calculate area of the whole plant
    sand_area_accumulated_total(i) = trapz(x,sand_elevation)-trapz(x,original_sand_no_grass); %calculate area of all sand that has been added to the system since the time loop started
    plant_area_exposed(i) = trapz(x(grass1_start:grass1_end),grass_elevation(grass1_start:grass1_end))-trapz(x(grass1_start:grass1_end),sand_elevation(grass1_start:grass1_end)); %calculate how much of the plant is not buried under the sand
    plant_area_covered_percent(i) = (plant_area_total(i)-plant_area_exposed(i))/plant_area_total(i); %calculate the percent of the plant buried in the sand in decimal form
    
    %Grow the Grass
    if grass_height(i-1)<=max_grass_height %if the grass is shorter than the max height, apply the following growth laws
        
        if plant_area_covered_percent(i)<=optimal_coverage_percent %and if the plant is covered less than or equal to the optimal burial percent, the bonus from being buried increased with additional coverage
            grass_growth(i) = (min_growth*grass_height(i-1)/max_grass_height*(1-grass_height(i-1)/max_grass_height)*(1+plant_area_covered_percent(i)/optimal_coverage_percent))*dt;%growth curve for grass G = c*G(1-B) scaled by covered percentage factor, c
            biomass_growth(i) = (min_growth*grass_height(i-1)/max_grass_height*(1-grass_height(i-1)/max_grass_height))*dt;%growth curve for grass G = G(1-B) only considering the biomass terms, not coverage by sand
            burial_growth_factor(i) = (1+plant_area_covered_percent(i)/optimal_coverage_percent); %the real time burial growth factor affecting the plant is the c term
        
        else %if the plant has gone past the optimal coverage percentage, the bonus from being buried starts to decrease with additional coverage
            grass_growth(i) = (min_growth*grass_height(i-1)/max_grass_height*(1-grass_height(i-1)/max_grass_height)*(1+plant_area_covered_percent(i)/optimal_coverage_percent-2*(plant_area_covered_percent(i)-optimal_coverage_percent)/optimal_coverage_percent))*dt;%growth curve for grass G = c*G(1-B) scaled by covered percentage factor, c
            biomass_growth(i) = (min_growth*grass_height(i-1)/max_grass_height*(1-grass_height(i-1)/max_grass_height))*dt;%growth curve for grass G = G(1-B) only considering the biomass terms, not coverage by sand
            burial_growth_factor(i) = (1+plant_area_covered_percent(i)/optimal_coverage_percent-2*(plant_area_covered_percent(i)-optimal_coverage_percent)/optimal_coverage_percent);%the real time burial growth factor affecting the plant is the c term
        
        end %end the coverage based loop
    
    else %if the grass is trying to grow larger than the max height, stop it
        grass_growth(i) = 0; %no growth past max height
        biomass_growth(i) = 0; %no growth past max height
        burial_growth_factor(i) = 0; %no growth past max height
   
    end %end the grass growth loop
    
    %Update Grass Elevation, Grass Height, and other Grass Parameters with Growth
    grass_elevation(grass1_start:grass1_end) = original_sand_no_grass(grass1_start:grass1_end)+grass_height(i-1)+grass_growth(i); %update grass elevation with growth this timestep
    grass_height(i) = grass_elevation(grass1_start)-original_sand_no_grass(grass1_start); %update grass height with grass growth this time step
    grass1_end_elevation = grass_elevation(grass1_end); % find the elevation of the grass at the end of the tuft
    grass1_end_effective_blocking_height = sand_elevation(grass1_start)+block_factor*grass_height(i); % find the elevation of the effective blocking height of the grass
    grass1_shadow_end = ceil((grass1_end+10/dx*grass_height(i))); % grass wind shadow extends ten times the height of the tuft behind the grass
    grass1_shadow_length = length(x(grass1_end:grass1_shadow_end)); %this is the length of the wind shadow in the array
    grass1_shadow_end_elevation = original_sand_no_grass(grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    %Update Wind Disturbance Height based on New Grass Height
    wind_disturbance_height = grass1_shadow_end_elevation-((grass1_shadow_end*dx*(grass1_shadow_end_elevation-grass1_end_effective_blocking_height)/((grass1_shadow_end*dx)-grass1_start*dx)))+((grass1_shadow_end_elevation-grass1_end_effective_blocking_height)/((grass1_shadow_end*dx)-grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    
    %For plotting percentage covered vs growth factor over a predetermined
    %range of percentages
    if range_of_percents(i)<= optimal_coverage_percent*100 %if the percent covered is less than or equal to the optimal percent
        range_of_growth_factors(i) = (1+0.01*range_of_percents(i)/optimal_coverage_percent); %growth factor is a increasing with additional coverage
    else % if the percent covered is greater than the optimal coverage percent
        range_of_growth_factors(i)=(1+0.01*range_of_percents(i)/optimal_coverage_percent-2*(0.01*range_of_percents(i)-optimal_coverage_percent)/optimal_coverage_percent); %growth factor decreased with addional coverage
    end %end this percentage loop
    
    %Plot the Results Each Time Step   
    if (rem(t(i),tplot)==0) %only plot if it is a multiple of tplot
        figure(1) 
        clf
        
        plot(x,grass_elevation, 'k') %plot the grass elevation
        hold all
        plot(x,sand_elevation,'k') %plot sand elevation on top of it

        %Color the grass green, do this first so sand can cover grass
        uu = [x,x];        % repeat x values
        gg = [grass_elevation,bottomline];   % vector of upper & lower boundaries
        fill(uu,gg,[.643,.686,.416]) %fill the polygon created

        %Color the sand beige, do this second so sand can cover the grass
        ss = [sand_elevation,bottomline];   % vector of upper & lower boundaries
        fill(uu,ss,[.91,.89,.87]) %fill the polygon created

        %Plot formatting
        title('Beach Evolution Through Time');
        xlabel('Distance Along Profile (m)');
        ylabel('Elevation (m)');
        set(gca,'fontsize',16,'fontname','arial')
        ht=text(4/5*xmax,.8,['  ',num2str(floor(round(t(i)/7))), ' weeks  '],'fontsize',18); %print time in animation
        axis([0 xmax 0 1.3])
        
        pause(0.02)
        hold off
    end
end

%% Finalize

figure(2) %These are plots for measureable values that change with time
clf
%Total Sand Accumulated in System
subplot(3,1,1)
plot(t,sand_area_accumulated_total)
xlabel('Time (days)');
ylabel('(m^2/unit width)');
title('Total Accumulated Sand (m^2/unit width)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 2])

%Plant Area Exposed
subplot(3,1,2)
plot(t,plant_area_exposed)
xlabel('Time (days)');
ylabel('(m^2/unit width)');
title('Plant Area Exposed (m^2/unit width)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 .5])

%Grass Height
subplot(3,1,3)
plot(t,grass_height)
xlabel('Time (days)');
ylabel('(m)');
title('Grass Height (m)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 1.1])



figure(3) %these are plots showing aspects of the sand burial growth factor relationship
clf

%Plant Area Percent Covered
subplot(3,1,1)
plot(t,plant_area_covered_percent*100)
xlabel('Time (days)');
ylabel('(%)');
title('Grass Area Buried (%)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 15])

%Plant Growth Curves
subplot(3,1,2)
plot(t,grass_growth*1000,t,biomass_growth*1000,'--k')
hold on
plot(t,range_of_zeros_for_time,'r')
xlabel('Time (days)');
ylabel('(mm/day)');
title('Vertical Grass Growth (mm/day)');
legend('Total Grass Growth per day','Grass Growth without Burial Growth Factor')
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax -.1 2])

%Plant Growth Factors
subplot(3,1,3)
plot(range_of_percents,range_of_growth_factors)
hold all
plot(range_of_percents,range_of_zeros_for_time,'r')
xlabel('Grass Area Buried (% of total grass area)');
ylabel('Burial Growth Factor');
title('Effect of Burial on Grass Growth')
set(gca,'fontsize',16,'fontname','arial')
axis([0 20 -.1 2.1])

