% Grasses on a Sandy Beach: Various Relationships with Burial
% A Comparison
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
% In this particular version of the model, the three relationships with
% burial (beneficial, neutral, and harmful) are compared visually against
% one another. For details on each particular relationship, see the
% introduction section for that relationship's specific model.

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

% Initialize all three relationships here, should be same parameters for
% all three, variable arrays will be unique to each relationship
% Constants
    % Topography Constants
    beach_length = 15; % length of the beach in meters
    beach_slope = .003; % slope of beach (for linear topography)
    
    % Grass Constants
    max_grass_height = 1; %maximum height grass can grow
    min_growth = 0.015; %grass grows at a minimum of this amount per timestep, m
    b_grass_height = 0.05; % height of grass tuft in meters, initialize with some starting height
    n_grass_height = 0.05; % height of grass tuft in meters, initialize with some starting height
    h_grass_height = 0.05; % height of grass tuft in meters, initialize with some starting height
    block_factor = 0.3; % effective height of grass that acts as a wind block
    optimal_coverage_percent = 0.08; %*100 = percent covered for best growth
    plant_death_coverage_percent = 0.2; %*100 = percent covered for death of plant
    
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
    b_sand_elevation = zeros(1,N); %create array for elevation profile of earth
    b_sand_elevation(1:N) = .1+beach_slope.*x; % linear valley profile
    original_sand_no_grass = b_sand_elevation; %make a copy of the sand array that will not change within the loop
    
    % Grass Arrays for Beneficial Relationship
    % Grass Location and Dimensions
    b_grass1_start = (2*xmax/10)/dx; %the grass starts at this position in the array
    b_grass1_end = (2*xmax/10+.5)/dx; %the grass ends at this position in the array
    b_grass_length = (b_grass1_end-b_grass1_start)*dx; %the length of the grass in meters
    b_grass1_end_effective_blocking_height = b_sand_elevation(b_grass1_start)+block_factor*b_grass_height; % find the elevation of the effective blocking height of the grass
    
    % Grass Elevations
    b_grass_elevation = b_sand_elevation; %create array for elevation including grass
    b_grass_elevation(b_grass1_start:b_grass1_end) = b_sand_elevation(b_grass1_start)+b_grass_height; % grass tuft is elevated above sand
    b_grass1_end_elevation = b_grass_elevation(b_grass1_end); % find the elevation of the grass at the end of the tuft
    
    % Wind Shadow Behind Grass
    b_grass1_shadow_end = (b_grass1_end+10/dx*b_grass_height); % grass wind shadow extends ten times the height of the tuft behind the grass
    b_grass1_shadow_length = length(x(b_grass1_end:b_grass1_shadow_end)); %this is the length of the wind shadow in the array
    b_grass1_shadow_end_elevation = b_grass_elevation(b_grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    % Wind Disturbance Arrays
    % Treat Wind Disturbance as a Linear Trajectory from Effective Blocking
    % Height to End of Shadow Sand Elevation, simple y = mx+b formula
    b_wind_disturbance_height = b_grass1_shadow_end_elevation-((b_grass1_shadow_end*dx*(b_grass1_shadow_end_elevation-b_grass1_end_effective_blocking_height)/((b_grass1_shadow_end*dx)-b_grass1_start*dx)))+((b_grass1_shadow_end_elevation-b_grass1_end_effective_blocking_height)/((b_grass1_shadow_end*dx)-b_grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    %Flux Arrays
    b_Q_diffusion = zeros(1,N); %preallocate Q array for diffusion
    b_dQdx = zeros(1,N-2); %preallocate dQdx array
    
    % Grass Arrays for Neutral Relationship
    % Sand Arrays
    n_sand_elevation = zeros(1,N); %create array for elevation profile of earth
    n_sand_elevation(1:N) = .1+beach_slope.*x; % linear valley profile
     
    % Grass Location and Dimensions
    n_grass1_start = (2*xmax/10)/dx; %the grass starts at this position in the array
    n_grass1_end = (2*xmax/10+.5)/dx; %the grass ends at this position in the array
    n_grass_length = (n_grass1_end-n_grass1_start)*dx; %the length of the grass in meters
    n_grass1_end_effective_blocking_height = n_sand_elevation(n_grass1_start)+block_factor*n_grass_height; % find the elevation of the effective blocking height of the grass
    
    % Grass Elevations
    n_grass_elevation = n_sand_elevation; %create array for elevation including grass
    n_grass_elevation(n_grass1_start:n_grass1_end) = n_sand_elevation(n_grass1_start)+n_grass_height; % grass tuft is elevated above sand
    n_grass1_end_elevation = n_grass_elevation(n_grass1_end); % find the elevation of the grass at the end of the tuft
    
    % Wind Shadow Behind Grass
    n_grass1_shadow_end = (n_grass1_end+10/dx*n_grass_height); % grass wind shadow extends ten times the height of the tuft behind the grass
    n_grass1_shadow_length = length(x(n_grass1_end:n_grass1_shadow_end)); %this is the length of the wind shadow in the array
    n_grass1_shadow_end_elevation = n_grass_elevation(n_grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    % Wind Disturbance Arrays
    % Treat Wind Disturbance as a Linear Trajectory from Effective Blocking
    % Height to End of Shadow Sand Elevation, simple y = mx+b formula
    n_wind_disturbance_height = n_grass1_shadow_end_elevation-((n_grass1_shadow_end*dx*(n_grass1_shadow_end_elevation-n_grass1_end_effective_blocking_height)/((n_grass1_shadow_end*dx)-n_grass1_start*dx)))+((n_grass1_shadow_end_elevation-n_grass1_end_effective_blocking_height)/((n_grass1_shadow_end*dx)-n_grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
   
    %Flux Arrays
    n_Q_diffusion = zeros(1,N); %preallocate Q array for diffusion
    n_dQdx = zeros(1,N-2); %preallocate dQdx array
    
    % Grass Arrays for Harmful Relationship
    % Sand Arrays
    h_sand_elevation = zeros(1,N); %create array for elevation profile of earth
    h_sand_elevation(1:N) = .1+beach_slope.*x; % linear valley profile
    
    % Grass Location and Dimensions
    h_grass1_start = (2*xmax/10)/dx; %the grass starts at this position in the array
    h_grass1_end = (2*xmax/10+.5)/dx; %the grass ends at this position in the array
    h_grass_length = (h_grass1_end-h_grass1_start)*dx; %the length of the grass in meters
    h_grass1_end_effective_blocking_height = b_sand_elevation(h_grass1_start)+block_factor*h_grass_height; % find the elevation of the effective blocking height of the grass
    
    % Grass Elevations
    h_grass_elevation = h_sand_elevation; %create array for elevation including grass
    h_grass_elevation(h_grass1_start:h_grass1_end) = h_sand_elevation(h_grass1_start)+h_grass_height; % grass tuft is elevated above sand
    h_grass1_end_elevation = h_grass_elevation(h_grass1_end); % find the elevation of the grass at the end of the tuft
    
    % Wind Shadow Behind Grass
    h_grass1_shadow_end = (h_grass1_end+10/dx*h_grass_height); % grass wind shadow extends ten times the height of the tuft behind the grass
    h_grass1_shadow_length = length(x(h_grass1_end:h_grass1_shadow_end)); %this is the length of the wind shadow in the array
    h_grass1_shadow_end_elevation = h_grass_elevation(h_grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    % Wind Disturbance Arrays
    % Treat Wind Disturbance as a Linear Trajectory from Effective Blocking
    % Height to End of Shadow Sand Elevation, simple y = mx+b formula
    h_wind_disturbance_height = h_grass1_shadow_end_elevation-((h_grass1_shadow_end*dx*(h_grass1_shadow_end_elevation-h_grass1_end_effective_blocking_height)/((h_grass1_shadow_end*dx)-h_grass1_start*dx)))+((h_grass1_shadow_end_elevation-h_grass1_end_effective_blocking_height)/((h_grass1_shadow_end*dx)-h_grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    %Flux Arrays
    h_Q_diffusion = zeros(1,N); %preallocate Q array for diffusion
    h_dQdx = zeros(1,N-2); %preallocate dQdx array
    
%Plotting Arrays
    bottomline=zeros(N:1);%Create a bottom line array for use in filling plots, needs to be same length as N
    bottomline(1:N)=-200000; %arbitrary very negative number
    range_of_percents = zeros(1,imax-1);
    range_of_percents = 0:20/(imax-1):20;
    range_of_zeros = zeros(1,imax-1);
    range_of_zeros_for_time = zeros(1,imax);
    
%% Run
for i = 2:imax %this is the time loop
    
    % Beneficial Relationship
    % Wind Speed Dependent Suspended Sediments
        % Calculate How Much Sediment to Deposit Along Wind Shadow
        b_percent_max_windspeed = ones(1,N-1); %initialize windspeed factor so that all wind is max windspeed

        for j = 1:(b_grass1_shadow_length) %loop through the wind shadow (zone of deposition), this is the only place we are going to alter the wind speed

            if b_sand_elevation(b_grass1_end+j)<=b_wind_disturbance_height(b_grass1_end+j) %if the sand isn't taller than the grass in a given unit cell
                % Linear Wind Speed Profile
                b_linear_wind_profile = 1/(b_grass1_shadow_end-b_grass1_end); %the slope of the line that goes from the end of the tuft to the end of the wind shadow
                b_percent_max_windspeed(b_grass1_end+j) = b_linear_wind_profile*x(b_grass1_end+j)-b_linear_wind_profile*b_grass1_end; %the line with values from 1->0 over the length of the wind shadow representing how much sand is going to fall out of suspension as a fraction of the total possible amount to fall out of suspension
            end

        end
        % Net Flux of Sand from Suspended Sediments
        b_dQdx_suspended_sand = suspended_sand_factor*(1 - b_percent_max_windspeed); %net flux into each box from suspended load being deposited is dependent on the wind speed profile as affected by the grass tuft
    
    % Diffusion of Sand
        %Flux of Sand from Diffusion
        b_Q_diffusion(1:N-1) = k.*diff(b_sand_elevation(1:N))./dx; %carry soil out of each dx box according to flux eq. for hillslope diffusion
        
        % Net Flux of Sand from Diffusion
        b_dQdx_diffusion(1:N-2) = (1/rho_sand)*(diff(b_Q_diffusion(1:N-1))./dx); %figure out the net flux in each box from sand diffusion, independent of deposited sand
        b_dQdX_diffusion(1:N-1) = [b_dQdx_diffusion(1) b_dQdx_diffusion]; %flux sand out of left hand side of box, ultimate downslope position
    
    % Combined Net Fluxes of Sand into Each Cell: Suspended Sediment + Diffusion
    b_dQdx(1:N-1) = b_dQdx_suspended_sand(1:N-1)+b_dQdX_diffusion(1:N-1);%+dQdx_diffusion(1:N-1); % diff Q to get dQdx
    b_dQdx(N) = b_dQdx(N-1); % net flux out of last box is same as box before it, allows model to drain...boundary condition

    % Sand Height Changes
    b_dHdt = b_dQdx; % conservation statement: change in sand height is a function of net flux of sand into box, no source or sink terms
    b_sand_elevation = b_sand_elevation + b_dHdt.*dt; % update total height of sand
    b_Hneg = find(b_sand_elevation<0); % find negative sand heights, shouldn't be any but just in case
    b_sand_elevation(b_Hneg)=0; % set negative sand heights to zero... no such thing as negative sand.
    
    %Calculate Area of Plant that is Covered by Sand using Trapezoidal
    %Integration
    b_plant_area_total(i) = trapz(x(b_grass1_start:b_grass1_end),b_grass_elevation(b_grass1_start:b_grass1_end))-trapz(x(b_grass1_start:b_grass1_end),original_sand_no_grass(b_grass1_start:b_grass1_end)); % Calculate area of the whole plant
    b_sand_area_accumulated_total(i) = trapz(x,b_sand_elevation)-trapz(x,original_sand_no_grass); %calculate area of all sand that has been added to the system since the time loop started
    b_plant_area_exposed(i) = trapz(x(b_grass1_start:b_grass1_end),b_grass_elevation(b_grass1_start:b_grass1_end))-trapz(x(b_grass1_start:b_grass1_end),b_sand_elevation(b_grass1_start:b_grass1_end)); %calculate how much of the plant is not buried under the sand
    b_plant_area_covered_percent(i) = (b_plant_area_total(i)-b_plant_area_exposed(i))/b_plant_area_total(i); %calculate the percent of the plant buried in the sand in decimal form
    
    %Grow the Grass
    if b_grass_height(i-1)<=max_grass_height %if the grass is shorter than the max height, apply the following growth laws
        
        if b_plant_area_covered_percent(i)<=optimal_coverage_percent %and if the plant is covered less than or equal to the optimal burial percent, the bonus from being buried increased with additional coverage
            b_grass_growth(i) = (min_growth*b_grass_height(i-1)/max_grass_height*(1-b_grass_height(i-1)/max_grass_height)*(1+b_plant_area_covered_percent(i)/optimal_coverage_percent))*dt;%growth curve for grass G = c*G(1-B) scaled by covered percentage factor, c
            b_biomass_growth(i) = (min_growth*b_grass_height(i-1)/max_grass_height*(1-b_grass_height(i-1)/max_grass_height))*dt;%growth curve for grass G = G(1-B) only considering the biomass terms, not coverage by sand
            b_burial_growth_factor(i) = (1+b_plant_area_covered_percent(i)/optimal_coverage_percent); %the real time burial growth factor affecting the plant is the c term
        
        else %if the plant has gone past the optimal coverage percentage, the bonus from being buried starts to decrease with additional coverage
            b_grass_growth(i) = (min_growth*b_grass_height(i-1)/max_grass_height*(1-b_grass_height(i-1)/max_grass_height)*(1+b_plant_area_covered_percent(i)/optimal_coverage_percent-2*(b_plant_area_covered_percent(i)-optimal_coverage_percent)/optimal_coverage_percent))*dt;%growth curve for grass G = c*G(1-B) scaled by covered percentage factor, c
            b_biomass_growth(i) = (min_growth*b_grass_height(i-1)/max_grass_height*(1-b_grass_height(i-1)/max_grass_height))*dt;%growth curve for grass G = G(1-B) only considering the biomass terms, not coverage by sand
            b_burial_growth_factor(i) = (1+b_plant_area_covered_percent(i)/optimal_coverage_percent-2*(b_plant_area_covered_percent(i)-optimal_coverage_percent)/optimal_coverage_percent);%the real time burial growth factor affecting the plant is the c term
        
        end %end the coverage based loop
    
    else %if the grass is trying to grow larger than the max height, stop it
        b_grass_growth(i) = 0; %no growth past max height
        b_biomass_growth(i) = 0; %no growth past max height
        b_burial_growth_factor(i) = 0; %no growth past max height
   
    end %end the grass growth loop
    
    %Update Grass Elevation, Grass Height, and other Grass Parameters with Growth
    b_grass_elevation(b_grass1_start:b_grass1_end) = original_sand_no_grass(b_grass1_start:b_grass1_end)+b_grass_height(i-1)+b_grass_growth(i); %update grass elevation with growth this timestep
    b_grass_height(i) = b_grass_elevation(b_grass1_start)-original_sand_no_grass(b_grass1_start); %update grass height with grass growth this time step
    b_grass1_end_elevation = b_grass_elevation(b_grass1_end); % find the elevation of the grass at the end of the tuft
    b_grass1_end_effective_blocking_height = b_sand_elevation(b_grass1_start)+block_factor*b_grass_height(i); % find the elevation of the effective blocking height of the grass
    b_grass1_shadow_end = ceil((b_grass1_end+10/dx*b_grass_height(i))); % grass wind shadow extends ten times the height of the tuft behind the grass
    b_grass1_shadow_length = length(x(b_grass1_end:b_grass1_shadow_end)); %this is the length of the wind shadow in the array
    b_grass1_shadow_end_elevation = original_sand_no_grass(b_grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    %Update Wind Disturbance Height based on New Grass Height
    b_wind_disturbance_height = b_grass1_shadow_end_elevation-((b_grass1_shadow_end*dx*(b_grass1_shadow_end_elevation-b_grass1_end_effective_blocking_height)/((b_grass1_shadow_end*dx)-b_grass1_start*dx)))+((b_grass1_shadow_end_elevation-b_grass1_end_effective_blocking_height)/((b_grass1_shadow_end*dx)-b_grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    
    %For plotting percentage covered vs growth factor over a predetermined
    %range of percentages
    if range_of_percents(i)<= optimal_coverage_percent*100 %if the percent covered is less than or equal to the optimal percent
        b_range_of_growth_factors(i) = (1+0.01*range_of_percents(i)/optimal_coverage_percent); %growth factor is a increasing with additional coverage
    else % if the percent covered is greater than the optimal coverage percent
        b_range_of_growth_factors(i)=(1+0.01*range_of_percents(i)/optimal_coverage_percent-2*(0.01*range_of_percents(i)-optimal_coverage_percent)/optimal_coverage_percent); %growth factor decreased with addional coverage
    end %end this percentage loop
    
    % Neutral Relationship
    % Wind Speed Dependent Suspended Sediments
        % Calculate How Much Sediment to Deposit Along Wind Shadow
        n_percent_max_windspeed = ones(1,N-1); %initialize windspeed factor so that all wind is max windspeed

        for j = 1:(n_grass1_shadow_length) %loop through the wind shadow (zone of deposition), this is the only place we are going to alter the wind speed

            if n_sand_elevation(n_grass1_end+j)<=n_wind_disturbance_height(n_grass1_end+j) %if the sand isn't taller than the grass in a given unit cell
                % Linear Wind Speed Profile
                n_linear_wind_profile = 1/(n_grass1_shadow_end-n_grass1_end); %the slope of the line that goes from the end of the tuft to the end of the wind shadow
                n_percent_max_windspeed(n_grass1_end+j) = n_linear_wind_profile*x(n_grass1_end+j)-n_linear_wind_profile*n_grass1_end; %the line with values from 1->0 over the length of the wind shadow representing how much sand is going to fall out of suspension as a fraction of the total possible amount to fall out of suspension
            end

        end
        % Net Flux of Sand from Suspended Sediments
        n_dQdx_suspended_sand = suspended_sand_factor*(1 - n_percent_max_windspeed); %net flux into each box from suspended load being deposited is dependent on the wind speed profile as affected by the grass tuft
    
    % Diffusion of Sand
        %Flux of Sand from Diffusion
        n_Q_diffusion(1:N-1) = k.*diff(n_sand_elevation(1:N))./dx; %carry soil out of each dx box according to flux eq. for hillslope diffusion
        
        % Net Flux of Sand from Diffusion
        n_dQdx_diffusion(1:N-2) = (1/rho_sand)*(diff(n_Q_diffusion(1:N-1))./dx); %figure out the net flux in each box from sand diffusion, independent of deposited sand
        n_dQdX_diffusion(1:N-1) = [n_dQdx_diffusion(1) n_dQdx_diffusion]; %flux sand out of left hand side of box, ultimate downslope position
    
    % Combined Net Fluxes of Sand into Each Cell: Suspended Sediment + Diffusion
    n_dQdx(1:N-1) = n_dQdx_suspended_sand(1:N-1)+n_dQdX_diffusion(1:N-1);%+dQdx_diffusion(1:N-1); % diff Q to get dQdx
    n_dQdx(N) = n_dQdx(N-1); % net flux out of last box is same as box before it, allows model to drain...boundary condition

    % Sand Height Changes
    n_dHdt = n_dQdx; % conservation statement: change in sand height is a function of net flux of sand into box, no source or sink terms
    n_sand_elevation = n_sand_elevation + n_dHdt.*dt; % update total height of sand
    n_Hneg = find(n_sand_elevation<0); % find negative sand heights, shouldn't be any but just in case
    n_sand_elevation(n_Hneg)=0; % set negative sand heights to zero... no such thing as negative sand.
    
    %Calculate Area of Plant that is Covered by Sand using Trapezoidal
    %Integration
    n_plant_area_total(i) = trapz(x(n_grass1_start:n_grass1_end),n_grass_elevation(n_grass1_start:n_grass1_end))-trapz(x(n_grass1_start:n_grass1_end),original_sand_no_grass(n_grass1_start:n_grass1_end)); % Calculate area of the whole plant
    n_sand_area_accumulated_total(i) = trapz(x,n_sand_elevation)-trapz(x,original_sand_no_grass); %calculate area of all sand that has been added to the system since the time loop started
    n_plant_area_exposed(i) = trapz(x(n_grass1_start:n_grass1_end),n_grass_elevation(n_grass1_start:n_grass1_end))-trapz(x(n_grass1_start:n_grass1_end),n_sand_elevation(n_grass1_start:n_grass1_end)); %calculate how much of the plant is not buried under the sand
    n_plant_area_covered_percent(i) = (n_plant_area_total(i)-n_plant_area_exposed(i))/n_plant_area_total(i); %calculate the percent of the plant buried in the sand in decimal form
    
    %Grow the Grass
    if n_grass_height(i-1)<=max_grass_height %if the grass is shorter than the max height, apply the following growth laws
            n_grass_growth(i) = (min_growth*n_grass_height(i-1)/max_grass_height*(1-n_grass_height(i-1)/max_grass_height))*dt;%growth curve for grass G = G(1-B) based on biomass
    else %if the grass is trying to grow larger than the max height, stop it
        n_grass_growth(i) = 0; %no growth past max height
   
    end %end the grass growth loop
    
    %Update Grass Elevation, Grass Height, and other Grass Parameters with Growth
    n_grass_elevation(n_grass1_start:n_grass1_end) = original_sand_no_grass(n_grass1_start:n_grass1_end)+n_grass_height(i-1)+n_grass_growth(i); %update grass elevation with growth this timestep
    n_grass_height(i) = n_grass_elevation(n_grass1_start)-original_sand_no_grass(n_grass1_start); %update grass height with grass growth this time step
    n_grass1_end_elevation = n_grass_elevation(n_grass1_end); % find the elevation of the grass at the end of the tuft
    n_grass1_end_effective_blocking_height = n_sand_elevation(n_grass1_start)+block_factor*n_grass_height(i); % find the elevation of the effective blocking height of the grass
    n_grass1_shadow_end = ceil((n_grass1_end+10/dx*n_grass_height(i))); % grass wind shadow extends ten times the height of the tuft behind the grass
    n_grass1_shadow_length = length(x(n_grass1_end:n_grass1_shadow_end)); %this is the length of the wind shadow in the array
    n_grass1_shadow_end_elevation = original_sand_no_grass(n_grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    %Update Wind Disturbance Height based on New Grass Height
    n_wind_disturbance_height = n_grass1_shadow_end_elevation-((n_grass1_shadow_end*dx*(n_grass1_shadow_end_elevation-n_grass1_end_effective_blocking_height)/((n_grass1_shadow_end*dx)-n_grass1_start*dx)))+((n_grass1_shadow_end_elevation-n_grass1_end_effective_blocking_height)/((n_grass1_shadow_end*dx)-n_grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    %For plotting percentage covered vs growth factor over a predetermined
    %range of percentages
    n_range_of_growth_factors(i) = 1; %growth factor is a increasing with additional coverage
    
    % Harmful Relationship
    % Wind Speed Dependent Suspended Sediments
        % Calculate How Much Sediment to Deposit Along Wind Shadow
        h_percent_max_windspeed = ones(1,N-1); %initialize windspeed factor so that all wind is max windspeed

        for j = 1:(h_grass1_shadow_length) %loop through the wind shadow (zone of deposition), this is the only place we are going to alter the wind speed

            if h_sand_elevation(h_grass1_end+j)<=h_wind_disturbance_height(h_grass1_end+j) %if the sand isn't taller than the grass in a given unit cell
                % Linear Wind Speed Profile
                h_linear_wind_profile = 1/(h_grass1_shadow_end-h_grass1_end); %the slope of the line that goes from the end of the tuft to the end of the wind shadow
                h_percent_max_windspeed(h_grass1_end+j) = h_linear_wind_profile*x(h_grass1_end+j)-h_linear_wind_profile*h_grass1_end; %the line with values from 1->0 over the length of the wind shadow representing how much sand is going to fall out of suspension as a fraction of the total possible amount to fall out of suspension
            end

        end
        % Net Flux of Sand from Suspended Sediments
        h_dQdx_suspended_sand = suspended_sand_factor*(1 - h_percent_max_windspeed); %net flux into each box from suspended load being deposited is dependent on the wind speed profile as affected by the grass tuft
    
    % Diffusion of Sand
        %Flux of Sand from Diffusion
        h_Q_diffusion(1:N-1) = k.*diff(h_sand_elevation(1:N))./dx; %carry soil out of each dx box according to flux eq. for hillslope diffusion
        
        % Net Flux of Sand from Diffusion
        h_dQdx_diffusion(1:N-2) = (1/rho_sand)*(diff(h_Q_diffusion(1:N-1))./dx); %figure out the net flux in each box from sand diffusion, independent of deposited sand
        h_dQdX_diffusion(1:N-1) = [h_dQdx_diffusion(1) h_dQdx_diffusion]; %flux sand out of left hand side of box, ultimate downslope position
    
    % Combined Net Fluxes of Sand into Each Cell: Suspended Sediment + Diffusion
    h_dQdx(1:N-1) = h_dQdx_suspended_sand(1:N-1)+h_dQdX_diffusion(1:N-1);%+dQdx_diffusion(1:N-1); % diff Q to get dQdx
    h_dQdx(N) = h_dQdx(N-1); % net flux out of last box is same as box before it, allows model to drain...boundary condition

    % Sand Height Changes
    h_dHdt = h_dQdx; % conservation statement: change in sand height is a function of net flux of sand into box, no source or sink terms
    h_sand_elevation = h_sand_elevation + h_dHdt.*dt; % update total height of sand
    h_Hneg = find(h_sand_elevation<0); % find negative sand heights, shouldn't be any but just in case
    h_sand_elevation(h_Hneg)=0; % set negative sand heights to zero... no such thing as negative sand.
    
    %Calculate Area of Plant that is Covered by Sand using Trapezoidal
    %Integration
    h_plant_area_total(i) = trapz(x(h_grass1_start:h_grass1_end),h_grass_elevation(h_grass1_start:h_grass1_end))-trapz(x(h_grass1_start:h_grass1_end),original_sand_no_grass(h_grass1_start:h_grass1_end)); % Calculate area of the whole plant
    h_sand_area_accumulated_total(i) = trapz(x,h_sand_elevation)-trapz(x,original_sand_no_grass); %calculate area of all sand that has been added to the system since the time loop started
    h_plant_area_exposed(i) = trapz(x(h_grass1_start:h_grass1_end),h_grass_elevation(h_grass1_start:h_grass1_end))-trapz(x(h_grass1_start:h_grass1_end),h_sand_elevation(h_grass1_start:h_grass1_end)); %calculate how much of the plant is not buried under the sand
    h_plant_area_covered_percent(i) = (h_plant_area_total(i)-h_plant_area_exposed(i))/h_plant_area_total(i); %calculate the percent of the plant buried in the sand in decimal form
    
    h_plant_buried_time = find(h_plant_area_covered_percent>.99);
    
    %Correct for above ground coverage when the plant dies
    h_false_growth1 = find(h_plant_area_covered_percent<=0); %find negative coverage
    h_false_growth2 = find(h_plant_area_covered_percent>.99); %find coverage >100%
    h_plant_area_covered_percent(h_false_growth1) = 1; %replace with 100% coverage
    h_plant_area_covered_percent(h_false_growth2) = 1; %replace with 100% coverage
    
    %Grow the Grass
    if h_plant_area_covered_percent(i) == 1 %if the plant is totally buried
        h_grass_growth(i) = 0; %it's not growing
        h_biomass_growth(i)=0; %no growth
        h_burial_growth_factor(i) = 0; %no growth factor either
    else
        if h_grass_height(i-1)<=max_grass_height %if the grass is shorter than the max height, apply the following growth laws

                h_grass_growth(i) = (min_growth*h_grass_height(i-1)/max_grass_height*(1-h_grass_height(i-1)/max_grass_height)*(1-h_plant_area_covered_percent(i)/plant_death_coverage_percent))*dt;%growth curve for grass G = c*G(1-B) scaled by covered percentage factor, c
                h_biomass_growth(i) = (min_growth*h_grass_height(i-1)/max_grass_height*(1-h_grass_height(i-1)/max_grass_height))*dt;%growth curve for grass G = G(1-B) only considering the biomass terms, not coverage by sand
                h_burial_growth_factor(i) = (1-h_plant_area_covered_percent(i)); %the real time burial growth factor affecting the plant is the c term, here it hurts the plant with increasing burial


        else %if the grass is trying to grow larger than the max height, stop it
                h_grass_growth(i) = 0; %no growth past max height
                h_biomass_growth(i) = 0; %no growth past max height
                h_burial_growth_factor(i) = 0; %no growth past max height

        end %end the grass growth loop
    end
    
    %Update Grass Elevation, Grass Height, and other Grass Parameters with Growth
    h_grass_elevation(h_grass1_start:h_grass1_end) = original_sand_no_grass(h_grass1_start:h_grass1_end)+h_grass_height(i-1)+h_grass_growth(i); %update grass elevation with growth this timestep
    h_grass_height(i) = h_grass_elevation(h_grass1_start)-original_sand_no_grass(h_grass1_start); %update grass height with grass growth this time step
    h_grass1_end_elevation = h_grass_elevation(h_grass1_end); % find the elevation of the grass at the end of the tuft
    h_grass1_end_effective_blocking_height = h_sand_elevation(h_grass1_start)+block_factor*h_grass_height(i); % find the elevation of the effective blocking height of the grass
    h_grass1_shadow_end = ceil((h_grass1_end+10/dx*h_grass_height(i))); % grass wind shadow extends ten times the height of the tuft behind the grass
    h_grass1_shadow_length = length(x(h_grass1_end:h_grass1_shadow_end)); %this is the length of the wind shadow in the array
    h_grass1_shadow_end_elevation = original_sand_no_grass(h_grass1_shadow_end); %this is the elevation of sand at the end of the wind shadow
    
    %Update Wind Disturbance Height based on New Grass Height
    h_wind_disturbance_height = h_grass1_shadow_end_elevation-((h_grass1_shadow_end*dx*(h_grass1_shadow_end_elevation-h_grass1_end_effective_blocking_height)/((h_grass1_shadow_end*dx)-h_grass1_start*dx)))+((h_grass1_shadow_end_elevation-h_grass1_end_effective_blocking_height)/((h_grass1_shadow_end*dx)-h_grass1_start*dx))*x; %disturbance to wind pattern is a linear decrease over the length of the wind shadow from the tip of the grass to where it intersects the sand
    
    
    %For plotting percentage covered vs growth factor over a predetermined
    %range of percentages
    h_range_of_growth_factors(i) = (1-0.01*range_of_percents(i)/plant_death_coverage_percent); %growth factor is a increasing with additional coverage
    
    
    %Plotting the Animation
    %Plot the Results Each Time Step   
    if (rem(t(i),tplot)==0) %only plot if it is a multiple of tplot
        figure(1) 
        clf
        
        % Beneficial Relationship with Burial
        subplot(3,1,1)
        plot(x,b_grass_elevation, 'k') %plot the grass elevation
        hold all
        plot(x,b_sand_elevation,'k') %plot sand elevation on top of it

        %Color the grass green, do this first so sand can cover grass
        b_uu = [x,x];        % repeat x values
        b_gg = [b_grass_elevation,bottomline];   % vector of upper & lower boundaries
        fill(b_uu,b_gg,[.643,.686,.416]) %fill the polygon created

        %Color the sand beige, do this second so sand can cover the grass
        b_ss = [b_sand_elevation,bottomline];   % vector of upper & lower boundaries
        fill(b_uu,b_ss,[.91,.89,.87]) %fill the polygon created

        %Plot formatting
        title('Benefical Grass Burial Beach Evolution');
        %xlabel('Distance Along Profile (m)');
        ylabel('Elevation (m)');
        set(gca,'fontsize',16,'fontname','arial')
        ht=text(4/5*xmax,.8,['  ',num2str(floor(round(t(i)/7))), ' weeks  '],'fontsize',18); %print time in animation
        axis([0 xmax 0 1.3])
       
        %Neutral Relationship with Burial
        subplot(3,1,2)
        plot(x,n_grass_elevation, 'k') %plot the grass elevation
        hold all
        plot(x,n_sand_elevation,'k') %plot sand elevation on top of it

        %Color the grass green, do this first so sand can cover grass
        n_uu = [x,x];        % repeat x values
        n_gg = [n_grass_elevation,bottomline];   % vector of upper & lower boundaries
        fill(n_uu,n_gg,[.643,.686,.416]) %fill the polygon created

        %Color the sand beige, do this second so sand can cover the grass
        n_ss = [n_sand_elevation,bottomline];   % vector of upper & lower boundaries
        fill(n_uu,n_ss,[.91,.89,.87]) %fill the polygon created

        %Plot formatting
        title('Neutral Grass Burial Beach Evolution');
        %xlabel('Distance Along Profile (m)');
        ylabel('Elevation (m)');
        set(gca,'fontsize',16,'fontname','arial')
        ht=text(4/5*xmax,.8,['  ',num2str(floor(round(t(i)/7))), ' weeks  '],'fontsize',18); %print time in animation
        axis([0 xmax 0 1.3])
        
        % Harmful Relationship with Burial
        subplot(3,1,3)
        plot(x,h_grass_elevation, 'k') %plot the grass elevation
        hold all
        plot(x,h_sand_elevation,'k') %plot sand elevation on top of it

        %Color the grass green, do this first so sand can cover grass
        h_uu = [x,x];        % repeat x values
        h_gg = [h_grass_elevation,bottomline];   % vector of upper & lower boundaries
        fill(h_uu,h_gg,[.643,.686,.416]) %fill the polygon created

        %Color the sand beige, do this second so sand can cover the grass
        h_ss = [h_sand_elevation,bottomline];   % vector of upper & lower boundaries
        fill(h_uu,h_ss,[.91,.89,.87]) %fill the polygon created

        %Plot formatting
        title('Harmful Grass Burial Beach Evolution');
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
%Clean up Plots when Plant Dies
for i = h_plant_buried_time(3):imax
        h_grass_growth(i) = 0; %it's not growing
        h_biomass_growth(i)=0; %no growth
        h_burial_growth_factor(i) = 0; %no growth factor either
        h_grass_height(i) = 0;
end

figure(2) %These are plots for measureable values that change with time
clf

%Total Sand Accumulated in System
subplot(3,1,1)
plot(t,b_sand_area_accumulated_total,'b',t,n_sand_area_accumulated_total,'k',t,h_sand_area_accumulated_total,'r')
hold all
% plot(t,n_sand_area_accumulated_total,'k')
% plot(t,h_sand_area_accumulated_total,'r')
legend('Beneficial Relationship with Burial','Neutral Relationship with Burial','Harmful Relationship with Burial','Location','northwest')
xlabel('Time (days)');
ylabel('(m^2/unit width)');
title('Total Accumulated Sand (m^2/unit width)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 2])

%Plant Area Exposed
subplot(3,1,2)
plot(t,b_plant_area_exposed,'b')
hold all
plot(t,n_plant_area_exposed,'k')
plot(t,h_plant_area_exposed,'r')
xlabel('Time (days)');
ylabel('(m^2/unit width)');
title('Plant Area Exposed (m^2/unit width)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 .5])

%Grass Height
subplot(3,1,3)
plot(t,b_grass_height,'b')
hold all
plot(t,n_grass_height,'k')
plot(t,h_grass_height,'r')
xlabel('Time (days)');
ylabel('(m)');
title('Grass Height (m)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 1.1])


figure(3) %these are plots showing aspects of the sand burial growth factor relationship
clf

%Plant Area Percent Covered
subplot(3,1,1)
plot(t,b_plant_area_covered_percent*100,'b',t,n_plant_area_covered_percent*100,'k',t,h_plant_area_covered_percent*100,'r')
hold all
% plot(t,n_plant_area_covered_percent*100,'k')
% plot(t,h_plant_area_covered_percent*100,'r')
legend('Beneficial Relationship with Burial','Neutral Relationship with Burial','Harmful Relationship with Burial','Location','northwest')
xlabel('Time (days)');
ylabel('(%)');
title('Grass Area Buried (%)');
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax 0 100])

%Plant Growth Curves
subplot(3,1,2)
plot(t,b_grass_growth*1000,'b')
hold all
plot(t,n_grass_growth*1000,'k')
plot(t,h_grass_growth*1000,'r')
plot(t,range_of_zeros_for_time,'--g')
xlabel('Time (days)');
ylabel('(mm/day)');
title('Vertical Grass Growth (mm/day)');
%legend('Total Grass Growth per day','Grass Growth without Burial Growth Factor')
set(gca,'fontsize',16,'fontname','arial')
axis([0 tmax -.5 2])

%Plant Growth Factors
subplot(3,1,3)
plot(range_of_percents,b_range_of_growth_factors,'b')
hold all
plot(range_of_percents,n_range_of_growth_factors,'k')
plot(range_of_percents,h_range_of_growth_factors,'r')
xlabel('Grass Area Buried (% of total grass area)');
ylabel('Burial Growth Factor');
title('Effect of Burial on Grass Growth')
set(gca,'fontsize',16,'fontname','arial')
axis([0 20 -.1 2.1])

% Final Dune Profile Comparison
figure(4)
clf
plot(x,b_sand_elevation,'b',x,n_sand_elevation,'k',x,h_sand_elevation,'r',x,original_sand_no_grass,'g')
xlabel('Distance Along Profile (m)');
ylabel('Sand Elevation (m)');
title('Final Dune Shapes')
set(gca,'fontsize',16,'fontname','arial')
legend('Beneficial Relationship with Burial','Neutral Relationship with Burial','Harmful Relationship with Burial','Original Beach Profile','Location','northwest')
axis([0 xmax 0 0.5])