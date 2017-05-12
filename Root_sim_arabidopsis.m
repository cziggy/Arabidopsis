% Initial conditions
g = 0.55; % growth rate
%d = 2; % branching distance above root tip
Initial_length = 0.1;
alpha = 2/3; % Branching factor
n = 1; %Number of branches
Initial_width = 0.8; % Width of new roots
Thickening_rate = 0.1; % Thickening rate
Max_width = 2; % Maximum width
Apical_length = 2; %Length of apical zone
Basal_length = 0; %Length of basal zone
Max_length = 1000;

Max_branches = 50; %Maximum number of branches
Time = 0; %Time counter
Time_max = 15;

B = [0.38 0]; % Vector of branching rates by root order (i.e B(i) is the branching rate for a root of order i)
BR = zeros(1, Max_branches); % Vector of branching rates for each root (i.e BR(i) is the branching rate for root i)
L = Initial_length; % Vector of lengths for all branches

%Define initial branch
Branch(1).Time = 0; %Time of root creation
Branch(1).Angle = 3*pi/2; %Branch angle taken anticlockwise from x axis. 
Branch(1).Length = Initial_length; %Length of initial root
Branch(1).Width = 0.01; %Initial diameter of root
Branch(1).Order = 1; %Order of initial root
Branch(1).Branching_rate = B(Branch(1).Order); %Branching rate of root
Branch(1).Max_width = Max_width;
Branch(1).Initial_width = Initial_width;
Branch(1).Thickening_rate = Thickening_rate;
Branch(1).Number_of_daughters = 0; 
Branch(1).x1 = 0;
Branch(1).y1 = 0;
BR(1) = Branch(1).Branching_rate;
Branch(1).Daughter_branches = zeros(1,Max_branches);
Branch(1).Inter_branch_distances = zeros(1,Max_branches);

while Time < Time_max
    
    %Find time until next branching event using Gillespie algorithm
    B1 = BR; %Create vector of branching rates B1 for branches that can branch. 
 
    i = 1;
    while i < n + 1
    branch_length = Branch(i).Length;
    if branch_length < Apical_length + Basal_length
        B1(i) = 0;
    end
    i = i + 1;
    end
    
    R = rand(4,1); % Generate random numbers
     
    a = sum(B1); % Calculate sum of branching rates
    
    if a > 0

        deltat = -log(R(1))./a; % Time until next reaction

        Time = Time + deltat; % Increase time counter

        %Update properties of root system over time deltat
        i = 1; 
        while i < n + 1
            %Length:
            L0 = Branch(i).Length; 
            %L = L0*exp(g*deltat)
            if L0 > Max_length
                L = Max_length;

            else
                L = L0 + g*deltat;
            end
            Branch(i).Length = L;

            %Width:
            a = -log(Branch(i).Max_width/Branch(i).Initial_width -1);
            Branch(i).Width = Branch(i).Max_width./(1+exp(-(a + Branch(i).Thickening_rate*((Time-Branch(i).Time))))); 
            i = i+1;
        end

        %Find which meristem will branch
        x = R(2)*sum(B1);

        c = 0; 
        i = 0;

        while c < x
            i = i + 1;
            c = c + B1(i);
        end

        Branching_meristem = i; % The branching meristem
        
        %Create new branch
        Branch(n+1).Start_time = Time;

        %Branching direction related to angle of mother branch
        Mother_angle = Branch(Branching_meristem).Angle;

        if Mother_angle > 3*pi/2
            if R(3) > 0.5
                Branch(n+1).Angle = Mother_angle + (2*pi-Mother_angle)/2;
            else
                Branch(n+1).Angle = Mother_angle - (Mother_angle -3*pi/2)/2;
            end

        elseif Mother_angle < 3*pi/2
            if R(3) < 0.5
                Branch(n+1).Angle = Mother_angle - (Mother_angle -pi)/2;
            else
                Branch(n+1).Angle = Mother_angle + (3*pi/2-Mother_angle)/2;
            end
        else
            if R(3) > 0.5
                Branch(n+1).Angle = 7*pi/4;
            else
                Branch(n+1).Angle = 5*pi/4;
            end
        end   

        Branch(n+1).Length = Initial_length;
        Branch(n+1).Width = Initial_width;
        Branch(n+1).Mother_branch = Branching_meristem;
        Branch(n+1).Order = Branch(Branching_meristem).Order + 1;
        Branch(n+1).Branching_rate = BR(Branch(n+1).Order);
        Branch(n+1).Mother_length_at_branching = Branch(Branching_meristem).Length;
        Branch(n+1).Max_width = Max_width;
        Branch(n+1).Initial_width = alpha*Branch(Branching_meristem).Width;
        Branch(n+1).Thickening_rate = Thickening_rate;
        Branch(n+1).Time = Time;
        
        if B(Branch(n+1).Order) > 0
            Branch(n+1).Daughter_branches = zeros(1,Max_branches);
            Branch(n+1).Inter_branch_distances = zeros(1,Max_branches);
        end
        
        %Branch(n+1).x1 = Branch(Branching_meristem).x1 + Branch(Branching_meristem).Length * cos(Branch(Branching_meristem).Angle);
        %Branch(n+1).y1 = Branch(Branching_meristem).y1 + Branch(Branching_meristem).Length * sin(Branch(Branching_meristem).Angle);
        R = lognrnd(1.8,1.45)/100;

        while R > 1
             R = lognrnd(1.8,1.45)/100;
        end

        %branching can occur anywhere along the branching zone
        Branch(n+1).x1 = Branch(Branching_meristem).x1 + (Basal_length + R*(Branch(Branching_meristem).Length - Basal_length - Apical_length))*cos(Branch(Branching_meristem).Angle);
        Branch(n+1).y1 = Branch(Branching_meristem).y1 + (Basal_length + R*(Branch(Branching_meristem).Length - Basal_length - Apical_length))*sin(Branch(Branching_meristem).Angle);
        Branch(n+1).Distance_along_mother = ((Branch(n+1).x1-Branch(Branching_meristem).x1)^2 + (Branch(n+1).y1 - Branch(Branching_meristem).y1)^2)^0.5;

        BR(n+1) = Branch(n+1).Branching_rate; %Update vector of branching rates
        
        % Update properties of mother branch
        Branch(Branching_meristem).Number_of_daughters = Branch(Branching_meristem).Number_of_daughters + 1;
        Branch(Branching_meristem).Daughter_branches(Branch(Branching_meristem).Number_of_daughters) = n+1;
        
        %Update inter-branch distances for mother branch
        
        Distance = Branch(n+1).Distance_along_mother
        
        if Branch(n+1).Distance_along_mother > sum(Branch(Branching_meristem).Inter_branch_distances) %If the new branch is the furthest along the mother
            j = Branch(Branching_meristem).Number_of_daughters
            
            if j > 1
                Branch(Branching_meristem).Inter_branch_distances(j) = Branch(n+1).Distance_along_mother - Branch(Branching_meristem).Inter_branch_distances(j - 1);
            else
                Branch(Branching_meristem).Inter_branch_distances(1) = Branch(n+1).Distance_along_mother;
            end
            
        else  
            
        x = 0;
        i = 1;
            % Find where the new branch should fit in the inter-branch
            % distances matrix
            while x < Branch(n+1).Distance_along_mother
                i = i + 1
                IBD = Branch(Branching_meristem).Inter_branch_distances
                x = x + Branch(Branching_meristem).Inter_branch_distances(i)

               
            end
            
            j = Branch(Branching_meristem).Number_of_daughters
            
            %Move all elements above i-1 up one
            while j > i - 1
                Branch(Branching_meristem).Inter_branch_distances(j+1) = Branch(Branching_meristem).Inter_branch_distances(j);
                a = Branch(Branching_meristem).Inter_branch_distances(j+1)
                j = j - 1;
            end
            
            %Insert IBD for new root
            Branch(Branching_meristem).Inter_branch_distances(i) = Branch(n+1).Distance_along_mother - Branch(Branching_meristem).Inter_branch_distances(i-1);
            b = Branch(Branching_meristem).Inter_branch_distances(i)
            % Update IBD for all roots below the new one
            
            j = i + 1; 
            
            while j < Branch(Branching_meristem).Number_of_daughters + 1
                Branch(Branching_meristem).Inter_branch_distances(j) = Branch(Branching_meristem).Inter_branch_distances(j) - Branch(Branching_meristem).Inter_branch_distances(j-1);
                c = Branch(Branching_meristem).Inter_branch_distances(j) 
                j = j + 1; 
            end
        end
        n = n+1; %Increase branch counter
    
    else
        % System grows until a root can branch
        
        Growth_times = zeros(n,1);
        % Find roots with a non-zero branching rate
        i = 1;
        while i < n + 1
            
            x = BR(i);
            
            if x > 0
                % Calculate the time for the root to grow long enough to
                % branch
                Growth_times(i) = (Basal_length + Apical_length)/g;
            end
            
            i = i+1;
        end
        
       deltat = min(Growth_times(Growth_times > 0)); %Find the shortest time until a root can branch
        
       Time = Time + deltat; % Increase time counter
       
        %Update properties of root system over time deltat
        i = 1; 
        while i < n+1
            %Length:
            L0 = Branch(i).Length; 
            %L = L0*exp(g*deltat)
            if L0 > Max_length
                L = Max_length;

            else
                L = L0 + g*deltat;
            end
            Branch(i).Length = L;

            %Width:
            a = -log(Branch(i).Max_width/Branch(i).Initial_width -1);
            Branch(i).Width = Branch(i).Max_width./(1+exp(-(a + Branch(i).Thickening_rate*((Time-Branch(i).Time))))); 
            i = i+1;
        end
        
    end
    
end

% Create vector of inter-branch distances for each mother root

%Find the end points for all branches and plot the system
i = 1; 
while i < n + 1
    Branch(i).x2 = Branch(i).x1 + Branch(i).Length * cos(Branch(i).Angle);
    Branch(i).y2 = Branch(i).y1 + Branch(i).Length * sin(Branch(i).Angle);
    
    x1 = Branch(i).x1;
    y1 = Branch(i).y1;
    x2 = Branch(i).x2;
    y2 = Branch(i).y2;
    
    x = linspace(x1, x2);
    y = ((y2-y1)/(x2-x1))*(x-x1)+y1;
    
    Thickness = 2*Branch(i).Width;
    plot(x,y, 'LineWidth', Thickness);
    hold on
    i = i + 1;
end

