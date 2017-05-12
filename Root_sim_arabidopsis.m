% Initial conditions
g = (1.4)*rand; % growth rate
b = (1.4)*rand;%d = 2; % branching distance above root tip
U = rand;

Initial_Length = 0.1;
alpha = 2/3; % Branching factor
n = 1; %Number of branches
Initial_Width = 0.8; % Width of new roots
Thickening_Rate = 0.1; % Thickening rate
Max_Width = 2; % Maximum width
Apical_length = 2; %Length of apical zone
Basal_length = 0; %Length of basal zone
Max_Length = 1000;

Number_of_Daughters_Data = [9 7 11 12 13 19 15 9 16 12]; %Update data for new time check
Length_Data = [15.94554349 13.03006591 13.97050182 14.11069989 16.70318919 18.23718774 15.32863698 14.59489812 16.34310175 13.90157061];
epsilon_number = 6;
epsilon_length = 6;

Time_max = 15;
Time = 0; %Time counter
m_max = 10;

%Root properties by root order
Branching_Rate = [b 0]; % Vector of branching rates by root order (i.e B(i) is the branching rate for a root of order i)


Max_Branches = ceil(2*Time_max*max(Branching_Rate)); %Maximum number of branches that can be created

%Preallocate vectors
Branch_Start_Time = zeros(Max_Branches, m_max);
Branch_Angle = zeros(Max_Branches,m_max);
Branch_Length = zeros(Max_Branches,m_max);
Branch_Width = zeros(Max_Branches,m_max);
Branch_Order = zeros(Max_Branches,m_max);
Branch_Number_of_Daughters = zeros(nnz(Branching_Rate),m_max);
Branch_Mother_Branch = zeros(Max_Branches,m_max);
Branch_Mother_Length_at_Branching = zeros(Max_Branches,m_max);
Branch_Distance_Along_Mother = zeros(Max_Branches,m_max);
Branch_x1 = zeros(Max_Branches,m_max);
Branch_y1 = zeros(Max_Branches,m_max);
Branch_Daughter_Branches = zeros(nnz(Branching_Rate),Max_Branches,m_max);
Branch_Inter_Branch_Distances = zeros(nnz(Branching_Rate),Max_Branches,m_max);

BR = zeros(Max_Branches,m_max); % Vector of branching rates for each root (i.e BR(i) is the branching rate for root i)

m = 1;
while m < m_max + 1
    Time = 0;
    n = 1;
    
    %Define initial branch
    Branch_Angle(1,m) = 3*pi/2; %Branch angle taken anticlockwise from x axis. 
    Branch_Length(1,m) = 4.38; %Length of initial root
    Branch_Width(1,m) = 0.01; %Initial diameter of root
    Branch_Order(1,m) = 1;
    
    BR(1,m) = Branching_Rate(Branch_Order(1,1));
    
    while Time < Time_max

        %Find time until next branching event using Gillespie algorithm
        B1 = BR(:,m); %Create vector of branching rates B1 for branches that can branch. 

        i = 1;
        while i < n + 1

            if Branch_Length(i,m) < Apical_length + Basal_length
                B1(i) = 0;
            end
        i = i + 1;
        end

        R = rand(6,1); % Generate random numbers

        a = sum(B1); % Calculate sum of branching rates

        if a > 0

            deltat = -log(R(1))./a; % Time until next reaction

            Time = Time + deltat; % Increase time counter

            %Update properties of root system over time deltat
            i = 1; 
            while i < n + 1
                %Length:
                if Branch_Length(i,m) < Max_length
                    Branch_Length(i,m) = Branch_Length(i,m) + g*deltat;
                end

                %Width:
                %Branch_Width(i,1) = Max_Width(Branch_Order(i,1))./(1+exp(-(-log(Max_Width(Branch_Order(i,1))/Initial_Width(Branch_Order(i,1)) -1)+ Branch(i).Thickening_rate*((Time-Branch(i).Time))))); 
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
            Branch_Start_Time(n+1,m) = Time;

            %Branching direction related to angle of mother branch
            Mother_angle = Branch_Angle(Branching_meristem,m);

            if Mother_angle > 3*pi/2
                if R(3) > 0.5
                    Branch_Angle(n+1,m) = Mother_angle + (2*pi-Mother_angle)/2;
                else
                    Branch_Angle(n+1,m) = Mother_angle - (Mother_angle -3*pi/2)/2;
                end

            elseif Mother_angle < 3*pi/2
                if R(3) < 0.5
                    Branch_Angle(n+1,m) = Mother_angle - (Mother_angle -pi)/2;
                else
                    Branch_Angle(n+1,m) = Mother_angle + (3*pi/2-Mother_angle)/2;
                end
            else
                if R(3) > 0.5
                    Branch_Angle(n+1,m) = 7*pi/4;
                else
                    Branch_Angle(n+1,m) = 5*pi/4;
                end
            end   

            Branch_Length(n+1,m) = Initial_length;
            Branch_Width(n+1,m) = Initial_width;
            Branch_Mother_Branch(n+1,m) = Branching_meristem;
            Branch_Order(n+1,m) = Branch_Order(Branching_meristem,m) + 1;
            Branch_Mother_Length_at_Branching(n+1,m) = Branch_Length(Branching_meristem,m);

            if U<0.5
                d = rand;
                uniform = 1;
                
            else
                d = gamrnd(1.7786,0.194183);
                uniform = 0;
                            
            end
            
            %branching can occur anywhere along the branching zone
            Branch_x1(n+1,m) = Branch_x1(Branching_meristem,m) + (Basal_length + d*(Branch_Length(Branching_meristem,m) - Basal_length - Apical_length))*cos(Branch_Angle(Branching_meristem,m));
            Branch_y1(n+1,m) = Branch_y1(Branching_meristem,m) + (Basal_length + d*(Branch_Length(Branching_meristem,m) - Basal_length - Apical_length))*sin(Branch_Angle(Branching_meristem,m));
            Branch_Distance_Along_Mother(n+1,m) = ((Branch_x1(n+1,m)-Branch_x1(Branching_meristem,m))^2 + (Branch_y1(n+1,m) - Branch_y1(Branching_meristem,m))^2)^0.5;

            % Update properties of mother branch
            Branch_Number_of_Daughters(Branching_meristem,m) = Branch_Number_of_Daughters(Branching_meristem,m) + 1;
            Branch_Daughter_Branches(Branch_Number_of_Daughters(Branching_meristem,m),Branching_meristem,m) = n+1;

            %Update inter-branch distances for mother branch
            if Branch_Distance_Along_Mother(n+1,m) > sum(Branch_Inter_Branch_Distances(Branching_meristem,:,m)) %If the new branch is the furthest along the mother
                j = Branch_Number_of_Daughters(Branching_meristem,m);

                if j > 1
                    x = sum(Branch_Inter_Branch_Distances(Branching_meristem,:,m));

                    Branch_Inter_Branch_Distances(Branching_meristem,j,m) = Branch_Distance_Along_Mother(n+1,m) - x;
                else

                Branch_Inter_Branch_Distances(Branching_meristem,1,m) = Branch_Distance_Along_Mother(n+1,m);

                end  

            else  

                y = 0;
                i = 0;
                % Find where the new branch should fit in the inter-branch
                % distances matrix

                while y < Branch_Distance_Along_Mother(n+1,m)
                    i = i + 1;
                    y = y + Branch_Inter_Branch_Distances(Branching_meristem,i,m);

                end

                j = Branch_Number_of_Daughters(Branching_meristem,m);

                %Move all elements above i-1 up one
                while j > i
                    Branch_Inter_Branch_Distances(Branching_meristem,j,m) = Branch_Inter_Branch_Distances(Branching_meristem,j-1,m);
                    j = j - 1;
                end

                %Insert IBD for new root
                if i>1
                    %Find the sum of all IBDs before new branch
                    k = 1;
                    y = 0; 

                    while k < i
                        y = y + Branch_Inter_Branch_Distances(Branching_meristem,k,m);
                        k = k + 1;
                    end

                    Branch_Inter_Branch_Distances(Branching_meristem,i,m) = Branch_Distance_Along_Mother(n+1,m) - y;
                else
                    Branch_Inter_Branch_Distances(Branching_meristem,i,m) = Branch_Distance_Along_Mother(n+1,m);
                end


             % Update IBD for root below the new one

             Branch_Inter_Branch_Distances(Branching_meristem,i+1,m) =Branch_Inter_Branch_Distances(Branching_meristem,i+1,m) - Branch_Inter_Branch_Distances(Branching_meristem,i,m);

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
            while i < n + 1
                %Length:
                if Branch_Length(i,m) < Max_length
                    Branch_Length(i,m) = Branch_Length(i,m) + g*deltat;
                end

                %Width:
                %Branch_Width(i,1) = Max_Width(Branch_Order(i,1))./(1+exp(-(-log(Max_Width(Branch_Order(i,1))/Initial_Width(Branch_Order(i,1)) -1)+ Branch(i).Thickening_rate*((Time-Branch(i).Time))))); 
                i = i+1;
            end


        end

    end
    % Only increase m if results fit with data
    x = abs(Branch_Number_of_Daughters(1,m)-Number_of_Daughters_Data(1));
                if x < epsilon_number %Check number of branches

                    y = abs(Branch_Length(1,m) - Length_Data(m));
                    
                    if y < epsilon_length %Check lengths
                        %m_increases = m_increases+1;
                        m = m + 1; %Only increase m if last m was a hit
                    else
                        m = m_max + 1;
                    end
                else
                    m = m_max + 1;
                end
               
end
            

