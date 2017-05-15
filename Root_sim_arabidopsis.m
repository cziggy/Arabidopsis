%Initial conditions
g = 0.55; % growth rate
b = 0.38;
%d = 2; % branching distance above root tip
Initial_Length = 0.1;
alpha = 2/3; % Branching factor
n = 1; %Number of branches
Initial_Width = 0.8; % Width of new roots
Thickening_Rate = 0.1; % Thickening rate
Max_Width = 2; % Maximum width
Apical_length = 2; %Length of apical zone
Basal_length = 0; %Length of basal zone
Max_Length = 1000; 

Time_max = 15;
Time = 0; %Time counter
m_max = 1;

Element_Length = 0.5;
Branching_Rate = b*ones(1,g*Time_max/Element_Length);

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

%Define initial branch
Branch_Angle(1,:) = 3*pi/2; %Branch angle taken anticlockwise from x axis. 
Branch_Length(1,:) = 4.38; %Length of initial root
Branch_Width(1,:) = 0.01; %Initial diameter of root
Branch_Order(1,:) = 1;

m = 1;

while Time < Time_max
    
    %Find time until next branching event using Gillespie algorithm
       
    R = rand(6,1); % Generate random numbers
     
    a = Branch_Length(1,m); % Calculate sum of branching rates
    
    if (a > Apical_length + Basal_length) || (a == Apical_length + Basal_length)

        deltat = -log(R(1))./sum(Branching_Rate); % Time until next reaction

        Time = Time + deltat; % Increase time counter

        %Update properties of root system over time deltat
        i = 1; 
        while i < n + 1
            %Length:
            if Branch_Length(i,m) < Max_length
                Branch_Length(i,m) = Branch_Length(i,m) + g*deltat;
                Number_of_Elements = Branch_Length(1,m)/Element_Length;
                
                if Number_of_Elements > nnz(Branching_Rate)
                    i = nnz(Branching_Rate);
                    while i < Number_of_Elements
                        i = i + 1;
                        Branching_Rate(i) = b;
                    end
                end
            end

            %Width:
            %Branch_Width(i,1) = Max_Width(Branch_Order(i,1))./(1+exp(-(-log(Max_Width(Branch_Order(i,1))/Initial_Width(Branch_Order(i,1)) -1)+ Branch(i).Thickening_rate*((Time-Branch(i).Time))))); 
            i = i+1;
        end

        %Find which meristem will branch
        x = R(2)*sum(Branching_Rate);

        c = 0; 
        i = 0;

        while c < x
            i = i + 1;
            c = c + Branching_Rate(i);
        end

        Branching_Element = i; % The branching element
        
        %Create new branch
        Branch_Start_Time(n+1,m) = Time;

        %Branching direction related to angle of mother branch
        Mother_angle = Branch_Angle(1,m);

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
        Branch_Mother_Branch(n+1,m) = 1;
        Branch_Order(n+1,m) = Branch_Order(Branching_meristem,m) + 1;
        Branch_Mother_Length_at_Branching(n+1,m) = Branch_Length(Branching_meristem,m);
        
        %Branch(n+1).x1 = Branch(Branching_meristem).x1 + Branch(Branching_meristem).Length * cos(Branch(Branching_meristem).Angle);
        %Branch(n+1).y1 = Branch(Branching_meristem).y1 + Branch(Branching_meristem).Length * sin(Branch(Branching_meristem).Angle);
        R = lognrnd(1.8,1.45)/100;

        while R > 1
             R = lognrnd(1.8,1.45)/100;
        end

        %branching can occur anywhere along the branching zone
        Branch_x1(n+1,m) = Branch_x1(Branching_meristem,m) + (Basal_length + R*(Branch_Length(Branching_meristem,m) - Basal_length - Apical_length))*cos(Branch_Angle(Branching_meristem,m));
        Branch_y1(n+1,m) = Branch_y1(Branching_meristem,m) + (Basal_length + R*(Branch_Length(Branching_meristem,m) - Basal_length - Apical_length))*sin(Branch_Angle(Branching_meristem,m));
        Branch_Distance_Along_Mother(n+1,m) = Branch_Length(1,m);
        
        % Update properties of mother branch
        Branch_Number_of_Daughters(1,m) = Branch_Number_of_Daughters(1,m) + 1;
        Branch_Daughter_Branches(Branch_Number_of_Daughters(1,m),m) = n+1;
        
        %Update inter-branch distances for mother branch
        if Branch_Distance_Along_Mother(n+1,m) > sum(Branch_Inter_Branch_Distances(:,m)) %If the new branch is the furthest along the mother
            j = Branch_Number_of_Daughters(1,m);

            if j > 1
                x = sum(Branch_Inter_Branch_Distances(:,m));

                Branch_Inter_Branch_Distances(j,m) = Branch_Distance_Along_Mother(n+1,m) - x;
            else

            Branch_Inter_Branch_Distances(1,m) = Branch_Distance_Along_Mother(n+1,m);
            
            end  
        
        else  

            y = 0;
            i = 0;
            % Find where the new branch should fit in the inter-branch
            % distances matrix
        
            while y < Branch_Distance_Along_Mother(n+1,m)
                i = i + 1;
                y = y + Branch_Inter_Branch_Distances(i,m);

            end

            j = Branch_Number_of_Daughters(1,m);

            %Move all elements above i-1 up one
            while j > i
                Branch_Inter_Branch_Distances(j,m) = Branch_Inter_Branch_Distances(j-1,m);
                j = j - 1;
            end

            %Insert IBD for new root
            if i>1
                %Find the sum of all IBDs before new branch
                k = 1;
                y = 0; 

                while k < i
                    y = y + Branch_Inter_Branch_Distances(k,m);
                    k = k + 1;
                end

                Branch_Inter_Branch_Distances(i,m) = Branch_Distance_Along_Mother(n+1,m) - y;
            else
                Branch_Inter_Branch_Distances(i,m) = Branch_Distance_Along_Mother(n+1,m);
            end


         % Update IBD for root below the new one
         
         Branch_Inter_Branch_Distances(Branching_meristem,i+1,m) =Branch_Inter_Branch_Distances(Branching_meristem,i+1,m) - Branch_Inter_Branch_Distances(Branching_meristem,i,m);
         
        end

        n = n+1; %Increase branch counter

        
    
    else
        % System grows until a root can branch
        delta_t = ((Basal_length + Apical_length)-Branch_Length(1,m))/g;
        
        Time = Time + delta_t;
        
        Branch_Length(1,m) = Apical_length + Basal_length;
       
    end
    
end
