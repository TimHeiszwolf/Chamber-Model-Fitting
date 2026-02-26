function sys_set = setParameterBounds(sys_init, parameters, margin_factor)
    %SETPARAMETERBOUNDS Sets parameter bounds and 'Free' status for a system.
    %
    %   sys_set = setParameterBounds(sys_init, parameters, margin_factor)
    %
    %   Updates the Minimum, Maximum, and Free properties of the parameters
    %   in an idgrey model. The bounds are calculated by dividing/multiplying
    %   the nominal parameter values by the provided margin factors.
    %
    %   Inputs:
    %       sys_init      - The initial idgrey model object.
    %       parameters    - An Nx2 cell array where column 2 contains the
    %                       nominal parameter values.
    %       margin_factor - Defines the width of the bounds. Can be:
    %           * Scalar: Applies the same factor to all parameters.
    %             (Min = Val/Factor, Max = Val*Factor).
    %           * Vector (Nx1): Applies a specific factor per parameter.
    %             (Min = Val/Factor(i), Max = Val*Factor(i)).
    %           * Matrix (Nx2): Applies specific lower and upper factors.
    %             (Min = Val/Factor(i,1), Max = Val*Factor(i,2)).
    %
    %   Note: If a margin factor is exactly 1, the parameter is fixed (Free = false).

    arguments
        sys_init
        parameters
        margin_factor = 10%For example [100; 100; 100; 100; 100; 100; 2; 2; 10; 10; 10; 10]; or [100, 100; 100, 100; 100, 100; 100, 100; 100, 100; 100, 100; inf, 2; inf, 2; 10, 10; 10, 10; 10, 10; 10, 10]; (where the first value is the minimum and the second the maximum)
    end
    
    number_of_parameters = size(parameters, 1);

    for i = 1:number_of_parameters

        if isscalar(margin_factor) 
            % Case 1: Scalar - Apply same factor to all
            sys_init.Structure.Parameters(i).Free = true;
            sys_init.Structure.Parameters(i).Minimum = parameters{i,2}/margin_factor;
            sys_init.Structure.Parameters(i).Maximum = parameters{i,2}*margin_factor;

        elseif (size(margin_factor, 2)) == 1 && (size(margin_factor, 1) == number_of_parameters)
            % Case 2: Vector (Nx1) - Apply specific factor to both Min and Max
            sys_init.Structure.Parameters(i).Minimum = parameters{i,2}/margin_factor(i);
            sys_init.Structure.Parameters(i).Maximum = parameters{i,2}*margin_factor(i);
            
            % If factor is one, explicitly set the parameters to not be free.
            if margin_factor(i)==1
                sys_init.Structure.Parameters(i).Free = false;
            else
                sys_init.Structure.Parameters(i).Free = true;
            end

        elseif (size(margin_factor, 2)) == 2 && (size(margin_factor, 1) == number_of_parameters)
            % Case 3: Matrix (Nx2) - Apply different factors to Min and Max
            sys_init.Structure.Parameters(i).Minimum = parameters{i, 2}/margin_factor(i, 1);
            sys_init.Structure.Parameters(i).Maximum = parameters{i, 2}*margin_factor(i, 2);
            
            % If both factors are one, explicitly set the parameters to not be free.
            if (margin_factor(i, 1)==1) && (margin_factor(i, 2)==1)
                sys_init.Structure.Parameters(i).Free = false;
            else
                sys_init.Structure.Parameters(i).Free = true;
            end
        else
            error('Margin factor must be a scalar or match the number of parameters.');
        end
    end

    sys_set = sys_init;
end