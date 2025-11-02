function shape_factor = nye_shape_factor(type_valley,half_width,thickness)
% Date: 2024/02/02

ratio = half_width / thickness;
switch type_valley
    case 'parabola'
        if rario<=1.0
            shape_factor = 0.445;
        elseif ratio>1.0 && ratio<=2.0
            shape_factor = 0.646;
        elseif ratio>2.0 && ratio<=3.0
            shape_factor = 0.746;
        elseif ratio>3.0 && ratio<=4.0
            shape_factor = 0.806;
        else
            shape_factor = 1.0;
        end
    case 'semi_ellipse'
        if rario<=1.0
            shape_factor = 0.500;
        elseif ratio>1.0 && ratio<=2.0
            shape_factor = 0.709;
        elseif ratio>2.0 && ratio<=3.0
            shape_factor = 0.799;
        elseif ratio>3.0 && ratio<=4.0
            shape_factor = 0.849;
        else
            shape_factor = 1.0;
        end
    case 'rectangle'
        if rario<=1.0
            shape_factor = 0.558;
        elseif ratio>1.0 && ratio<=2.0
            shape_factor = 0.789;
        elseif ratio>2.0 && ratio<=3.0
            shape_factor = 0.884;
        else
            shape_factor = 1.0;
        end
    otherwise
        disp('The type of valley shape must be parabola, semi_ellipse, rectangle')
end

end
