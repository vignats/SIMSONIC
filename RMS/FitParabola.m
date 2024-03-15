function [surface] = FitParabola(binary_image, boundary_endost_init, print)
% FITPARABOLA Fits a parabolic curve to the provided boundary points of the endosteum.
%
%   surface = FitParabola(binary_image, boundary_endost_init) fits a parabolic curve to
%   the initial boundary points of the endosteum obtained from the binary image. It returns
%   the coordinates of the fitted parabola.
%
%   Inputs:
%   - binary_image: Binary image containing the bone structure.
%   - boundary_endost_init: Initial boundary points of the endosteum.
%
%   Output:
%   - surface: Matrix containing the coordinates of the fitted parabola.
%
%   Example:
%       surface = FitParabola(binary_image, boundary_endost_init)
%
%   See also: fit

    model = fit(boundary_endost_init(:, 1),  boundary_endost_init(:, 2), 'poly2');
    y = feval(model, boundary_endost_init(:, 1));  % Y coordinate of the parabola
    surface = [boundary_endost_init(:, 1), y];     % Associate the Y coordinate of the parabola with the X coordinate
    
    % Plot final curve
    if print
        figure
        imshow(binary_image);
        hold on;
        plot(boundary_endost_init(:, 1), boundary_endost_init(:, 2), 'b.');  % Trace le contour
        hold on;
        plot(boundary_endost_init(:, 1), y, 'r-');  % Trace la parabole ajustée
        xlabel('X');
        ylabel('Y');
        legend('Contour', 'Parabole');
        title('Fitted parabola to the endost boundary')
    end
end