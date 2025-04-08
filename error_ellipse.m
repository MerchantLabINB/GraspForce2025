function rotation_ang = error_ellipse(C, mu, varargin)
    % Plots an error ellipse based on a covariance matrix (C) and a mean value (mu). 
    % Additional parameters set confidence level.
    p = inputParser;
    addOptional(p, 'conf', 0.95);
    addOptional(p, 'plt', 1);
    parse(p, varargin{:});
    conf = p.Results.conf;
    plt = p.Results.plt;

    % Estimates major and minor axes of the ellipse
    [V, D] = eig(C);
    [d, idx] = sort(diag(D), 'descend');
    D = diag(D);
    D = D(idx);
    V = V(:, idx);

    % Rotation angle
    theta = atan2(V(2,1), V(1,1));

    % Transforms rotation angle to to degrees
    rotation_ang = rad2deg(theta);

    % Estimate ellipse parameters
    chi2_val = sqrt(chi2inv(conf, 2));
    a = chi2_val * sqrt(D(1,1));
    b = chi2_val * sqrt(D(2,1));

    % Ellipse parameters in polar coordinates
    t = linspace(0, 2*pi, 100);
    X = a * cos(t);
    Y = b * sin(t);

    % Rotate and displace the ellipse
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    ellipse = R * [X; Y];

    if plt == 1
        plot(ellipse(1,:) + mu(1), ellipse(2,:) + mu(2), 'r-');
    end
end
