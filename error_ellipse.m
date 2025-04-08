function rotation_ang = error_ellipse(C, mu, varargin)
    % Graficar una elipse de error basada en una matriz de covarianza (C) y una media (mu)
    % Los parámetros adicionales especifican el nivel de confianza.
    p = inputParser;
    addOptional(p, 'conf', 0.95);
    addOptional(p, 'plt', 1);
    parse(p, varargin{:});
    conf = p.Results.conf;
    plt = p.Results.plt;

    % Calcular el eje mayor y menor de la elipse
    [V, D] = eig(C);
    [d, idx] = sort(diag(D), 'descend');
    D = diag(D);
    D = D(idx);
    V = V(:, idx);

    % Calcular el ángulo de rotación
    theta = atan2(V(2,1), V(1,1));

    % Transformar el ángulo de rotación a grados
    rotation_ang = rad2deg(theta);

    % Calcular los parámetros de la elipse
    chi2_val = sqrt(chi2inv(conf, 2));
    a = chi2_val * sqrt(D(1,1));
    b = chi2_val * sqrt(D(2,1));

    % Parámetros de la elipse en coordenadas polares
    t = linspace(0, 2*pi, 100);
    X = a * cos(t);
    Y = b * sin(t);

    % Rotar y desplazar la elipse
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    ellipse = R * [X; Y];

    if plt == 1
        plot(ellipse(1,:) + mu(1), ellipse(2,:) + mu(2), 'r-');
    end
end
