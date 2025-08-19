function Algomat = ADImethod(u, f, D, params)

% u is the input matrix representing the solution at the previous time step,
% and we use Algomat function to update the current value of u by the
% alternating direction implicit (ADI) algorithm for nonlinear Reaction-Diffusion Equations
% with Neumann boundary conditions (the method easily adapt to the linear equation
% with no difficulty since the nonlinear ones possibly include linear ones
% as the starting point).
% f is the right-hand side representing the reaction part for the equation.
% u and f have the same size.
% D is the diffusion coefficient, and must be nonzero in this
% implementation.

% Here, we use the same spacing for x and y direction, and the number of
% grid point in two directions are the same. We then use the same notation
% to denote the number of grid points and the spacing discretization in
% both directions.
nx = params.grid_size(1);
dx = params.dx;
dt = params.dt1;
u_half = zeros(nx);
Algomat = zeros(nx);

% The spacing on the x direction is j = 1:nx.
% apply second order central difference in y direction to approximate
% second order derivative in y direction
u_pad_y = [u(2, :); u; u(end - 1, :)];
du2dy = u_pad_y(1:end - 2, :) - 2 * u_pad_y(2:end - 1, :) + u_pad_y(3:end, :);

% update y direction explicitly
rhs = u + (dt / 2) * (D * du2dy / dx ^ 2 + f);

% solve x direction implicitly
mu_x = D * dt / dx ^ 2;
for i = 1:nx
    u_half(i, :) = Tridiagonal(rhs(i, :), nx, mu_x);
end

% The spacing on the y direction is i = 1:nx.
% apply second order central difference in x direction to approximate
% second order derivative in x direction
u_pad_x = [u_half(:, 2) u_half u_half(:, end - 1)];
du2dx = u_pad_x(:, 1:end - 2) - 2 * u_pad_x(:, 2:end - 1) + u_pad_x(:, 3:end);

% update x direction explicitly
rhs = u_half + (dt / 2) * (D * du2dx / dx ^ 2 + f);

% solve y direction implicitly
mu_y = D * dt / dx ^ 2;
for j = 1:nx
    Algomat(:, j) = Tridiagonal(rhs(:, j), nx, mu_y);
end
end

function mat = Tridiagonal(b, n, mu, d, a ,c)
% d is the vector for diagonal elements, a is the vector for subdiagonal
% elements, c is the vector for superdiagonal elements, and b is the vector
% for the right-hand side. n is the dimension of the linear tridiagonal
% system. The input b and n are necessary.

if nargin == 3
    d = (1 + mu) * ones(1, n);
    a = zeros(1, n - 1);
    c = zeros(1, n - 1);
    for i = 1:n - 2
        a(i) = - mu / 2;
    end
    a(n - 1) = - mu;
    for i = 2:n - 1
        c(i) = - mu / 2;
    end
    c(1) = - mu;
end

mat = zeros(1, n);
for i = 2:n
    d(i) = d(i) - (a(i - 1) / d(i - 1)) * c(i - 1);
    b(i) = b(i) - (a(i - 1) / d(i - 1)) * b(i - 1);
end

mat(n) = b(n) / d(n);
for i = n - 1:-1:1
    mat(i) = (b(i) - c(i) * mat(i + 1)) / d(i);
end

end