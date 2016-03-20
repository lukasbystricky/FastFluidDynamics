function [u_out, v_out, dsolution] = lid_driven_cavity(N, dt, tf, visc)

%% FFD : Fast Fluid Dynamics
%
% This code solve the lid driven cavity using FFD as described in the paper
% "Stable Fluids" by Jos Stam. 

nsteps = floor(tf/dt);

u = zeros((N-2)^2, 1);
v = zeros((N-2)^2, 1);

u_out = zeros(N,N,nsteps);
v_out = zeros(N,N,nsteps);

h = 1/(N-1);

%prepare matrices for diffusion and projection
N2 = length(u);

e = -1*ones(N2,1);
ep1 = -1*ones(N2,1);
ep1(sqrt(N2)+1:sqrt(N2):end) = 0;

em1 = -1*ones(N2,1);
em1(sqrt(N2):sqrt(N2):end) = 0;

d = (4 + h^2/(dt*visc))*ones(N2,1);

A_diffuse = (dt*visc/h^2)*spdiags([e, em1, d, ep1, e], [-sqrt(N2), -1, 0 1, sqrt(N2)], N2, N2);

N2 = N^2;
ep1 = 1*ones(N2,1);
ep1(sqrt(N2)+1:sqrt(N2):end) = 0;
ep1(2:sqrt(N2):end) = 2;

epn = 1*ones(N2,1);
epn(1:sqrt(N2)+sqrt(N2)) = 2;

em1 = 1*ones(N2,1);
em1(sqrt(N2):sqrt(N2):end) = 0;
em1(sqrt(N2)-1:sqrt(N2):end) = 2;

emn = 1*ones(N2,1);
emn(N2-2*sqrt(N2)+1:N2) = 2;

d = -4*ones(N2,1);

A_project = 1/(h^2)*spdiags([emn, em1, d, ep1, epn], [-sqrt(N2), -1, 0 1, sqrt(N2)], N2, N2);  

A_project(end,:) = 0;
A_project(end,end) = 1;


%enter time loop
t = 0;
tstep = 0;

u_old = u;
v_old = v;
tol = 1e-10;

while t < tf
    t = t+dt;
    tstep = tstep + 1;
    
    [u,v] = diffuse(u, v, h, dt, visc, A_diffuse);
    [u,v] = project(u, v, h, A_project);
    
    [u,v] = advect(u,v,h,dt);
    [u,v] = project(u,v,h, A_project);
    
    dsolution = norm(([u,;v]-[u_old;v_old])/N^2);
    if (dsolution < tol)
        disp(['steady state reached at t=',num2str(t)]);
        break;
    end
    
    u_old = u;
    v_old = v;
    
    [utmp,vtmp] = add_bc(u,v);

    u_out(:,:,tstep) = rot90(utmp);
    v_out(:,:,tstep) = rot90(vtmp);
    
end

end


%% diffusion step - solve diffusion equation using backwards Euler
function [u_new, v_new] = diffuse(u, v, h, dt, visc, A)

N2 = length(u);    
u(1:sqrt(N2)) = u(1:sqrt(N2)) + dt*visc/h^2;

u_new = A\u;
v_new = A\v;

end


%% projection step - solve Poisson equation
function [u_new, v_new] = project(u, v, h, A)

N2 = (sqrt(length(u)) + 2)^2;
[u_w_bc, v_w_bc] = add_bc(u,v);

%calculate div(velocity)
dv = zeros(sqrt(N2),sqrt(N2));

%interior
for i = 2:sqrt(N2)-1
    for j = 2:sqrt(N2)-1      
        dv(i,j) = (u_w_bc(i-1,j) - u_w_bc(i+1,j))/(2*h) + ...
            (v_w_bc(i,j-1)-v_w_bc(i,j+1))/(2*h);
    end
end

%left and right boundaries
for j = 2:sqrt(N2)-1
    dv(1,j) = -u_w_bc(2,j)/h;
    dv(end,j) = u_w_bc(end-1,j)/h;
end

%top and bottom boundaries
for i = 2:sqrt(N2)-1    
    dv(i,1) = -v_w_bc(i,2)/h;
    dv(end,j) =  v_w_bc(i,end-1)/h;
end

dv = reshape(dv, N2, 1);

%solve poisson equation for q
q = A\dv;

%calculate grad(q)
q = reshape(q, sqrt(N2), sqrt(N2));
qx = zeros(sqrt(N2),sqrt(N2));
qy = zeros(sqrt(N2),sqrt(N2));

for i = 2:sqrt(N2) - 1
    for j = 2:sqrt(N2) - 1
       qx(i,j) = (q(i-1,j) - q(i+1,j))/(2*h);
       qy(i,j) = (q(i,j-1) - q(i,j+1))/(2*h);
    end
end

%left and right boundaries, qx = 0 from boundary conditions
for j = 2:sqrt(N2)-1
    qy(1,j) = (q(1,j-1) - q(1,j+1))/(2*h);
    qy(end,j) = (q(end,j-1) - q(end,j+1))/(2*h);
end

%top and bottom boundaries, qy = 0 from boundary conditions
for i = 2:sqrt(N2)-1    
    qx(i,1) = (q(i-1,1) - q(i+1,1))/(2*h);
    qx(i,end) = (q(i-1,end) - q(i+1,1))/(2*h);
end

u_new = u_w_bc - qx;
v_new = v_w_bc - qy;

u_new = reshape(u_new(2:end-1,2:end-1), size(u));
v_new = reshape(v_new(2:end-1,2:end-1), size(v));

end

%% advection step
function [u_new, v_new] = advect(u, v, h, dt)

[u_w_bc, v_w_bc] = add_bc(u,v);

u_new = zeros(size(u_w_bc));
v_new = zeros(size(v_w_bc));

[N,~] = size(u_w_bc);

for i = 2:N-1
    for j = 2:N-1
        
        %find position at t-dt
        x0 = (i-1)*h;
        y0 = (j-1)*h;
        
        x1 = min(max(x0 - dt*u_w_bc(i,j),0),1);
        y1 = min(max(y0 - dt*v_w_bc(i,j),0),1);
        
        %identify cell 
        iL = min(floor(x1/h)+1,N-1);
        iR = iL + 1;
        jB = min(floor(y1/h)+1,N-1);
        jT = jB + 1;
        
        xL = (iL-1)*h;
        xR = (iR-1)*h;
        yB = (jB-1)*h;
        yT = (jT-1)*h;
        
        uTL = u_w_bc(iL,jT);
        uTR = u_w_bc(iR,jT);
        uBL = u_w_bc(iL,jB);
        uBR = u_w_bc(iR,jB);
        
        vTL = v_w_bc(iL,jT);
        vTR = v_w_bc(iR,jT);
        vBL = v_w_bc(iL,jB);
        vBR = v_w_bc(iR,jB);
        
        %bilinear interpolation
        A = [1, xR, yB, xR*yB;
             1, xR, yT, xR*yT;
             1, xL, yB, xL*yB;
             1, xL, yT, xL*yT];
        u_vec = [uBR; uTR; uBL; uTL];
        v_vec = [vBR; vTR; vBL; vTL];
        
        coeffs_u = A\u_vec;
        coeffs_v = A\v_vec;
        
        u_new(i,j) = coeffs_u(1) + coeffs_u(2)*x1 + coeffs_u(3)*y1 + coeffs_u(4)*x1*y1;
        v_new(i,j) = coeffs_v(1) + coeffs_v(2)*x1 + coeffs_v(3)*y1 + coeffs_v(4)*x1*y1;       
    end
end

u_new = reshape(u_new(2:end-1,2:end-1), size(u));
v_new = reshape(v_new(2:end-1,2:end-1), size(v));  

end


function [u_w_bc, v_w_bc] = add_bc(u,v)

N2 = length(u);

u = reshape(u, sqrt(N2), sqrt(N2));
v = reshape(v, sqrt(N2), sqrt(N2));

u_w_bc = zeros(sqrt(N2) + 2, sqrt(N2) + 2);
v_w_bc = zeros(sqrt(N2) + 2, sqrt(N2) + 2);

u_w_bc(:,1) = 1;
u_w_bc(2:end-1,2:end-1) = u;

v_w_bc(2:end-1,2:end-1) = v;

end
