funcprot(0)

function y = f(x)
    y = log(x + 1);
endfunction

h = 1e-4;        
n_max = 5;

function d = derivee(f, x0, n, h)
    m = ceil(n/2) + 1;
    kk = -m:m;

    A = zeros(2*m+1, 2*m+1);
    for i = 0:2*m
        A(i+1, :) = kk.^i;
    end

    B = zeros(2*m+1, 1);
    B(n+1) = factorial(n);

    c = A \ B;

    d = 0;
    for i = 1:length(kk)
        d = d + c(i) * f(x0 + kk(i)*h);
    end

    d = d / h^n;
endfunction

function coeff = calc_DL(f, n, h)
    coeff = zeros(1, n+1);
    for k = 0:n
        coeff(k+1) = derivee(f, 0, k, h) / factorial(k);
    end
endfunction

function y = DL(x, coeff)
    y = zeros(x);
    for k = 0:length(coeff)-1
        y = y + coeff(k+1) * (x.^k);
    end
endfunction

coeff = calc_DL(f, n_max, h);

printf("k | DL\n");
for k = 0:n_max
    printf("%d | %.8f\n", k, coeff(k+1));
end

x = linspace(-0.5, 0.5, 100);
y_exact = f(x);

scf();
plot(x, y_exact, 'k');

colors = ['r','g','b','m','y','c'];

for n = 1:n_max
    coeff_n = calc_DL(f, n, h);
    y_DL = DL(x, coeff_n);
    plot(x, y_DL, colors(n));
end

legend('Exact','n=1','n=2','n=3','n=4','n=5');
title('Convergence du developpement limite');

x_test = 0.3;
err = zeros(1, n_max);

for n = 1:n_max
    coeff_n = calc_DL(f, n, h);
    approx = DL(x_test, coeff_n);
    exact = f(x_test);
    err(n) = abs(exact - approx);
end

disp('Erreurs :');
disp(err);

scf();
plot(1:n_max, err, '-o');
xlabel('Ordre n');
ylabel('Erreur');
title('Erreur vs ordre du DL');
// partie integrale
a = 0;
b = 1;
n = 10;
h = (b - a) / n;

I_trap = 0;

for i = 1:n-1
    x = a + i*h;
    I_trap = I_trap + f(x);
end

I_trap = h * ( (f(a) + f(b))/2 + I_trap );

I_simp = f(a) + f(b);

for i = 1:n-1
    x = a + i*h;
    
    if modulo(i,2) == 0 then
        I_simp = I_simp + 2*f(x);
    else
        I_simp = I_simp + 4*f(x);
    end
end

I_simp = I_simp * h / 3;

printf("Trapèzes = %f\n", I_trap);
printf("Simpson  = %f\n", I_simp);

I_exact = (2*log(2) - 1);

printf("Valeur exacte = %f\n", I_exact);

err_trap = abs(I_exact - I_trap);
err_simp = abs(I_exact - I_simp);

printf("Erreur Trapèzes = %e\n", err_trap);
printf("Erreur Simpson  = %e\n", err_simp);

if err_simp < err_trap then
    disp("Simpson est plus précise");
else
    disp("Trapèzes est plus précise");
end
n_values = 10:10:100;
err_trap_tab = [];
err_simp_tab = [];

for n = n_values
    
    h = (b - a)/n;
    
    I_trap = 0;
    for i = 1:n-1
        I_trap = I_trap + f(a + i*h);
    end
    I_trap = h*((f(a)+f(b))/2 + I_trap);
    
    I_simp = f(a) + f(b);
    for i = 1:n-1
        if modulo(i,2)==0 then
            I_simp = I_simp + 2*f(a+i*h);
        else
            I_simp = I_simp + 4*f(a+i*h);
        end
    end
    I_simp = I_simp*h/3;
    
    err_trap_tab($+1) = abs(I_exact - I_trap);
    err_simp_tab($+1) = abs(I_exact - I_simp);
    
end

scf();
loglog(n_values, err_trap_tab, '-o');
loglog(n_values, err_simp_tab, '-x');
legend("Trapèzes", "Simpson");
xlabel("n");
ylabel("Erreur");
title("Convergence des méthodes");
