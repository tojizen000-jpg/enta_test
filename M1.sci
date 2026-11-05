function y = f(x)
    y = log(x + 1);
endfunction

x_test = 0.5;
erreurs = zeros(1,5);

// Q1
for n = 1:5
    DL = 0;
    for k = 1:n
        DL = DL + (-1)^(k+1)/k * x_test^k;
    end
    erreurs(n) = abs(f(x_test) - DL);
    printf("n=%d | coeff a%d = %10.6f\n", n, n, (-1)^(n+1)/n);
end

// Q2
printf("\nCoefficients du DL :\n");
for n = 1:5
    printf("a%d = %.6f\n", n, (-1)^(n+1)/n);
end

// Q3
printf("\nErreur |f(x) - DL(x)| a x=0.5 :\n");
for n = 1:5
    printf("n=%d | erreur = %.2e\n", n, erreurs(n));
end

// Q4
figure(1);
plot(1:5, erreurs, 'bo-');
xlabel('n');
ylabel('Erreur');
title('Erreur DL de ln(x+1) a x=0.5');
xgrid();

a = 0; b = 1;
exacte = 2*log(2) - 1;

function I = trapezes(a, b, n)
    h = (b-a)/n;
    x = linspace(a, b, n+1);
    y = log(x+1);
    I = h * (y(1)/2 + sum(y(2:n)) + y(n+1)/2);
endfunction

function I = simpson(a, b, n)
    if modulo(n,2) <> 0 then n = n+1; end
    h = (b-a)/n;
    x = linspace(a, b, n+1);
    y = log(x+1);
    I = (h/3) * (y(1) + 4*sum(y(2:2:n)) + 2*sum(y(3:2:n-1)) + y(n+1));
endfunction

// Q1
eT10 = abs(exacte - trapezes(a,b,10));
eS10 = abs(exacte - simpson(a,b,10));
printf("\nTrapezes (n=10) = %.8f\n", trapezes(a,b,10));
printf("Simpson  (n=10) = %.8f\n", simpson(a,b,10));

// Q2
printf("\nErreur Trapezes = %.2e\n", eT10);
printf("Erreur Simpson  = %.2e\n", eS10);
if eS10 < eT10 then
    printf("=> Simpson est plus precis\n");
end

n_vals = 2:2:50;
eT = zeros(1, length(n_vals));
eS = zeros(1, length(n_vals));
for i = 1:length(n_vals)
    eT(i) = abs(exacte - trapezes(a, b, n_vals(i)));
    eS(i) = abs(exacte - simpson(a,  b, n_vals(i)));
end

// Q3
figure(2);
plot2d("ll", n_vals, eT, style=2);
plot2d("ll", n_vals, eS, style=5);
legend(['Trapezes'; 'Simpson'], 2);
xlabel('n (log)');
ylabel('Erreur (log)');
title('Log-log erreurs integration');
xgrid();

pT = (log(eT($)) - log(eT(1))) / (log(n_vals($)) - log(n_vals(1)));
pS = (log(eS($)) - log(eS(1))) / (log(n_vals($)) - log(n_vals(1)));
printf("\nPente Trapezes = %.2f\n", pT);
printf("Pente Simpson  = %.2f\n", pS);

// Q4
printf("\nn    | Err Trapezes | Err Simpson\n");
for n = 10:10:100
    printf("%-5d| %-13.2e| %.2e\n", n, abs(exacte-trapezes(a,b,n)), abs(exacte-simpson(a,b,n)));
end
