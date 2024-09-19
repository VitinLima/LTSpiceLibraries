R1 = 1e3;
R2 = 10e3;
R3 = 100e3;

Vh = 2.5;
Vl = 0.6;

A = [5-Vh -Vh; 5-Vl -Vl];
B = [Vh/R3+(Vh-5)/(R1+R2); Vl*(1/R2+1/R3)];
1./(A\B)