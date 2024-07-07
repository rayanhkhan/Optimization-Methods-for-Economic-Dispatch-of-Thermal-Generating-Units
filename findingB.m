basemva = 100;  
accuracy = 0.0001;  
maxiter = 100;

%-----------------------------------------busdata--------------------------------------------------

%        Bus Bus  Voltage Angle   ----Load---- ----------Generator----------   Injected
%        No  code Mag.    Degree  MW     Mvar     MW      Mvar   Qmin  Qmax    Mvar
busdata=[1   2    1.06    0.0     0       0       114.17 -16.9   0     10      0
         2   2    1.045   0.0     21.7    12.7    40.0    0     -42.0  50      0
         3   0    1.01    0.0     94.2    19.1    0       0      23.4  40      0
         4   0    1.00    0.0     47.8    -3.9    0       0      0     0       0
         5   0    1.00    0.0     7.6     1.6     0       0      0     0       0
         6   0    1.00    0.0     11.2    7.5     0       0      0     0       0            
         7   1    1.00    0.0     0       0       0       0      0     0       0       
         8   0    1.00    0.0     0       0       0       0      0     0       0      
         9   0    1.00    0.0     29.5    16.6    0       0      0     0       0.19             
         10  0    1.00    0.0     9.0     5.8     0       0      0     0       0           
         11  0    1.00    0.0     3.5     1.8     0       0      0     0       0            
         12  0    1.00    0.0     6.1     1.6     0       0      0     0       0                                                                        
         13  0    1.00    0.0     13.8    5.8     0       0      0     0       0             
         14  0    1.00    0.0     14.9    5.0     0       0      0     0       0];

     
%-----------------------------------------linedata--------------------------------------------------

%                                               Line code
%         Bus bus    R        X       1/2 B     = 1 for lines
%         nl  nr    p.u.     p.u.      p.u.     > 1 or < 1 tr. tap at bus nl
linedata=[1   2   0.01938   0.05917   0.02640     1
          1   5   0.05403   0.22304   0.02190     1
          2   3   0.04699   0.19797   0.01870     1
          2   4   0.05811   0.17632   0.02460     1
          2   5   0.05695   0.17388   0.01700     1
          3   4   0.06701   0.17103   0.01730     1
          4   5   0.01335   0.04211   0.00640     1
          4   7   0         0.20912   0           0.978
          4   9   0         0.55618   0           0.969
          5   6   0         0.25202   0           0.932
          6   11  0.09498   0.1989    0           1
          6   12  0.12291   0.25581   0           1
          6   13  0.06615   0.13027   0           1
          7   8   0         0.17615   0           1
          7   9   0         0.11001   0           1
          9   10  0.03181   0.0845    0           1
          9   14  0.12711   0.27038   0           1
          10  11  0.08205   0.19207   0           1
          12  13  0.22092   0.19988   0           1
          13  14  0.17093   0.34802   0           1];
      
      
%------------------------------------Formation of Bus admittance matrix--------------------------------------------------

nl = linedata(:,1);            %from BUS
nr = linedata(:,2);            %to BUS
R = linedata(:,3);             %Line resistance
X = linedata(:,4);             %Line Reactance
Bc = 1i*linedata(:,5);         %Half line charging susceptance
a = linedata(:, 6);            %tap setting value
nbr = length(linedata(:,1));   %number of lines
nbus = max(max(nl), max(nr));  %number of bus
Z = R + 1i*X;                  %Branch impedance
y = ones(nbr,1)./Z;            %Branch admittance

for n = 1:nbr
    if a(n) <= 0
        a(n) = 1;
    else
    end
    Ybus = zeros(nbus,nbus);     %initialize Ybus to zero
    
%Formation of the off-diagonal elements
    for k = 1:nbr
       Ybus(nl(k),nr(k)) = Ybus(nl(k),nr(k)) - y(k)/a(k);
       Ybus(nr(k),nl(k)) = Ybus(nl(k),nr(k));
    end
end

%Formation of the diagonal elements
for  n = 1:nbus
     for k = 1:nbr
         if nl(k) == n
            Ybus(n,n) = Ybus(n,n) + y(k)/(a(k)^2) + Bc(k);
         elseif nr(k) == n
            Ybus(n,n) = Ybus(n,n) + y(k) + Bc(k);
         else
         end
     end
end

% The Bus admittance, Ybus is formed in a 14x14 matrix


%------------------------------------Power-flow solution--------------------------------------------------

ns = 0; 
ng = 0; 
Vm = 0; 
delta = 0; 
yload = 0; 
deltad = 0;
nbus = length(busdata(:,1));

for k = 1:nbus
    n = busdata(k,1);           %Bus no.
    kb(n) = busdata(k,2);       %Bus code
    Vm(n) = busdata(k,3);       %Voltage Magnitude
    delta(n) = busdata(k, 4);   %Angle
    Pd(n) = busdata(k,5);       %Load MW
    Qd(n) = busdata(k,6);       %Load Mvar
    Pg(n) = busdata(k,7);       %Gen MW
    Qg(n) = busdata(k,8);       %Gen Mvar
    Qmin(n) = busdata(k, 9);    %Qmin
    Qmax(n) = busdata(k, 10);   %Qmax
    Qsh(n) = busdata(k, 11);    %Injected Mvar
    
    delta(n) = pi/180*delta(n);                        %converting to degrees
    V(n) = Vm(n)*(cos(delta(n)) + 1i*sin(delta(n)));   %complex conjugate
    P(n)=(Pg(n)-Pd(n))/basemva;                        %converting real power to p.u.
    Q(n)=(Qg(n)-Qd(n)+ Qsh(n))/basemva;                %converting reactive power to p.u.
    S(n) = P(n) + 1i*Q(n);                             %Apparent Power
end

for k=1:nbus
    if kb(k) == 1
        ns = ns+1;
    else
    end
    if kb(k) == 2 
        ng = ng+1; 
    else
    end
    ngs(k) = ng;
    nss(k) = ns;
end

Ym = abs(Ybus); 
t = angle(Ybus);
m = 2*nbus - ng - 2*ns;
maxerror = 1; 
converge = 1;
iter = 0;

%Newton Raphson iteration
while maxerror >= accuracy && iter <= maxiter              % Test for max. power mismatch
    for i = 1:m
        for k = 1:m
           A(i,k)=0;      %Initializing Jacobian matrix
        end
    end
    iter = iter+1;
    for n = 1:nbus
        nn = n - nss(n);
        lm = nbus + n-ngs(n) - nss(n) - ns;
        J11=0;
        J22=0;
        J33=0;
        J44=0;
           for i = 1:nbr
             if nl(i) == n || nr(i) == n
                if nl(i) == n
                    l = nr(i);
                end
                if nr(i) == n
                    l = nl(i); 
                end
                J11 = J11 + Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l) - delta(n) + delta(l));
                J33 = J33 + Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l) - delta(n) + delta(l));
                if kb(n) ~= 1
                    J22 = J22 + Vm(l)*Ym(n,l)*cos(t(n,l) - delta(n) + delta(l));
                    J44 = J44 + Vm(l)*Ym(n,l)*sin(t(n,l) - delta(n) + delta(l));
                else
                end
                if kb(n) ~= 1  && kb(l) ~= 1
                    lk = nbus + l - ngs(l) - nss(l) - ns;
                    ll = l - nss(l);
                                                %off-diagonal elements of J1
                    A(nn, ll) = -Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l) - delta(n) + delta(l));
                    if kb(l) == 0               %off-diagonal elements of J2
                        A(nn, lk) = Vm(n)*Ym(n,l)*cos(t(n,l) - delta(n) + delta(l));
                    end
                    if kb(n) == 0               %off-diagonal elements of J3
                        A(lm, ll) = -Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l) - delta(n)+delta(l)); 
                    end
                    if kb(n) == 0 && kb(l) == 0  %off-diagonal elements of J4
                        A(lm, lk) = -Vm(n)*Ym(n,l)*sin(t(n,l) - delta(n) + delta(l));
                    end
                else
                end
             else
             end
           end
           Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;
           Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
           if kb(n) == 1                         %Swing bus P
               P(n) = Pk; 
               Q(n) = Qk; 
           end   
             if kb(n) == 2  
                 Q(n) = Qk;
                 if Qmax(n) ~= 0
                     Qgc = Q(n)*basemva + Qd(n) - Qsh(n);
                   if iter <= 7                          % Between the 2th & 6th iterations
                      if iter > 2                        % the Mvar of generator buses are
                        if Qgc  < Qmin(n)                % tested. If not within limits Vm(n)
                            Vm(n) = Vm(n) + 0.01;        % is changed in steps of 0.01 pu to
                        elseif Qgc  > Qmax(n)            % bring the generator Mvar within
                            Vm(n) = Vm(n) - 0.01;        % the specified limits.
                        end 
                      else
                      end
                   else
                   end
                 else
                 end
             end
           if kb(n) ~= 1
             A(nn,nn) = J11;                                  %diagonal elements of J1
             DC(nn) = P(n) - Pk;
           end
           if kb(n) == 0
             A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n)) + J22;    %diagonal elements of J2
             A(lm,nn) = J33;                                  %diagonal elements of J3
             A(lm,lm) = -2*Vm(n)*Ym(n,n)*sin(t(n,n)) - J44;   %diagonal of elements of J4
             DC(lm) = Q(n) - Qk;
           end
    end
    DX = A\DC';
    for n = 1:nbus
          nn = n - nss(n);
          lm = nbus + n - ngs(n) - nss(n) - ns;
          if kb(n) ~= 1
            delta(n) = delta(n) +DX(nn);
          end
          if kb(n) == 0
            Vm(n) = Vm(n) + DX(lm);
          end
     end
    maxerror = max(abs(DC));
    if iter == maxiter && maxerror > accuracy 
       fprintf('\nWARNING: Iterative solution did not converged after ')
       fprintf('%g', iter), fprintf(' iterations.\n\n')
       fprintf('Press Enter to terminate the iterations and print the results \n')
       converge = 0; 
       pause
    else
    end
end


if converge ~= 1
   tech= ('                      ITERATIVE SOLUTION DID NOT CONVERGE'); 
else
   tech=('                   Power Flow Solution by Newton-Raphson Method');
end   
V = Vm.*cos(delta) + 1i*Vm.*sin(delta);
deltad = 180/pi*delta;
k  = 0;


for n = 1:nbus
     if kb(n) == 1
         k = k+1;
         S(n) = P(n) + 1i*Q(n);
         Pg(n) = P(n)*basemva + Pd(n);
         Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
         Pgg(k) = Pg(n);
         Qgg(k) = Qg(n);     
     elseif kb(n) ==2
         k = k+1;
         S(n) = P(n) + 1i*Q(n);
         Qg(n) = Q(n)*basemva + Qd(n) - Qsh(n);
         Pgg(k) = Pg(n);
         Qgg(k) = Qg(n);  
     end
     yload(n) = (Pd(n)- 1i*Qd(n) + 1i*Qsh(n))/(basemva*Vm(n)^2);
end


busdata(:,3) = Vm'; 
busdata(:,4) = deltad';
Pgt = sum(Pg);  
Qgt = sum(Qg); 
Pdt = sum(Pd); 
Qdt = sum(Qd); 
Qsht = sum(Qsh);

%------------------------------------Forming the new bus--------------------------------------------------

disp(tech)
fprintf('                      Maximum Power Mismatch = %g \n', maxerror)
fprintf('                             No. of Iterations = %g \n\n', iter)
head =['    Bus  Voltage  Angle    ------Load------    ---Generation---   Injected'
       '    No.  Mag.     Degree     MW       Mvar       MW       Mvar       Mvar '
       '                                                                          '];
disp(head)
for n=1:nbus
     fprintf(' %5g', n)
     fprintf(' %7.3f', Vm(n))
     fprintf(' %8.3f', deltad(n))
     fprintf(' %9.3f', Pd(n))
     fprintf(' %9.3f', Qd(n))
     fprintf(' %9.3f', Pg(n))
     fprintf(' %9.3f ', Qg(n))
     fprintf(' %8.3f\n', Qsh(n))
end
    fprintf('      \n')
    fprintf('    Total              ')
    fprintf(' %9.3f', Pdt)
    fprintf(' %9.3f', Qdt)
    fprintf(' %9.3f', Pgt)
    fprintf(' %9.3f', Qgt)
    fprintf(' %9.3f\n\n', Qsht)

%------------------------------------Finding the Loss(B) Coefficients--------------------------------------------------

Zbus = inv(Ybus);
ngg = 0;
I = -1/basemva*(Pd - 1i*Qd)./conj(V); %new
ID = sum(I);  %new

for k = 1:nbus
  if kb(k) == 0
  else
      ngg = ngg + 1;  
  end
   if kb(k) == 1 
       ks = k; 
   else
   end
end

d1=I/ID;
DD = sum(d1.*Zbus(ks,:));  %new
kg=0; 
kd=0;

for k = 1:nbus
    if kb(k) ~= 0
        kg = kg + 1;
        t1(kg) = Zbus(ks,k)/DD;   %new
    else
        kd = kd+1;
        d(kd) = I(k)/ID;
    end
end

nd = nbus - ngg;
C1g = zeros(nbus, ngg);
kg = 0;

for k = 1:nbus
  if kb(k) ~= 0
      kg = kg+1;
    for m = 1:ngg
       if kb(m) ~= 0
           C1g(k, kg) = 1;
       else
       end
    end
  else
  end
end

C1gg = eye(ngg,ngg);
C1D = zeros(ngg,1);
C1 = [C1g,conj(d1)'];
C2gD = [C1gg; -t1];
CnD = [C1D;-t1(1)];
C2 = [C2gD,CnD];
C = C1*C2;
kg = 0;

for k = 1:nbus
  if kb(k)~=0
    kg=kg+1;
    al(kg) = (1 - 1i*((Qg(k) + Qsh(k))/Pg(k)))/conj(V(k));  %new
  else
  end
end

alp = [al, -V(ks)/Zbus(ks,ks)];

for k = 1:ngg + 1
    for m = 1:ngg + 1
        if k == m
            alph(k,k) = alp(k);
        else
            alph(k,m) = 0;
        end
    end
end

T = alph*conj(C)'*real(Zbus)*conj(C)*conj(alph);
BB = 0.5*(T + conj(T));

for k = 1:ngg
    for m=1:ngg
        B(k,m) = BB(k,m);
    end 
    B0(k) = 2*BB(ngg + 1, k);
end

B00 = BB(ngg + 1, ngg + 1);
B, B0, B00
PL = Pgg*(B/basemva)*Pgg'+B0*Pgg' + B00*basemva;
%fprintf('Total system loss = %g MW \n', PL)

%--------------------------------------------------------------------------------------------




































