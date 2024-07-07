clear all
close all
clc
image1 = imread('image0.png');
image2 = imread('image1.png');

subplot(1, 2, 1);
imshow(image1);
title('Line data and cost coefficients');

subplot(1, 2, 2);
imshow(image2);
title('Bus data');
disp('                       Power System III Lab')
disp('                         EEE 4732 Project')
disp('Optimization Methods for Economic Dispatch of Thermal Generating Units')
disp('                         ')
disp('Rayan Hossain Khan     190021214')
disp('Shadman Saqib          190021228')
disp('Md. Sifat Aziz         190021131')
disp('                         ')

a = 'We are considering the "IEEE 14 BUS SYSTEM" with 3 generators';
a = [a newline 'The Fuel Cost Functions for the three generators are given by:'];
a = [a newline 'C1 = 105 + 2.45*P1 + 0.005*(P1)^2'];
a = [a newline 'C2 = 44.10 + 3.51*P2 + 0.005*(P2)^2'];
a = [a newline 'C3 = 40.60 + 3.89*P3 + 0.005*(P3)^2'];
a = [a newline];
a = [a newline 'The generation limits of the three generators are given by:'];
a = [a newline '10 <= P1 <= 160'];
a = [a newline '20 <= P2 <= 80'];
a = [a newline '20 <= P3 <= 50'];
a = [a newline];
a = [a newline 'The load demand is taken as PD = 250 W'];
disp(a)
disp('                         ')

chr = 'Select 1 for ELD without limits and losses';
chr = [chr newline 'Select 2 for ELD with limits ignoring losses'];
chr = [chr newline 'Select 3 for ELD with limits and losses'];
chr = [chr newline];
ib=input(chr);

if ib==1
    lambda = (250 + 2.45/(2*0.005) + 3.51/(2*0.005) + 3.89/(2*0.005))/(1/(2*0.005) + 1/(2*0.005) + 1/(2*0.005));
    P1 = (lambda - 2.45)/(2*0.005);
    P2 = (lambda - 3.51)/(2*0.005);
    P3 = (lambda - 3.89)/(2*0.005);
    fprintf('\nIncremental cost, Lambda = %.2f $/MWh \n',lambda)
    fprintf('The optimal generation is: \nP1 = %.2f MW \nP2 = %.2f MW \nP3 = %.2f MW \n', P1, P2, P3)
    CT = 105 + 2.45*P1 + 0.005*(P1)^2 + 44.10 + 3.51*P2 + 0.005*(P2)^2 + 40.60 + 3.89*P3 + 0.005*(P3)^2;
    fprintf('Total generated power = %.2f\n',P1+P2+P3)
    fprintf('The total cost is %.2f $/h\n', CT)
elseif ib==2
    b2
elseif ib==3
    b3
end
























