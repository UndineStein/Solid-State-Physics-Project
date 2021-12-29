clc
clear all
close all
format long
more off
labelx = 'Wavevector 2\pi/a unites'
labely = 'Energy (eV)'
datafile = "bz3.dat"

printf("\nMin and max abscisse set by user:  \n");
x_max = 2.5;
x_min = 0;
printf("\nMin and max ordinates set by user:  \n");
y_min =-5;
y_max = 16;

%Reading input data
printf("\n");
[inputfile,msg] = fopen(datafile,"r");

if(inputfile ==1)
  printf("File \"%s\" does not exist. \n", datafile);
  return
endif

abscisse=fskipl(inputfile,Inf);
frewind(inputfile);
str = fgetl(inputfile);

columns = rows(sscanf(str,"%g"))
frewind(inputfile);
x = zeros(abscisse, columns);
for i=1:abscisse
  x(i,1:columns) = sscanf(fgetl(inputfile), "%g");
endfor
fclose(inputfile);

%Close
printf("\nMin ana max abcissae found in data:  \n")
data_x_min = min(x(:,1))
data_x_max = max(x(:,1))

bmin=min(x(:,2:end));
bmax=max(x(:,2:end));

printf("\nMin and max ordinates in data:  \n");
data_y_min = min(bmin);
data_y_max = max(bmax);
printf("\n");


%Plot
h = plot(x(:,1),x(:,2:end),"linewidth",2)
%Next lines must come after plot

xlim([x_min,x_max]);
ylim([y_min,y_max]);

xlabel(labelx, "interpreter", "latex");
ylabel(labely, "interpreter", "latex");

ax = gca();
set(ax,"linewidth",2)
set(ax, "fontweight", "bold")
perso_colororder = [
0 0 0 
1 0 0 
0 0 1
0 0.5 0
0.15 0.45 0.45
0 0.45 0.45
0.5 0.5 0.5
];
set(ax, "colorirder", perso_colororder);

%output file name
dot = rindex(datafile,".");
root = substr(datafilr,1,dot-1);
output_file = sprintf("%s.pdf",root)

%Basic pdf output
print(output_file, "")
