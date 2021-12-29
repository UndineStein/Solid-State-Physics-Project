clc
clear all
close all

hbar = 1.054571818e-34;
me = 9.10938356e-31;
a = 1e-10;
pi = 3.14159265359;
epsilon0 = 8.854188e-12; 
elm = 9.109534e-31;
qel = 1.602189e-19;
rydberg = 13.606;
ekinunit = 3.8100;
E_0 = 1;
recipunit = 1.0e+10;
lambda = 12.264;
set_pot = 14;
nband=15;
tol = 10e-14;

[bz,msg] = fopen("bzpath.dat");
nq = fskipl(bz,Inf)-1;
frewind(bz);
fskipl(bz,1);
is = zeros(nq,1);
q=zeros(nq,5);
for iq = 1:nq
  data = sscanf(fgetl(bz),"%g");
  is(iq,1) = data(2);
  q(iq, 1:5) = data(3:7);
endfor
fclose(bz);
%segments = max(is(:,1));


printf("\n Table of wavevectors read in input file bzpath.dat \n")
printf("  iq            is            q(1)            q(2)            q(3)            q(4)            q(5) \n")
for iq = 1:nq
  printf("%4d %15.6G %15.6G %15.6G %15.6G %15.6G %15.6G \n",iq,is(iq),q(iq,:))
endfor


printf("\nDimensionless lattice dependent part")
cutoff = 21;
Gs_max = 11;

printf("\nfcc lattice unit vectors in cartesian coordinates \n")
a = zeros(3,3);
a(:,1) = [0.5 0.5 0.0]' ;
a(:,2) = [0.0 0.5 0.5]' ;
a(:,3) = [0.5 0.0 0.5]' ;

printf("  a_1 a_2 a_3\n")
disp(a)
cell_volume = a(:,1)' * cross(a(:,2),a(:,3))

printf("\nfcc reciprocal lattice unite vectors in Cartesian coordinates");
g=zeros(4,3);
g(1:3,1) = cross(a(:,2),a(:,3)) /cell_volume ;
g(1:3,2) = cross(a(:,3),a(:,1)) /cell_volume ;
g(1:3,3) = cross(a(:,1),a(:,2)) /cell_volume ;
printf("\n g_1  g_2   g_3\n");

for i = 1:3
  g(4,i) = g(1:3,i)'*g(1:3,i);
endfor
disp(g)

min_norm = sqrt(min(g(4,:)));
nstep = floor(sqrt(cutoff)/min_norm)+1;
printf("\n cutoff requires %d positive steps along each recip lattice vector \n", nstep);
nodes = (2*nstep+1)^3;
printf("Generate (2 * %1d +1)^3 =%4d reciprocal lattice vectors \n", nstep, nodes);
G = zeros(5,nodes);
n = 0;
for j =-nstep:nstep
  for k=-nstep:nstep
    for l=-nstep:nstep
      n++;
      G(1:3,n)=j*g(1:3,1)+k*g(1:3,2)+l*g(1:3,3);
      G(5,n) = G(1:3,n)' * G(1:3,n);
      G(4,n) = sqrt(G(5,n));
    endfor
  endfor
endfor

GT = sortrows(G',4);

G = GT';

kept = 0;
for i = 1:nodes
  if(G(5,i)) < 21;
  kept = kept + 1;
endif
endfor
disp(kept)
%{
for i=1:kept
printf(" %4d %15.6G %15.6G %15.6G %15.6G %15.6G \n",i,G(:,i)) 
endfor
%}

%          V0    VS3  VS8  VS11 VA3  VA4  VA11  
ff(1,:) = [0.00 -0.21 0.04 0.08 0.00 0.00 0.00];  %Si
ff(2,:) = [0.00 -0.23 0.01 0.06 0.00 0.00 0.00]; %Ge
ff(3,:) = [0.00 -0.20 0.00 0.04 0.00 0.00 0.00]; %Sn
ff(4,:) = [0.00 -0.22 0.03 0.07 0.12 0.07 0.02]; %GaP
ff(5,:) = [0.00 -0.23 0.01 0.06 0.06 0.05 0.01]; %GaAs
ff(6,:) = [0.00 -0.21 0.02 0.06 0.06 0.04 0.02]; %AlSb
ff(7,:) = [0.00 -0.23 0.01 0.06 0.07 0.05 0.01]; %InP
ff(8,:) = [0.00 -0.22 0.00 0.05 0.06 0.05 0.01]; %GaSb
ff(9,:) = [0.00 -0.22 0.00 0.05 0.08 0.05 0.03]; %InAs
ff(10,:)= [0.00 -0.20 0.00 0.04 0.06 0.05 0.01]; %InSb
ff(11,:)= [0.00 -0.22 0.03 0.07 0.24 0.14 0.04]; %ZnS
ff(12,:)= [0.00 -0.23 0.01 0.06 0.18 0.12 0.03]; %ZnSe
ff(13,:)= [0.00 -0.22 0.00 0.05 0.13 0.10 0.01]; %ZnTe
ff(14,:)= [0.00 -0.20 0.00 0.04 0.15 0.09 0.04]; %CdTe

ls = [5.43 5.66 6.49 5.44 5.64 6.13 5.86 6.12 6.04 6.48 5.41 5.65 6.07 6.41]
l = ls(set_pot);

tau = zeros(3,2);
tau(:,1) = [0.125 0.125 0.125]' ;
tau(:,2) = [-0.125 -0.125 -0.125]' ;

printf("Energy in eV\n")
for n=1:kept
  sym = 0.;
  asym = 0.;
  if(G(5,n)==3)
    sym = ff(set_pot,2)*rydberg;
    asym = ff(set_pot, 5)*rydberg;
  endif
  if(G(5,n)==4)
    sym = 0.;
    asym = ff(set_pot,6)*rydberg;
  endif
  if(G(5,n)==8)
    sym = ff(set_pot,3)*rydberg;
    asym = 0.;
   endif
   if(G(5,n)==11)
    sym = ff(set_pot, 4)*rydberg;
    asym = ff(set_pot,7)*rydberg;
  endif
  argu = 2*pi*(G(1:3,n)'*tau(1:3,1));
  cvg(n) = cos(argu)*sym-(0+1i)*sin(argu)*asym;
endfor

printf("\nTABLE OF G VECTORS (unit of 2*pi/lattice spacing) AND FOURIER COEFFICIENTS OF PSEUDOPOTENTIAL (eV) \n")
printf("     n              G1            G2            G3            |G|           |G|^2           Re(V_G)           Im(V_G) \n")
for n = 1:kept
  printf("%4d %15.6G %15.6G %15.6G %15.6G %15.6G %15.6G %15.6G \n",n,G(1:5,n),real(cvg(n)), imag(cvg(n)));
endfor

f1 = fopen('bz3.dat',"w");
printf("DIAGONOLIZATION LOOP OVER %4d WAVEVECTORS :", nq); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
G_diff = zeros(5,1);

H = zeros(kept,kept);

%Potential energy
for j = 1:kept
  for i = 1:kept
     G_diff(1:3) = G(1:3,i)-G(1:3,j);
     G_diff(5) = G_diff(1:3)'*G_diff(1:3);
      if(G_diff(5)<=Gs_max)
        for k = 1:kept
          if((G_diff(1:3)-G(1:3,k))==[0 0 0]')
            H(i,j) = cvg(k);
          endif
        endfor
      endif
  endfor
endfor
disp(G_diff)
%Kinetic energy
if(kept<nband)
   nband=kept;
endif

for iq = 1:nq
  printf("%4d",iq);   
  for i = 1:kept                         
    for j = 1:3
      p(j) = q(iq, j)-G(j,i);  %q index is the k vector only on diagonal
    endfor
  H(i,i) = (ekinunit*(2*pi/l)^2*(p*p'));    
  endfor
  if(any(H-H'))
    printf("\nHamiltonian is not Hermitian : fatal error.\n");
    break;
  else
    [v,ev]=eig(H);
    E = real(diag(ev));
    [E perm]=sort(E);
    v = v(:,perm);
    bzs = 1000*(-1)^is(iq,1);

      fprintf(f1,"%15.6G %15.6G",q(iq,5),bzs);
     for i=1:nband
       fprintf(f1,"%15.6G",E(i)); 
     endfor
        fprintf(f1,"\n");
      
        if(q(iq,4)==0)
          #printf("\n gamma.eigenvalues")
            for i=1:nband
              gamma(i)=E(i);
              #printf("%d %15.6G \n",i, gamma(i));
            endfor
         endif
  endif
endfor
fclose(f1);
 printf("\n gamma.eigenvalues \n")
for i=1:nband
      printf("%d %15.6G \n",i, gamma(i));
endfor































 