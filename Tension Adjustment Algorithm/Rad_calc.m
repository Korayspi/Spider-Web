clear; clc;

n_anc=5;
n_rad=30;
L_anc=10;
L_fra=8;
spi_dist=0.2;

Fspi=10e-6;
Frad=100e-6;
Ffr=0e-6;
Fanc=0e-6;

[Nod, I, J, ele_type, nodo_central, ang_rad]=MakeWeb_nospi(n_anc,n_rad,L_anc,L_fra);

n=length(Nod);
m=length(ele_type);
ele_type=ele_type';

%% Thread Parameters
rad_no=sum(ele_type==2);
fr_no=sum(ele_type==3);
anc_no=sum(ele_type==4);

% Diameter of spi, rad, fr, anc
D_thread=[2.4e-6, 3.93e-6, 7.23e-6, 8.03e-6]'; 
D=D_thread(ele_type);

A_ele=pi*(D/2).^2;
A=pi*(D_thread/2).^2;
% 

%% Structure Info 

strain0=zeros(4,1);
stress0=zeros(4,1);

% Major Ampullate
while A(2)*stress0(2) < Frad
    strain0(2) = strain0(2)+1e-5;
    stress0(2) = get_stress_real(strain0(2),2);
end
while A(3)*stress0(3) < Ffr
    strain0(3) = strain0(3)+1e-5;
    stress0(3) = get_stress_real(strain0(3),2);
end
while A(4)*stress0(4) < Fanc
    strain0(4) = strain0(4)+1e-5;
    stress0(4) = get_stress_real(strain0(4),2);
end
% Minor Ampullate
while A(1)*stress0(1) < Fspi
    strain0(1) = strain0(1)+1e-4;
    stress0(1) = get_stress_real(strain0(1),1);
end

F=A.*stress0;
F0_ele=F(ele_type);
    
edge=I(1:n_anc); % Determining attachment points
Fapp=0*Nod;
moving=1:length(Nod);
moving(edge)=[]; % Attechment points don't move

%% Calculation of L0
L_ele=vecnorm(Nod(I,:)-Nod(J,:),2,2);
rads=(ele_type==2);
rad_tension=Frad;
rad_tension=rad_tension';
ele_L0=L_ele;
frs=(ele_type==3);
ele_L0(frs)=ele_L0(frs).*1.02;
%% Maj. Amp. Silk Tension Equalizing
t=0;
%     scl_un=0.9766; % Scale unit
scl_un=15;
dx_un=3.5290;
div_un=1.2;
% div_un=5;

%     scl_un=vars(1); % Scale unit
% dx_un=vars(2);
%     div_un=vars(2);
tic;

dx=dx_un*ones(n_rad,1);

[Fnod, F_ele, Nodm, max_val]=relax_6err(Nod, A_ele, ele_L0, Fapp, moving, I, J, ele_type,1e-7);    
val=F_ele(rads)-rad_tension;
sign_val=0*dx;

%     count=100;
max_err=[];

while sum(abs(val)>5e-7,'All')
%     for i=1:count
    valp=val;
    val=F_ele(rads)-rad_tension;
    ele_L0(rads)=ele_L0(rads)+dx.*val;
    ele_L0(ele_L0<0.0001)=0.0001;

    dx(val.*valp<0)=dx(val.*valp<0)./div_un;
    dx(val.*valp>=0)=dx(val.*valp>=0)+scl_un;

    [Fnod, F_ele, Nodm, max_val]=relax_6err(Nod, A_ele, ele_L0, Fapp, moving, I, J, ele_type,1e-7);    
%         sign_val=sign_val+sign(valp.*val);
%         if mod(t,10)==0
%                 sum(sign_val<-4)
%                 dx(sign_val>4)=dx(sign_val>4)+scl_un;
%                 dx(sign_val<-4)=dx(sign_val<-4)./div_un;
%                 sign_val=0*dx;
%         end
    t=t+1;
%         max_val(t)=max(F_ele);
%         min_val(t)=min(F_ele);
    prev_val=max(abs(val))

end
max_err(end+1)=max_val;
Nod=Nodm;
%     opt_val=max(abs(val));
strain_rads=L_ele(rads)./ele_L0(rads);
save('radial_pretension.mat')
%  
figure()
PlotWeb3D(Nod,I,J);
%     title(['Assymetry=',num2str(assym_val)])
F_elep=round(F_ele*1e6);
% 
% eles=(ele_type==2);
eles=1:m;
ele_pos=(Nod(I(eles),:)+Nod(J(eles),:))/2;
F_plt=F_elep(eles);
for i=1:length(ele_pos)
    text(ele_pos(i,1)-0.5,ele_pos(i,2),ele_pos(i,3),num2str(F_plt(i)),'FontSize', 12)
%         text(ele_pos(i,1),ele_pos(i,2),ele_pos(i,3),num2str(i),'FontSize', 12)
end
view(0,90)
axis off
