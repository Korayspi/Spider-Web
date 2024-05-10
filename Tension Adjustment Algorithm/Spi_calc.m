clear; clc;

A=load('radial_pretension.mat');
strain_rads=A.strain_rads;
nodo_central=A.nodo_central;
Nod=A.Nod;
I=A.I;
J=A.J;
ele_type=A.ele_type;
ele_L0=A.ele_L0;
ele_L0rad=ele_L0(ele_type==2);
ele_L0fr=ele_L0(ele_type==3);
ele_L0anc=ele_L0(ele_type==4);
pre_Li_rads=vecnorm(Nod(I(ele_type==2),:)-Nod(J(ele_type==2),:),2,2);

n_anc=5;
n_rad=30;
L_anc=10;
L_fra=8;
spi_dist=0.2;

hub_pos=0;
rmv_spi=4;

[Nod, I, J, ele_type, spi_Nod_no]=MakeWeb_add_spis(Nod,I,J,ele_type,n_rad,spi_dist,rmv_spi,nodo_central);



n=length(Nod);
m=length(ele_type);
ele_type=ele_type';

%% Thread Parameters
spi_no=sum(ele_type==1);
rad_no=sum(ele_type==2);
fr_no=sum(ele_type==3);

% Diameter of spi, rad, fr, anc
D_thread=[2.4e-6, 3.93e-6, 7.23e-6, 8.03e-6]'; 

Fspi=10e-6;
D=D_thread(ele_type);

A_ele=pi*(D/2).^2;
A=pi*(D_thread/2).^2;
% 

%% Structure Info 

% Minor Ampullate
stress0=0;
strain0=0;
while A(1)*stress0 < Fspi
    strain0 = strain0+1e-4;
    stress0 = get_stress_real(strain0,1);
end

% F_elem=F0_ele;

edge=I(1:n_anc); % Determining attachment points
Fapp=0*Nod;
moving=1:length(Nod);
moving(edge)=[]; % Attechment points don't move

%% Calculation of L0
L_ele=vecnorm(Nod(I,:)-Nod(J,:),2,2);
rads=(ele_type==2);
fr=(ele_type==3);
anc=(ele_type==4);

L_rads=L_ele(rads);

t=1;
spi_Nod_no=spi_Nod_no-rmv_spi+1;

for i=1:n_rad
    L_rads(t:t+spi_Nod_no(i)-1)=L_rads(t:t+spi_Nod_no(i)-1)*(ele_L0rad(i)/pre_Li_rads(i));
    t=t+spi_Nod_no(i);
end
t=1;
edge_rad=zeros(n_rad,1);
center_rad=zeros(n_rad,1);
for i=1:n_rad
    edge_rad(i)=t+spi_Nod_no(i)-1;
    center_rad(i)=t;
    t=t+spi_Nod_no(i);
end

ele_L0=L_ele;

ele_L0(rads)=L_rads;
ele_L0(fr)=ele_L0fr;
ele_L0(anc)=ele_L0anc;
%% Radial relaxation
[Fnod, F_ele, Nodm, max_val]=relax_6err_spis(Nod, A_ele, ele_L0, Fapp, moving, I, J, ele_type, 1e-7);


%% Min. Amp. Silk Tension Equalizing
spis=(ele_type==1);
L_spis=L_ele(spis);

t=0;
scl_un=1e-4; % Scale unit
dx_un=2e-1;
div_un=2;

dx=dx_un*ones(spi_no,1).*ele_L0(spis); 
val=F_ele(spis)-Fspi;
sign_val=0*dx;

count=300;
max_val=[];
min_val=[];

%% %%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%% %%
% while sum(abs(val)>1e-7,'All')
for i=1:count
    t=t+1
    valp=val;
    val=F_ele(spis)-Fspi;
    direction=sign(val);
    ele_L0(spis)=ele_L0(spis)+dx.*direction;
    ele_L0(ele_L0<0.0001)=0.0001;
    dx(val.*valp<0)=dx(val.*valp<0)./div_un;
    dx(val.*valp>0)=dx(val.*valp>0)+scl_un;
    
    [~, F_ele, ~, ~]=relax_6err_spis(Nod, A_ele, ele_L0, Fapp, moving, I, J, ele_type,1e-7);        
 
    max_val(t)=max(abs(val));
    if max_val(t)<1.0e-6
        break
    end
end
[Fnod, F_ele, Nodm, F_val_nod]=relax_6err_spis(Nod, A_ele, ele_L0, Fapp, moving, I, J, ele_type,1e-7);
Nod=Nodm;

%% FE divide and final relaxation
ele_no_crs=length(I); %Length 
I_crs=I;
J_crs=J;

ele_type0=ele_type;
rads=(ele_type==2);
fr=(ele_type==3);
anc=(ele_type==4);
Irads=I(rads);
Jrads=J(rads);
Ifr=I(fr);
Jfr=J(fr);
Ianc=I(anc);
Janc=J(anc);

Ispis=I(spis);
Jspis=J(spis);
fe_len=0.1;
moving=1:length(Nod);
moving(nodo_central)=[];
Nod_abaq=Nod;
Nod_abaq(:,3)=Nod_abaq(:,3)+(-1+2*rand(length(Nod),1)); % This is for 3D relaxation of the web
[Fnod, F_ele, Nod_abaq, dx]=relax_zero_crs(Nod_abaq, A_ele, ele_L0, Fapp, moving, I, J, ele_type,0);
Nod_abaq(:,3)=Nod_abaq(:,3)-Nod_abaq(nodo_central,3);
the_error_spi=max(F_ele);

[Nod_abaq2, Ispis, Jspis, ele_L0spi, stp_spi]=FE_divide2(Nod_abaq, Ispis, Jspis, ele_L0(ele_type==1), fe_len);
[Nod_abaq2, Irads, Jrads, ele_L0rad, stp_rad]=FE_divide2(Nod_abaq2, Irads, Jrads, ele_L0(ele_type==2), fe_len);
[Nod_abaq2, Ifr, Jfr, ele_L0fr, stp_fr]=FE_divide2(Nod_abaq2, Ifr, Jfr, ele_L0(ele_type==3), fe_len);
[Nod_abaq2, Ianc, Janc, ele_L0anc, stp_anc]=FE_divide2(Nod_abaq2, Ianc, Janc, ele_L0(ele_type==4), fe_len);
% 
I=[Ianc,Ifr,Irads,Ispis];
J=[Janc,Jfr,Jrads,Jspis];
ele_L0=[ele_L0anc';ele_L0fr';ele_L0rad';ele_L0spi'];
ele_type=[4*ones(length(Ianc),1);3*ones(length(Ifr),1);2*ones(length(Irads),1);1*ones(length(Ispis),1)];
A_ele_new=0*ele_L0;
A_ele_new(ele_type==1)=A(1);
A_ele_new(ele_type==2)=A(2);
A_ele_new(ele_type==3)=A(3);
A_ele_new(ele_type==4)=A(4);

%% This is for the allocation of nod positions to stressed state of the web.
all_stps=[stp_spi;stp_rad;stp_fr;stp_anc];
Nod0=FE_divide_crs(Nod,I_crs,J_crs,ele_type0,all_stps);
save('spider_web_model.mat')