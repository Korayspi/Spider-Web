function [Nod, I, J, ele_type, nodo_central, ang_rad, envi_nods]=MakeWeb_full_nospi(n_anc,n_rad,L_anc,L_fra)

%% Anchor Threads
    thet=2*pi/n_anc; %
    phi=(n_anc-2)*pi/(2*n_anc);
    thet_rad = 2*pi/n_rad;
    phi_rad = thet_rad/2;
    I=[];
    J=[];
    ang_anc=zeros(1,n_anc);
    for i=1:n_anc
        ang_anc(i)=phi+i*thet;
        Nod(i,:) = L_anc * [cos(phi+i*thet),sin(phi+i*thet),0];
        Nod(i+n_anc,:) = L_fra * [cos(phi+i*thet),sin(phi+i*thet),0];
        I=[I, i];
        J=[J, i+n_anc];
    end

% Deciding radial attachment nods on the frame
    alpha=(n_anc-2)*pi/(2*n_anc);
    ang_rad=zeros(1,n_rad);
    for i=1:n_rad
%        Nod(2*n_anc+i,:)=(L_fra*sin(alpha) / sin(pi-alpha-mod(-phi+phi_rad+i*thet_rad,2*pi/n_anc))) * [cos(phi_rad+i*thet_rad),sin(phi_rad+i*thet_rad),0];
        Nod(2*n_anc+i,:)=(L_fra*sin(alpha) / sin(pi-alpha-mod(-phi+i*thet_rad,2*pi/n_anc))) * [cos(i*thet_rad),sin(i*thet_rad),0];
        ang_rad(i)=i*thet_rad;
    end

    for i=1:n_anc
        d_ang=[abs(ang_anc(i)-ang_rad-2*pi); abs(ang_anc(i)-ang_rad)];
        d_ang=min(d_ang,[],1);
        [M_ang, I_ang]=mink(d_ang,2);
        I=[I, n_anc+i, n_anc+i];
        J=[J, 2*n_anc+I_ang(1),2*n_anc+I_ang(2)];
    end
    % Center Nod
    envi_nods=Nod;
    nodo_central=length(Nod)+1;
    Nod=[Nod; 0, 0, 0];

    rem=zeros(2,n_anc);
    rem(1:end)=[J(end:-1:end-2*n_anc+1)];
    rem=[rem(1,:),rem(2,:);rem(2,:),rem(1,:)];
    ind=length(I);
    I=[I,zeros(1,2*n_rad)];
    J=[J,zeros(1,2*n_rad)];
    t=ind-1;

    for i=1:n_rad
        t=t+2;
        I(t)=2*n_anc+i;
        I(t+1)=2*n_anc+i;
        J(t)=2*n_anc+1+mod(i,n_rad);
        J(t+1)=length(Nod);
    end

    for i=1:length(rem)
        T=rem(:,i)==[I;J];
        r=find(sum(T)==2);
        I(r)=[];
        J(r)=[];
    end
    ele_type=zeros(1,length(I));
    ele_type(1:n_anc)=4;
    ele_type(J==length(Nod))=2;
    ele_type(ele_type==0)=3;
% PlotWeb3D(Nod,I,J)
% for i=1:length(Nod)
%     text(Nod(i,1),Nod(i,2),Nod(i,3),num2str(i),'FontSize', 12)
% end
% 
% ele_pos=(Nod(I,:)+Nod(J,:))/2;
% for i=1:length(I)
%     text(ele_pos(i,1),ele_pos(i,2),ele_pos(i,3),num2str(ele_type(i)),'FontSize', 12)
% end
end