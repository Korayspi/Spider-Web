function [Nod, I, J, ele_type, spi_Nod_no]=MakeWeb_full3(Nod,I,J,ele_type,n_rad,spider_size,rmv_spi,nodo_central)

    Ianc=I(ele_type==4);
    Janc=J(ele_type==4);
    Ifr=I(ele_type==3);
    Jfr=J(ele_type==3);
    Irad0=I(ele_type==2);
    Jrad0=J(ele_type==2);

    for i=1:n_rad
        rad(i).nod=[];
    end
    
    % Here I'm defining capture area
    len_cptr=vecnorm(Nod(Irad0,:),2,2);
    cptr_rad_diff=min(len_cptr)-len_cptr; % This ensures that the capture area is circular

    pos_rad=Nod(Irad0,:)-Nod(nodo_central,:);
    rad_len=vecnorm(pos_rad,2,2);

    vec_rad=pos_rad./rad_len; % This defines the position of the nodes on the radial
    rad_len=rad_len+cptr_rad_diff; % This ensures that the capture area is circular

    spi_Nod_no=floor((rad_len+0.00005-0.5)/spider_size); %How many spiral node?
    t=1;
    radNod=[];
    for i=1:n_rad
        for j=(1+rmv_spi):spi_Nod_no(i)
            radNod(end+1,:)=Nod(nodo_central,:)+j*spider_size*vec_rad(i,:);
            rad(i).nod(end+1)=t;
            t=t+1;
        end
        radNod_edge(i)=t;
%         Nod(end+1,:)=radNod(i,:);
%         rad(i).nod(end+1)=t;
%         t=t+1;
    end
    Nod_no=length(Nod);
    Nod=[Nod;radNod];
    A=max(spi_Nod_no);
%% Radials
Irad=[];
Jrad=[];
for i=1:n_rad
    % All radials connect to frame

    rad(i).nod=rad(i).nod+Nod_no;
    % All radials connect to hub
    Irad(end+1)=nodo_central;
    Jrad(end+1)=rad(i).nod(1);

    for j=1:length(rad(i).nod)-1
        Irad(end+1)=rad(i).nod(j);
        Jrad(end+1)=rad(i).nod(j+1);
    end
    Irad(end+1)=rad(i).nod(j+1);
    Jrad(end+1)=Irad0(i);
end

%% Spirals
Ispi=[];
Jspi=[];

for j=1:A-rmv_spi
    ind=find(spi_Nod_no>=A);
    t=1;
    while t<length(ind)
        k=0;
        if abs(ind(t)-ind(t+1))==1
            Ispi(end+1)=rad(ind(t)).nod(A-rmv_spi);
            Jspi(end+1)=rad(ind(t+1)).nod(A-rmv_spi);
        else
            k=k+1;
        end
        t=t+1;
        if k~=0
            if abs(ind(end)-ind(1))==n_rad-1
                Ispi(end+1)=rad(ind(end)).nod(A-rmv_spi);
                Jspi(end+1)=rad(ind(1)).nod(A-rmv_spi);
            end
        end
    end
    if length(ind)==n_rad
        Ispi(end+1)=rad(ind(end)).nod(A-rmv_spi);
        Jspi(end+1)=rad(ind(1)).nod(A-rmv_spi);
    end
    A=A-1;
end
I=[Ianc,Ifr,Irad,Ispi];
J=[Janc,Jfr,Jrad,Jspi];
type_inds=[length(Ianc),length([Ianc,Ifr]),length([Ianc,Ifr,Irad])];
ele_type=zeros(1,length(I));
ele_type(1:type_inds(1))=4;
ele_type(type_inds(1)+1:type_inds(2))=3;
ele_type(type_inds(2)+1:type_inds(3))=2;
ele_type(ele_type==0)=1;
% PlotWeb3D(Nod,Ispi,Jspi)
% %  PlotWeb3D(Nod,I,J)
% % % for i=1:length(Nod)
% % %     text(Nod(i,1),Nod(i,2),Nod(i,3),num2str(i),'FontSize', 12)
% % % end
% % 
% ele_pos=(Nod(I,:)+Nod(J,:))/2;
% for i=1:length(I)
%     text(ele_pos(i,1),ele_pos(i,2),ele_pos(i,3),num2str(ele_type(i)),'FontSize', 12)
% end
end