function Nod=FE_divide_crs(Nod,Icrs,Jcrs,ele_type,all_stp)
t=0;
for index=1:4    
    I0=Icrs(ele_type==index);
    J0=Jcrs(ele_type==index);
    ele_no = length(I0);
    Nod0=Nod;
    for i=1:ele_no
        t=t+1;
        stp=all_stp(t);
        ele_n1=I0(i);
        ele_n2=J0(i);
        nod_stp=(Nod0(ele_n1,:)-Nod0(ele_n2,:))/stp;
        new_nods=(1:stp-1)'*nod_stp+ones(stp-1,1)*Nod0(ele_n2,:);
        Nod=[Nod; new_nods];
    end
end
end