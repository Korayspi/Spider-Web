function [Nod, Itype, Jtype, ele_L0n, all_stp]=FE_divide2(Nod, Itype, Jtype, ele_L0t, fe_len)

ele_no = length(Itype);
ele_len = vecnorm(Nod(Itype,:)-Nod(Jtype,:),2,2);

I0=Itype;
J0=Jtype;
Nod0=Nod;
ele_L0n=[];
all_stp=zeros(ele_no,1);
for i=1:ele_no
    stp=round(ele_len(i)/fe_len);
    all_stp(i)=stp;

    ele_n1=I0(i);
    ele_n2=J0(i);
    Itype(1)=[];
    Jtype(1)=[];
    nod_stp=(Nod0(ele_n1,:)-Nod0(ele_n2,:))/stp;
    
    N_no=length(Nod); 
    new_nods=(1:stp-1)'*nod_stp+ones(stp-1,1)*Nod0(ele_n2,:);
    Nod=[Nod; new_nods];
    new_ele=[ele_n2, N_no+(1:stp-1), ele_n1];
    for j=1:stp
        Itype(end+1)=new_ele(j);
        Jtype(end+1)=new_ele(j+1);
        ele_L0n(end+1)=ele_L0t(i)/stp;
    end
end
end