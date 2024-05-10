function [Fnod, F_elem, Nodm, max_val]=relax_6err(Nodm, A, ele_L0, Fapp, moving, I, J, ele_type, error)
    
    dx_un=0.1;
    vars(1)=1.01;
    vars(2)=1e-5;
    vars(3)=10;
    div_un=vars(1);
    
    dx=dx_un*ones(length(moving),3);    
    sign_val=0*dx;
    
    t=1;
    
    [Fnod, F_elem] = ForceCalc_Z(Fapp, Nodm, A, ele_L0, I, J, ele_type);
    Nod_best=Nodm;
    Fnod_best=Fnod;
    F_elem_best=F_elem;
    Fnod(abs(Fnod)<error)=0;

count=100000;
max_val=max(vecnorm(Fnod(moving,:),2,2));
    while t<=count
        if max_val<1e-7
            break
        end
        t=t+1;
        Fnodp=Fnod;
        Nodm(moving,:) = Nodm(moving,:) + sign(Fnod(moving,:)).*dx;
        [Fnod, F_elem] = ForceCalc_Z(Fapp, Nodm, A, ele_L0, I, J, ele_type);
        Fnod(abs(Fnod)<error)=0;
        dx(Fnod(moving,:).*Fnodp(moving,:)<0)=dx(Fnod(moving,:).*Fnodp(moving,:)<0)/div_un;
        sign_val=sign_val+sign(Fnod(moving,:).*dx);
        if mod(t,vars(3))==0
            chk_val=abs(sign_val);
            dx(chk_val==vars(3))=dx(chk_val==vars(3))+vars(2);
            sign_val=0*dx;
        end
        max_valp=max_val;
        max_val=max(vecnorm(Fnod(moving,:),2,2));
        if max_val<max_valp
            Nod_best=Nodm;
            Fnod_best=Fnod;
            F_elem_best=F_elem;
        end
        if max_val<1e-4
            vars(2)=5e-7;
        end
        if max_val<1e-5
            vars(2)=1e-7;
        end
    end
Nodm=Nod_best;
Fnod=Fnod_best;
F_elem=F_elem_best;
end