function [Fnod, F_ele] = ForceCalc_Z(Fapp, Nod, A, ele_L0, I, J, ele_type)
    Fnod=Fapp;
    ele_vec = Nod(J,:)-Nod(I,:); %element between two nodes
    len = vecnorm(ele_vec,2,2);
    vec_xy = ele_vec./len; %x-y vector of the elements

    strain=log(len./ele_L0);
    push_force=strain<0;
%     strain(strain<0)=0;

    spis=(ele_type==1);
    non_spis=(ele_type~=1);
    
    F_ele=0*strain;
    F_ele(spis) = A(spis).*get_stress_real(strain(spis),1);
    F_ele(non_spis) = A(non_spis).*get_stress_real(strain(non_spis),2);
    F_ele(push_force)=real(F_ele(push_force));
    for i=1:length(I)
        Fnod(I(i),:) = Fnod(I(i),:) + F_ele(i,:).*vec_xy(i,:);
        Fnod(J(i),:) = Fnod(J(i),:) - F_ele(i,:).*vec_xy(i,:);
    end
end