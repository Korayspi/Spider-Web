clc; clear;
weight_data=load('results_m.mat');

web=load('assym0_0p2spi_100T.mat');
Nod=web.Nod;
I=web.I_crs;
J=web.J_crs;
A=web.A;
center_cons=find(I==web.nodo_central);

rad_vec=Nod(J(center_cons),[1,2])-Nod(web.nodo_central,[1,2]);

results_weight=weight_data.results;


sgn_rads=[4 5 6 7 1 30 29 26 25 24 23];
iassym=0;
weights=[0.001,250,500,750,1000]*1e-6;
tens_val=[100,125,150,175,200];

tens_name={'100','125','150','175','200'};
weight_names={'0','25','50','75','100'};
weight_names_axis={'0mg', '25mg', '50mg', '75mg', '100mg'};


diff_vals=zeros(1,30);
diff_vals_str=zeros(1,30);

for s=1:11
    for w=1:5
        name_data=results_weight(s+11*(w-1)).name
        for i=1:30
            pos_data=results_weight(s+11*(w-1)).leg_pos(i,:);                
            pos_data=pos_data-pos_data(1);

            str_data=results_weight(s+11*(w-1)).leg_str(i,:);                
            str_data=str_data-str_data(1);

            [fr1, amp1]=returnFFT((1:250)*(0.1/250),tukeywin(250,0.2)'.*pos_data(1:250),2500);
            [fr1_str, amp1_str]=returnFFT((1:250)*(0.1/250),tukeywin(250,0.2)'.*str_data(1:250),2500);

            max_val=amp1(21);
            diff_vals(i)=max_val;
            diff_vals_str(i)=amp1_str(41);
        end
        FigH = figure('Position', get(0, 'Screensize'));
        PlotWeb2D(Nod,I,J)
        axis off
        hold on
        amp=rad_vec.*diff_vals'*4e5;
        amp=[amp;amp(1,:)];
        plot(amp(:,1),amp(:,2),'LineWidth',3,'Color','b')
        hold off        
        ax = gca; 
        ax.FontSize = 22;
%         saveas(FigH, ['2Dweb_',num2str(sgn_rads(s)),'sgn_',weight_names{w},'m.eps'],'eps');
    end
end