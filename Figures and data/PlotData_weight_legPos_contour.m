clc; clear;
weight_data=load('results_m_tens.mat');
% weight_data=load('results_w_tens.mat');

results_weight=weight_data.results;


sgn_rads=[4 5 6 7 1 30 29 26 25 24 23];
iassym=0;
weights=[0.001,250,500,750,1000]*1e-6;
tens_val=[100,125,150,175,200];

tens_name={'100','125','150','175','200'};
weight_names={'0','25','50','75','100'};
weight_names_axis={'0mg', '25mg', '50mg', '75mg', '100mg'};

diff_vals=zeros(5);
for t=1
    for tens=1:5
        for w=1:5
%             t+1+3*(w-1)+15*(tens-1)
            pos_data=results_weight(t+1+3*(tens-1)+15*(w-1)).leg_pos(sgn_rads(9+t),:);
            name_data=results_weight(t+1+3*(tens-1)+15*(w-1)).name
            pos_data=pos_data-pos_data(1);
            [fr1, amp1]=returnFFT((1:250)*(0.1/250),tukeywin(250,0.2)'.*pos_data,2500);

            pos_data1=results_weight(t+1+3*(tens-1)+15*(w-1)).leg_pos(sgn_rads(9+t)-1,:);
            pos_data1=pos_data1-pos_data1(1);
            [fr2, amp2]=returnFFT((1:250)*(0.1/250),tukeywin(250,0.2)'.*pos_data1,2500);

            max_val=0.5*(amp1(21)+amp2(21));
            diff_vals(w,tens)=max_val;        
        end
    end
    
    FigH = figure('Position', get(0, 'Screensize'));
    contour(weights*1e5,tens_val,diff_vals',[0.4e-5, 0.75e-5, 0.85e-5],'LineWidth',3,'Color','k')
    ylabel('Tension on Radials (\mu N)');
    xlabel('Spider Weight')
    title(['Radial: ',num2str(sgn_rads(9+t))]);
%     set(gca,'xticklabel',weight_names_axis)
    ax = gca; 
    ax.FontSize = 22;
    axis square
%     saveas(FigH, ['Radial_',num2str(sgn_rads(9+t)),'legpos.eps'],'eps');
end