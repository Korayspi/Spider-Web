clc; clear;
weight_data=load('results_m.mat');
% weight_data=load('results_w.mat');

results_weight=weight_data.results;

sgn_rads=[4 5 6 7 1 30 29 26 25 24 23];
iassym=0;
weights=[0.001,250,500,750,1000]*1e-6;
tens_val=[100,125,150,175,200];

tens_name={'100','125','150','175','200'};
weight_names={'0','25','50','75','100'};
weight_names_axis={'0mg', '25mg', '50mg', '75mg', '100mg'};
wind_names_axis={'0m/s', '0.5m/s', '1m/s','2m/s','4m/s'};

diff_vals=zeros(1,5);
sgn2noise=zeros(1,5);
for s=10
    for w=1:5
        name_data=results_weight(s+11*(w-1)).name
        pos_data=results_weight(s+11*(w-1)).leg_pos(sgn_rads(s),:);
        pos_data=pos_data-pos_data(1);

        [fr1, amp1]=returnFFT((1:250)*(0.1/250),tukeywin(250,0.2)'.*pos_data,2500);
        
        if sgn_rads(s)==1
            pos_data1=results_weight(s+11*(w-1)).leg_pos(30,:);
        else
            pos_data1=results_weight(s+11*(w-1)).leg_pos(sgn_rads(s)-1,:);
        end
        mean(pos_data1)
        pos_data1=pos_data1-pos_data1(1);

        [fr2, amp2]=returnFFT((1:250)*(0.1/250),tukeywin(250,0.2)'.*pos_data1,2500);

        max_val=0.5*amp1(21)+0.5*amp2(21);
        diff_vals(w)=max_val;
    end
    FigH = figure('Position', get(0, 'Screensize'));
    bar(1:5,diff_vals,'LineWidth',3)
    set(gca,'xticklabel',weight_names_axis)
    ax = gca; 
    ax.FontSize = 22;
    
    ylabel('Spide (\mu N)');
    xlabel('Spider Weight')
    set(gca,'xticklabel',weight_names_axis)
%     saveas(FigH, ['Radial_',num2str(sgn_rads(s)),'_bargraph.eps'],'eps');
end
