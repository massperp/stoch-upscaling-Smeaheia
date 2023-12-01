

% Script for generating Fig. 6
% Set kmod (4, 3, 1) to choose among the cases of clay permeability 

kmod = 4;

u_evals = combvec([0:0.01:1],[0:0.05:1])';

ind_matrix = nchoosek([1:5],2);
Np1 = 51;
Np2 = 51;
st1 = 1/(Np1-1);
st2 = 1/(Np2-1);
labels = {'$K$','$P_c$','$S$', '$K^{r}_{w}$', '$K^{r}_{nw}$'};

figure('Position',[100,100,600,800]);
t = tiledlayout(10,5,'TileSpacing','Compact','Padding','Compact');
s_d_pts = 21;
s_d_ind = 20; %16 %20;
zlims = [0 1; 0 1e8; 0 1; 0 1; 0 1];
if s_d_ind == 20
    zlims = [0 1; 0 1e5; 0.6 1; 0.4 1; 0 1];
elseif s_d_ind == 2 || s_d_ind == 5
    zlims = [0 1; 0 5e12; 0 0.1; 0 0.1; 0.8 1]; 
end


if kmod == 1
    kc = 1;
elseif kmod == 3
    kc = 0.001;
    zlims(1,:) = [0 2];
elseif kmod == 4
    kc = 0.0001;
end

    %figure('Position',[100,100,1500,300]);


for i=1:size(ind_matrix,1)
    pathname = ['../src/flow_fun_evals/variables_', num2str(ind_matrix(i,1)), num2str(ind_matrix(i,2)), '_k', num2str(kmod), '/'];
    
    for j=1:5
    
        f = fopen([pathname, 'copula_evals_Y_',num2str(j),'.txt']);
        temp = textscan(f,'%f');
        fclose(f);
        copula_samples{j} = temp{1};
    
        nexttile
        h = surf([0:st1:1], [0:st2:1],reshape(copula_samples{j}(s_d_ind:s_d_pts:end),Np2,Np1));

        axis([0 1 0 1 zlims(j,:)])
        set(h,'edgecolor','none')
        if i == 1
            title(labels{j},'FontWeight','Normal','interpreter','latex', 'FontSize',15)
        end
        view(-50,30) %view(-50,11)
        xlabel(['U_',num2str(ind_matrix(i,1))])
        ylabel(['U_',num2str(ind_matrix(i,2))])
        xticks([0 1])
        yticks([0 1])
        axes_handle = gca;
        
        if j==2
            hZ=axes_handle.ZRuler;
            hZ.SecondaryLabel.Position(1:3) = [0.5 1 0.9e5];
        end

        axes_handle.XLabel.Position(3) = zlims(j,1) - 0.05*(zlims(j,2)-zlims(j,1)); %axes_handle.XLabel.Position(3) + 0.4*zlims(j,2);
        axes_handle.YLabel.Position(3) = zlims(j,1) - 0.05*(zlims(j,2)-zlims(j,1)); %axes_handle.YLabel.Position(3) + 0.4*zlims(j,2);
    end 
end


figure('Position',[100,100,600,300]);
%figure('Position',[100,100,750,1500]);
t = tiledlayout(3,5,'TileSpacing','Compact','Padding','Compact');
ind_matrix_red = ind_matrix([1,2,5],:);
for i=1:size(ind_matrix_red,1)
    pathname = ['../src/flow_fun_evals/variables_', num2str(ind_matrix_red(i,1)), num2str(ind_matrix_red(i,2)), '_k', num2str(kmod), '/'];
    
    for j=1:5
    
        f = fopen([pathname, 'copula_evals_Y_',num2str(j),'.txt']);
        temp = textscan(f,'%f');
        fclose(f);
        copula_samples{j} = temp{1};
    
        nexttile
        h = surf([0:st1:1], [0:st2:1],reshape(copula_samples{j}(s_d_ind:s_d_pts:end),Np2,Np1));
        axis([0 1 0 1 zlims(j,:)]);
        set(h,'edgecolor','none')
        %title(['Y_',num2str(j)]); %, ' as fcn of U_', num2str(ind_matrix(i,1)+1),', U_', num2str(ind_matrix(i,2)+1)]);
        if i == 1
            title(labels{j},'FontWeight','Normal','interpreter','latex', 'FontSize',15)
        end
        view(-50,30) %view(-50,11)
        xlabel(['U_',num2str(ind_matrix_red(i,1))])
        ylabel(['U_',num2str(ind_matrix_red(i,2))])
        xticks([0 1])
        yticks([0 1])
        axes_handle = gca;
        
        if j==2
            hZ=axes_handle.ZRuler;
            hZ.SecondaryLabel.Position(1:3) = [0.5 1 0.9e5];
        end
        axes_handle.XLabel.Position(3) = zlims(j,1) - 0.05*(zlims(j,2)-zlims(j,1)); %axes_handle.XLabel.Position(3) + 0.4*zlims(j,2);
        axes_handle.YLabel.Position(3) = zlims(j,1) - 0.05*(zlims(j,2)-zlims(j,1)); %axes_handle.YLabel.Position(3) + 0.4*zlims(j,2);
    end 
end