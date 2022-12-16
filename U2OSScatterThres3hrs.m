close all;
p53=[];
F1=[];
%Images are the wells we want to use to pick the threshold for p53. Using
%the log of the p53 signal works best. 
Images = [47:-6:5];
%Images = [[5:6:47]];
%Images = [[6:6:48]];
%Loop through  each image and store the p53 data in the p53 array
for i = 1:size(Images,2)
    p53 = [p53; log(handles.Measurements.Nuclei.Intensity_MeanIntensity_p53{Images(i)})];
    F1N = handles.Measurements.Nuclei.Intensity_MeanIntensity_Foxo1{Images(i)};
    F1C = handles.Measurements.Cytoplasm.Intensity_MeanIntensity_Foxo1{Images(i)};
    F1NCRatio = F1N./(F1N+F1C);
    F1= [F1; F1NCRatio];
end

update_to_branch = 5;
%Tried multithresh for the Foxo1 threshold but manually picking it at the
%trough between the two peaks worked better
%F1On = multithresh(F1);
F1On = .53;
p53On = -2.8;


%Display what the threshold is doing
figure();
[f,xi] = ksdensity(F1,[0:.01:1]);
[fp53,xip53] = ksdensity(p53,[-4:.01:-.5]);
subplot(2,2,2)
yline(F1On,'-.')
hold on
xline(p53On,'-.')
dscatter(p53,F1);
xlim([-4 -.5])
ylim([0 1])
title('All Data p53 v Foxo1')


subplot(2,2,1)
plot(xi,f)
hold on
xline(F1On,'-.')
camroll(-270)
title('ksDensity Foxo1')

subplot(2,2,4)
plot(xip53,fp53)
xlim([-4 -.5])
hold on
xline(p53On,'-.')
set(gca,'ydir','reverse')
title('ksDensity p53')


%Plot for 3 hours
Wells = [flipdim([6:6:48],2)];
ConcH2O2 = ['  0 uM';' 20 uM';' 40 uM';' 60 uM';' 80 uM';'100 uM';'150 uM';'200 uM'];
figure(100)
for i = 1:size(Wells,2)
    F1N = handles.Measurements.Nuclei.Intensity_MeanIntensity_Foxo1{Wells(i)};
    F1C = handles.Measurements.Cytoplasm.Intensity_MeanIntensity_Foxo1{Wells(i)};
    F1NCRatio = F1N./(F1N+F1C);
    if (size(F1NCRatio,1) >150)
        
        
        subplot(1,8,i);
        p53Temp = log(handles.Measurements.Nuclei.Intensity_MeanIntensity_p53{Wells(i)}); 
        dscatter(p53Temp,F1NCRatio);
        axis square
        hold on
        xline(p53On,'-.')
        yline(F1On,'-.')
        xlim([-3.8 -1.6])
        ylim([0 0.9])
        title(ConcH2O2(i,:))
        %title([temp{1} temp{2}])
        F1OnCells = (logical(F1NCRatio > F1On));
        p53OnCells = logical(p53Temp > p53On);
        F1Only(i) = sum(and(F1OnCells,~p53OnCells))/size(F1NCRatio,1);
        p53Only(i) = sum(and(~F1OnCells,p53OnCells))/size(F1NCRatio,1);
        Both(i) = sum(and(F1OnCells,p53OnCells))/size(F1NCRatio,1);
        Neither(i) = sum(and(~F1OnCells,~p53OnCells))/size(F1NCRatio,1);
        text(-3.7,.85,[num2str(round(F1Only(i)*100)) '%'],'FontSize',6);
        text(-2,.85,[num2str(round(Both(i)*100)) '%'],'FontSize',6);
        text(-2,.075,[num2str(round(p53Only(i)*100)) '%'],'FontSize',6);
        text(-3.7,.075,[num2str(round(Neither(i)*100)) '%'],'FontSize',6);
    end
end