function check_TropOMI_quality()
load('C:\Projects\TropOMI\data\NO2_output\OFFL\Downsview_ERA\transport\TropOMI_transport');
DU = 2.6870e+16;
%% no2
% pixel_bins = NaN(height(TropOMI),450);
% for i=1:450
%     TF = TropOMI.ground_pixel == i;
%     if sum(TF) > 0
%         pixel_bins(TF,i) = TropOMI.no2(TF);
%     end
% end
% figure;boxplot(pixel_bins);
% 
% %% no2 err
% pixel_bins = NaN(height(TropOMI),450);
% for i=1:450
%     TF = TropOMI.ground_pixel == i;
%     if sum(TF) > 0
%         pixel_bins(TF,i) = TropOMI.no2_err(TF);
%     end
% end
% figure;boxplot(pixel_bins);
% 
% %% amf
% pixel_bins = NaN(height(TropOMI),450);
% for i=1:450
%     TF = TropOMI.ground_pixel == i;
%     if sum(TF) > 0
%         pixel_bins(TF,i) = TropOMI.amf(TF);
%     end
% end
% figure;boxplot(pixel_bins);
    
%% no2
pixel_bins = table;
bin_size = 15;
pixel_bins.no2 = NaN(height(TropOMI),450/bin_size);
pixel_bins.no2_err = NaN(height(TropOMI),450/bin_size);
pixel_bins.amf = NaN(height(TropOMI),450/bin_size);
for i=1:450/bin_size
    TF = (TropOMI.ground_pixel > (i-1)*bin_size) & (TropOMI.ground_pixel <= (i)*bin_size);
    if sum(TF) > 0
        pixel_bins.no2(TF,i) = TropOMI.no2(TF);
        pixel_bins.no2_err(TF,i) = TropOMI.no2_err(TF);
        pixel_bins.amf(TF,i) = TropOMI.amf(TF);
        pixel_bins.amf_trop(TF,i) = TropOMI.amf_trop(TF);
    end
end
figure;boxplot(pixel_bins.no2./DU);    
ylabel('TropOMI NO_2 [DU]');
xlabel('pixel bins');

figure;boxplot(pixel_bins.no2_err./DU);    
ylabel('TropOMI NO_2 err [DU]');
xlabel('pixel bins');

figure;boxplot(pixel_bins.amf);    
ylabel('TropOMI NO_2 AMF [DU]');
xlabel('pixel bins');

figure;boxplot(pixel_bins.amf_trop);    
ylabel('TropOMI NO_2 trop. AMF');
xlabel('pixel bins');
    