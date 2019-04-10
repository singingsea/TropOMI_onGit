function residual_check()
load('C:\Projects\TropOMI\plots\4Vitali\scatter.mat');
DU = 2.6870e+16;

data = data_Downsview_NRTI;
residual = data.TropOMI_NO2./DU - 0.676.*data.Pandora_NO2;
err = std(residual);
disp(['Downsview NRTI err = ' num2str(err) ' DU']);
disp(['Downsview NRTI err = ' num2str(err*DU,'%e') ' molec/cm^2']);

data = data_Downsview_OFFL;
residual = data.TropOMI_NO2./DU - 0.883.*data.Pandora_NO2;
err = std(residual);
disp(['Downsview OFFL err = ' num2str(err) ' DU']);
disp(['Downsview OFFL err = ' num2str(err*DU,'%e') ' molec/cm^2']);

data = data_FortMcKay_NRTI;
residual = data.TropOMI_NO2./DU - 0.548.*data.Pandora_NO2;
err = std(residual);
disp(['FortMcKay NRTI err = ' num2str(err) ' DU']);
disp(['FortMcKay NRTI err = ' num2str(err*DU,'%e') ' molec/cm^2']);

data = data_FortMcKay_OFFL;
residual = data.TropOMI_NO2./DU - 0.750.*data.Pandora_NO2;
err = std(residual);
disp(['FortMcKay OFFL err = ' num2str(err) ' DU']);
disp(['FortMcKay OFFL err = ' num2str(err*DU,'%e') ' molec/cm^2']);
