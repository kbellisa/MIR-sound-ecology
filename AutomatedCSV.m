
%[metricskey_filename, metricskey_pathname] = uigetfile('*');
metricskey_fullpath = '\\datadepot.rcac.purdue.edu\depot\bpijanow\data\CGS_Projects\MetricsSoundBankData\MetricsKey_edited.xlsx' %[metricskey_pathname metricskey_filename];

key_time_range = 'B2:B60';
key_Time_range = 'C2:C60';
key_Day_range = 'D2:D60';
key_Site_range = 'E2:E60';


times = xlsread(metricskey_fullpath, key_time_range);
[numT, strT, Times] = xlsread(metricskey_fullpath, key_Time_range);
[numD, strD, Days] = xlsread(metricskey_fullpath, key_Day_range);
[numS, Sites, rawS] = xlsread(metricskey_fullpath, key_Site_range);


for nn = 1:length(Sites)
   pathnames(nn, :) = strcat(Sites(nn), '_', Days(nn), '_' , Times(nn));
end

arizona_pathnames = pathnames(1:27);
arizona_times = times(1:27);
laselva_pathnames = pathnames(28:44);
laselva_times = times(28:44);
wellspilger_pathnames = pathnames(45:59);
wellspilger_times = times(45:59);
    

'PATHNAMES LOADED'
%ARIZONA (SC, SK, SL, SP, VAR)

%SC

az_foldername = '\\datadepot.rcac.purdue.edu\depot\bpijanow\data\CGS_Projects\MetricsSoundBankData\METRICSsoundBank\AZ_FE\';
cr_foldername = '\\datadepot.rcac.purdue.edu\depot\bpijanow\data\CGS_Projects\MetricsSoundBankData\METRICSsoundBank\CR_FE\';
me_foldername = '\\datadepot.rcac.purdue.edu\depot\bpijanow\data\CGS_Projects\MetricsSoundBankData\METRICSsoundBank\ME_FE\';

for nn = 1:length(arizona_pathnames)
   arizona_sc_pathnames(nn, :) = strcat(az_foldername,'AZ_SC/', arizona_pathnames(nn), '_vamp_vamp-libxtract_spectral_centroid_spectral_centroid.csv');
   arizona_sk_pathnames(nn, :) = strcat(az_foldername,'AZ_Sk/', arizona_pathnames(nn), '_vamp_vamp-libxtract_spectral_skewness_spectral_skewness.csv');
   arizona_sl_pathnames(nn, :) = strcat(az_foldername,'AZ_Sl/', arizona_pathnames(nn), '_vamp_vamp-libxtract_spectral_slope_spectral_slope.csv');
   arizona_sp_pathnames(nn, :) = strcat(az_foldername,'AZ_Sp/', arizona_pathnames(nn), '_vamp_vamp-libxtract_spread_spread.csv');
   arizona_var_pathnames(nn, :) = strcat(az_foldername,'AZ_Var/', arizona_pathnames(nn), '_vamp_vamp-libxtract_spectral_variance_spectral_variance.csv');
end

for nn = 1:length(laselva_pathnames)
   laselva_sc_pathnames(nn, :) = strcat(cr_foldername,'CR_SC/', laselva_pathnames(nn), '_vamp_vamp-libxtract_spectral_centroid_spectral_centroid.csv');
   laselva_sk_pathnames(nn, :) = strcat(cr_foldername,'CR_Sk/', laselva_pathnames(nn), '_vamp_vamp-libxtract_spectral_skewness_spectral_skewness.csv');
   laselva_sl_pathnames(nn, :) = strcat(cr_foldername,'CR_Sl/', laselva_pathnames(nn), '_vamp_vamp-libxtract_spectral_slope_spectral_slope.csv');
   laselva_sp_pathnames(nn, :) = strcat(cr_foldername,'CR_Sp/', laselva_pathnames(nn), '_vamp_vamp-libxtract_spread_spread.csv');
   laselva_var_pathnames(nn, :) = strcat(cr_foldername,'CR_Var/', laselva_pathnames(nn), '_vamp_vamp-libxtract_spectral_variance_spectral_variance.csv');
end


for nn = 1:length(wellspilger_pathnames)
   wells_sc_pathnames(nn, :) = strcat(me_foldername,'ME_SC/', wellspilger_pathnames(nn), '_vamp_vamp-libxtract_spectral_centroid_spectral_centroid.csv');
   wells_sk_pathnames(nn, :) = strcat(me_foldername,'ME_Sk/', wellspilger_pathnames(nn), '_vamp_vamp-libxtract_spectral_skewness_spectral_skewness.csv');
   wells_sl_pathnames(nn, :) = strcat(me_foldername,'ME_Sl/', wellspilger_pathnames(nn), '_vamp_vamp-libxtract_spectral_slope_spectral_slope.csv');
   wells_sp_pathnames(nn, :) = strcat(me_foldername,'ME_Sp/', wellspilger_pathnames(nn), '_vamp_vamp-libxtract_spread_spread.csv');
   wells_var_pathnames(nn, :) = strcat(me_foldername,'ME_Var/', wellspilger_pathnames(nn), '_vamp_vamp-libxtract_spectral_variance_spectral_variance.csv');
end

'BEGINNING COMPUTATION'

   az_comp = [];
   cr_comp = [];
   me_comp = [];
   
for nn = 1:length(arizona_pathnames)
    %ARIZONA DATA
    az_sc_data = csvread(arizona_sc_pathnames{nn, :});
    az_sk_data = csvread(arizona_sk_pathnames{nn, :});
    az_sl_data = csvread(arizona_sl_pathnames{nn, :});
    az_sp_data = csvread(arizona_sp_pathnames{nn, :});
    az_var_data = csvread(arizona_var_pathnames{nn, :});
    
    %ARIZONA VALUES
    az_sc_times = az_sc_data(:, 1);
    az_sc_vals = az_sc_data(:, 2);
    
    az_sk_times = az_sk_data(:, 1);
    az_sk_vals = az_sk_data(:, 2);
    
    az_sl_times = az_sl_data(:, 1);
    az_sl_vals = az_sl_data(:, 2);
    
    az_sp_times = az_sp_data(:, 1);
    az_sp_vals = az_sp_data(:, 2);
    
    az_var_times = az_var_data(:, 1);
    az_var_vals = az_var_data(:, 2);
    
    sc_index = find(abs(az_sc_times - times(nn)) < 0.01);
        if sc_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            az_comp(nn, 1) = mean(az_sc_vals(lower_limit:upper_limit));

        elseif sc_index > length(az_sc_data) - 65
            lower_limit = length(az_sc_data) - 131;
            upper_limit = length(az_sc_data);
            az_comp(nn, 1) = mean(az_sc_vals(lower_limit:upper_limit));
        else
            lower_limit = sc_index - 65;
            upper_limit = sc_index + 65;
            az_comp(nn, 1) = mean(az_sc_vals(lower_limit:upper_limit));
        end
        
         sk_index = find(abs(az_sk_times - times(nn)) < 0.01);
        if sk_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            az_comp(nn, 2) = mean(az_sk_vals(lower_limit:upper_limit));

        elseif sk_index > length(az_sk_data) - 65
            lower_limit = length(az_sk_data) - 131;
            upper_limit = length(az_sk_data);
            az_comp(nn, 2) = mean(az_sk_vals(lower_limit:upper_limit));
        else
            lower_limit = sk_index - 65;
            upper_limit = sk_index + 65;
            az_comp(nn, 2) = mean(az_sk_vals(lower_limit:upper_limit));
        end
        
         sl_index = find(abs(az_sl_times - times(nn)) < 0.01);
        if sl_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            az_comp(nn, 3) = mean(az_sl_vals(lower_limit:upper_limit));

        elseif sl_index > length(az_sl_data) - 65
            lower_limit = length(az_sl_data) - 131;
            upper_limit = length(az_sl_data);
            az_comp(nn, 3) = mean(az_sl_vals(lower_limit:upper_limit));
        else
            lower_limit = sl_index - 65;
            upper_limit = sl_index + 65;
            az_comp(nn, 3) = mean(az_sl_vals(lower_limit:upper_limit));
        end
        
         sp_index = find(abs(az_sp_times - times(nn)) < 0.01);
        if sp_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            az_comp(nn, 4) = mean(az_sp_vals(lower_limit:upper_limit));

        elseif sc_index > length(az_sp_data) - 65
            lower_limit = length(az_sp_data) - 131;
            upper_limit = length(az_sp_data);
            az_comp(nn, 4) = mean(az_sp_vals(lower_limit:upper_limit));
        else
            lower_limit = sp_index - 65;
            upper_limit = sp_index + 65;
            az_comp(nn, 4) = mean(az_sp_vals(lower_limit:upper_limit));
        end
        
         var_index = find(abs(az_var_times - times(nn)) < 0.01);
        if var_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            az_comp(nn, 5) = mean(az_var_vals(lower_limit:upper_limit));

        elseif var_index > length(az_var_data) - 65
            lower_limit = length(az_var_data) - 131;
            upper_limit = length(az_var_data);
            az_comp(nn, 5) = mean(az_var_vals(lower_limit:upper_limit));
        else
            lower_limit = var_index - 65;
            upper_limit = var_index + 65;
            az_comp(nn, 5) = mean(az_var_vals(lower_limit:upper_limit));
        end   
end


for nn = 1:length(laselva_pathnames)
 %LASELVA DATA
    cr_sc_data = csvread(laselva_sc_pathnames{nn, :});
    cr_sk_data = csvread(laselva_sk_pathnames{nn, :});
    cr_sl_data = csvread(laselva_sl_pathnames{nn, :});
    cr_sp_data = csvread(laselva_sp_pathnames{nn, :});
    cr_var_data = csvread(laselva_var_pathnames{nn, :});
    
    %LASELVA VALUES
    cr_sc_times = cr_sc_data(:, 1);
    cr_sc_vals = cr_sc_data(:, 2);
    
    cr_sk_times = cr_sk_data(:, 1);
    cr_sk_vals = cr_sk_data(:, 2);
    
    cr_sl_times = cr_sl_data(:, 1);
    cr_sl_vals = cr_sl_data(:, 2);
    
    cr_sp_times = cr_sp_data(:, 1);
    cr_sp_vals = cr_sp_data(:, 2);
    
    cr_var_times = cr_var_data(:, 1);
    cr_var_vals = cr_var_data(:, 2);
    
        sc_index = find(abs(cr_sc_times - times(nn)) < 0.01);
        if sc_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            cr_comp(nn, 1) = mean(cr_sc_vals(lower_limit:upper_limit));

        elseif sc_index > length(cr_sc_data) - 65
            lower_limit = length(cr_sc_data) - 131;
            upper_limit = length(cr_sc_data);
            cr_comp(nn, 1) = mean(cr_sc_vals(lower_limit:upper_limit));
        else
            lower_limit = sc_index - 65;
            upper_limit = sc_index + 65;
            cr_comp(nn, 1) = mean(cr_sc_vals(lower_limit:upper_limit));
        end
        
         sk_index = find(abs(cr_sk_times - times(nn)) < 0.01);
        if sk_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            cr_comp(nn, 2) = mean(cr_sk_vals(lower_limit:upper_limit));

        elseif sk_index > length(cr_sk_data) - 65
            lower_limit = length(cr_sk_data) - 131;
            upper_limit = length(cr_sk_data);
            cr_comp(nn, 2) = mean(cr_sk_vals(lower_limit:upper_limit));
        else
            lower_limit = sk_index - 65;
            upper_limit = sk_index + 65;
            cr_comp(nn, 2) = mean(cr_sk_vals(lower_limit:upper_limit));
        end
        
         sl_index = find(abs(cr_sl_times - times(nn)) < 0.01);
        if sl_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            cr_comp(nn, 3) = mean(cr_sl_vals(lower_limit:upper_limit));

        elseif sl_index > length(cr_sl_data) - 65
            lower_limit = length(cr_sl_data) - 131;
            upper_limit = length(cr_sl_data);
            cr_comp(nn, 3) = mean(cr_sl_vals(lower_limit:upper_limit));
        else
            lower_limit = sl_index - 65;
            upper_limit = sl_index + 65;
            cr_comp(nn, 3) = mean(cr_sl_vals(lower_limit:upper_limit));
        end
        
         sp_index = find(abs(cr_sp_times - times(nn)) < 0.01);
        if sp_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            cr_comp(nn, 4) = mean(cr_sp_vals(lower_limit:upper_limit));

        elseif sc_index > length(cr_sp_data) - 65
            lower_limit = length(cr_sp_data) - 131;
            upper_limit = length(cr_sp_data);
            cr_comp(nn, 4) = mean(cr_sp_vals(lower_limit:upper_limit));
        else
            lower_limit = sp_index - 65;
            upper_limit = sp_index + 65;
            cr_comp(nn, 4) = mean(cr_sp_vals(lower_limit:upper_limit));
        end
        
         var_index = find(abs(cr_var_times - times(nn)) < 0.01);
        if var_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            cr_comp(nn, 5) = mean(cr_var_vals(lower_limit:upper_limit));

        elseif var_index > length(cr_var_data) - 65
            lower_limit = length(cr_var_data) - 131;
            upper_limit = length(cr_var_data);
            cr_comp(nn, 5) = mean(cr_var_vals(lower_limit:upper_limit));
        else
            lower_limit = var_index - 65;
            upper_limit = var_index + 65;
            cr_comp(nn, 5) = mean(cr_var_vals(lower_limit:upper_limit));
        end
        
    
end

for nn = 1:length(wellspilger_pathnames)
 
    %WELLS DATA
    me_sc_data = csvread(wells_sc_pathnames{nn, :});
    me_sk_data = csvread(wells_sk_pathnames{nn, :});
    me_sl_data = csvread(wells_sl_pathnames{nn, :});
    me_sp_data = csvread(wells_sp_pathnames{nn, :});
    me_var_data = csvread(wells_var_pathnames{nn, :});
    
    %WELLS VALUES
    me_sc_times = me_sc_data(:, 1);
    me_sc_vals = me_sc_data(:, 2);
    
    me_sk_times = me_sk_data(:, 1);
    me_sk_vals = me_sk_data(:, 2);
    
    me_sl_times = me_sl_data(:, 1);
    me_sl_vals = me_sl_data(:, 2);
    
    me_sp_times = me_sp_data(:, 1);
    me_sp_vals = me_sp_data(:, 2);
    
    me_var_times = me_var_data(:, 1);
    me_var_vals = me_var_data(:, 2);
    
        sc_index = find(abs(me_sc_times - times(nn)) < 0.01);
        if sc_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            me_comp(nn, 1) = mean(me_sc_vals(lower_limit:upper_limit));

        elseif sc_index > length(me_sc_data) - 65
            lower_limit = length(me_sc_data) - 131;
            upper_limit = length(me_sc_data);
            me_comp(nn, 1) = mean(me_sc_vals(lower_limit:upper_limit));
        else
            lower_limit = sc_index - 65;
            upper_limit = sc_index + 65;
            me_comp(nn, 1) = mean(me_sc_vals(lower_limit:upper_limit));
        end
        
         sk_index = find(abs(me_sk_times - times(nn)) < 0.01);
        if sk_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            me_comp(nn, 2) = mean(me_sk_vals(lower_limit:upper_limit));

        elseif sk_index > length(me_sk_data) - 65
            lower_limit = length(me_sk_data) - 131;
            upper_limit = length(me_sk_data);
            me_comp(nn, 2) = mean(me_sk_vals(lower_limit:upper_limit));
        else
            lower_limit = sk_index - 65;
            upper_limit = sk_index + 65;
            me_comp(nn, 2) = mean(me_sk_vals(lower_limit:upper_limit));
        end
        
         sl_index = find(abs(me_sl_times - times(nn)) < 0.01);
        if sl_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            me_comp(nn, 3) = mean(me_sl_vals(lower_limit:upper_limit));

        elseif sl_index > length(me_sl_data) - 65
            lower_limit = length(me_sl_data) - 131;
            upper_limit = length(me_sl_data);
            me_comp(nn, 3) = mean(me_sl_vals(lower_limit:upper_limit));
        else
            lower_limit = sl_index - 65;
            upper_limit = sl_index + 65;
            me_comp(nn, 3) = mean(me_sl_vals(lower_limit:upper_limit));
        end
        
         sp_index = find(abs(me_sp_times - times(nn)) < 0.01);
        if sp_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            me_comp(nn, 4) = mean(me_sp_vals(lower_limit:upper_limit));

        elseif sc_index > length(me_sp_data) - 65
            lower_limit = length(me_sp_data) - 131;
            upper_limit = length(me_sp_data);
            me_comp(nn, 4) = mean(me_sp_vals(lower_limit:upper_limit));
        else
            lower_limit = sp_index - 65;
            upper_limit = sp_index + 65;
            me_comp(nn, 4) = mean(me_sp_vals(lower_limit:upper_limit));
        end
        
         var_index = find(abs(me_var_times - times(nn)) < 0.01);
        if var_index <= 65
            lower_limit = 1;
            upper_limit = 131;
            me_comp(nn, 5) = mean(me_var_vals(lower_limit:upper_limit));

        elseif var_index > length(me_var_data) - 65
            lower_limit = length(me_var_data) - 131;
            upper_limit = length(me_var_data);
            me_comp(nn, 5) = mean(me_var_vals(lower_limit:upper_limit));
        else
            lower_limit = var_index - 65;
            upper_limit = var_index + 65;
            me_comp(nn, 5) = mean(me_var_vals(lower_limit:upper_limit));
        end
        
end 


csvwrite('azcompdata.csv', az_comp)
csvwrite('crcompdata.csv', cr_comp)
csvwrite('mecompdata.csv', me_comp)

total_data = [az_comp; cr_comp; me_comp];
csvwrite('compdata.csv', total_data)

