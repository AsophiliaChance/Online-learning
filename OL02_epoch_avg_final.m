clear;clc;

basedir = {'D:\online_learning\equiprobability_EEG\pre';'D:\online_learning\equiprobability_EEG\post'};
filt ='*_ica_clean.set';

deviant_types = {'S  1';'S  4';'S  7';'S 10'};

art_thresh = 150; % +/- uV
twin = [-0.1 0.8];
savedir='D:\online_learning\equiprobability_EEG\results';
%%
eeglab;
close(gcf);
%%
for md=1:2% md1的25(OL27)有问题md2的25(OL27)30(OL33)有问题
    cd(basedir{md});
    outputdir = basedir{md};
    files = dir(filt);
    savenam=strrep(basedir{md}, '\', '_');
    savenam=strrep(savenam, 'D:_online_learning_equiprobability_EEG_', '');
%%
    for curfile =25%:length(files)
        file = files(curfile).name;
        EEG = pop_loadset(file,pwd);
        [pth,nam,ext] = fileparts(file);
        
        nam = nam(1:8);
        fprintf('Working on %s\n',[nam ext]);
        %% 重参考
        %EEG=pop_chanedit(EEG, 'append',1,'changefield',{2,'labels',''},'insert',2,'insert',2,'delete',2,'delete',2,'delete',2,'insert',2,'changefield',{2,'labels','Cz'},'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BESA\\GSN-HydroCel-129.sfp');
        EEG=pop_chanedit(EEG, 'lookup','S:\\Program\\matlab2019toolbox\\eeglab_current\\eeglab2024.0\\plugins\\dipfit\\standard_BESA\\GSN-HydroCel-129.sfp');
        EEG = pop_reref(EEG, [57 100]); 
      
%%        
        for idev=1:length(deviant_types)
            %%           
            [EEG1,indices0] = pop_epoch( EEG, { deviant_types{idev}  }, [-0.1         0.80], 'valuelim', [-art_thresh  art_thresh]);
            EEG1 = pop_rmbase( EEG1, [-100 0] ,[]);
          
            EEG1 = eeglab2fieldtrip(EEG1, 'preprocessing');

            
            %%          
            cfg = [];
            cfg.keeptrials = 'no';
            ERP{idev}   = ft_timelockanalysis(cfg, EEG1);
            
            cfg.keeptrials = 'yes';
            single_trial{idev}   = ft_timelockanalysis(cfg, EEG1);
                         
        end
         save([savedir,'\',nam,'_erp'],'ERP');
         save([savedir,'\',nam,'_singletrial'],'single_trial');
    end       
    end
  
