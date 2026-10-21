# Analysis code: online vs offline Mandarin tone learning

This folder holds the data-processing, statistics and plotting scripts behind the results reported in the manuscript and the Supplementary Materials. It contains no participant data.

Scripts are numbered `s01`–`s28` in the order they should be run. MATLAB script names must start with a letter, so each number is prefixed with `s`.

## Software

- MATLAB with EEGLAB 2024.0 (EEG preprocessing, alpha power, topography) and the Statistics and Machine Learning Toolbox (`fitlme`, `fitglme`)
- Python 3 with NeuroKit2 (ECG/EDA processing), pandas, numpy and matplotlib
- Praat (audio segmentation and normalisation)

The scripts use the original absolute paths, for example `D:\online_learning\...`. Change these paths at the top of each script before running it.

## Scripts

### 01_EEG_preprocessing

| Script | Purpose | Main output |
|---|---|---|
| `s01_preprocess_training_EEG.m` | Training-session EEG preprocessing up to ICA | `*_ica_clean.set` |
| `s02_preprocess_resting_EEG.m` | Resting-state EEG preprocessing | `*_rest_ica_clean.set` |
| `s03_detect_vocal_onset_accelerometer.m` | Vocal onset/offset detection from the neck accelerometer | `*_with_vocal_events.set` |
| `s04_compute_training_alpha_power.m` | Posterior alpha power 0.3–2.3 s after each production offset | `alpha_VOFF_block_level_0p3_2p3.csv` |
| `s05_compute_resting_alpha_baseline.m` | Resting-state alpha baseline | `alpha_rest_baseline_2s.csv` |

### 02_alpha_statistics (main text: alpha results, Tables 1–2)

| Script | Purpose |
|---|---|
| `s06_alpha_LMM_simple_effects.m` | LMM `AlphaPower ~ Group*Block + baseline` with random-effects comparison; omnibus tests, adjusted means, group contrasts within each block, block contrasts within each group |
| `s07_alpha_LMM_subject_points.m` | The same model; exports the subject-level points used in Figure 2 |
| `s08_alpha_change_contrasts_effect_sizes.m` | Between-group differences in change (Table 1), within-group B3/B1 ratios (Table 2), approximate partial η², Holm correction |

### 03_alpha_topography (Figure 2, scalp map)

| Script | Purpose |
|---|---|
| `s09_alpha_topography_channel_averages.m` | Channel-level alpha averaged within participants (rest and training blocks), then across participants |
| `s10_alpha_topography_plot.m` | Scalp map with the posterior ROI (`alpha_topography_only_panel_c_600dpi.tif`) |

### 04_production_rating

| Script | Purpose |
|---|---|
| `s11_split_production_recordings.praat` | Splits the long production recordings |
| `s12_segment_syllables.praat` | Extracts single-syllable files |
| `s13_normalize_audio_70dB.praat` | Intensity normalisation to 70 dB before rating |
| `s14_merge_rater_sheets.m` | Merges each rater's rating sheets |
| `s15_interrater_agreement_fleiss_kappa.m` | Fleiss' κ across the three raters |
| `s16_final_tone_by_majority.m` | Final tone category by majority agreement (`all_rater_finalresult.xlsx`) |
| `s17_count_correct_tones.m` | Correct-tone counts per participant, session and tone (`tone_count_summary.xlsx`) |

### 05_production_statistics (main text: production results, Table 3)

| Script | Purpose |
|---|---|
| `s18_production_LMM_simple_effects.m` | LMM on empirical-logit production accuracy (`group*session + tone`), with adjusted means and simple effects |
| `s19_production_change_difference_check.m` | Earlier version of the change-difference contrasts; `s20` reads its output only to check that it reproduces the earlier estimates |
| `s20_production_change_contrasts_holm.m` | All six between-group differences in change (Holm-corrected, Table 3) and within-group Post−Pre changes |

### 06_identification (main text: tone identification and reaction time)

| Script | Purpose |
|---|---|
| `s21_identification_extract_and_model_selection.m` | Trial-level data extraction, accuracy GLMM and reaction-time LMM, random-effects comparison |
| `s22_identification_final_GLMM.m` | Final accuracy GLMM, Type III tests and predicted probabilities |

### 07_autonomic_supplementary (Supplementary S1)

| Script | Purpose |
|---|---|
| `s23_ECG_EDA_preprocessing_neurokit.ipynb` | ECG and EDA processing with NeuroKit2 (heart rate, RMSSD, SCL, SCR rate) and quality-control flags; cell outputs cleared |
| `s24_EDA_distribution_QC.py` | Distribution and quality-control checks for the EDA measures |
| `s25_ECG_LMM_HR_logRMSSD.m` | LMMs for heart rate and natural-log RMSSD |
| `s26_EDA_LMM_SCL_logSCRrate.m` | LMMs for SCL and natural-log SCR rate |

### 08_questionnaires_supplementary (Supplementary S4)

| Script | Purpose |
|---|---|
| `s27_questionnaire_analysis.m` | Group descriptives and comparisons for the background scales and the post-training ratings |

### 09_figures

| Script | Purpose |
|---|---|
| `s28_plot_figure2.py` | Figure 2: production and alpha bar charts combined with the scalp map from `s10`. It reads the outputs of `s07`, `s08`, `s10`, `s18` and `s20` |

Figure 1 (the experimental procedure) is a schematic drawing and has no script.

# Experiment materials: online vs offline Mandarin tone learning

This folder holds the PsychoPy experiment files and the stimulus materials used in the study. It contains no participant data.

The experiments were built with PsychoPy 2022.2.4. Every file sits in `psychopy_experiment/` with its original name, because the experiment files refer to stimuli by relative path.

## Experiment files

Each `.psyexp` file is the PsychoPy Builder source. The matching `_lastrun.py` file is the Python script PsychoPy generated from it.

| File | Session part | Condition file used |
|---|---|---|
| `301_resting_state` | Eyes-open/eyes-closed resting state and tone instruction | none |
| `302_1_pretest1` | Pre-test, part 1: attentive listening (Supplementary S3) | `ERPstimuli.xlsx` |
| `302_2_pretest2` | Pre-test, part 2: tone instruction, tone identification and tone production tests | `audio.xlsx`, `loopTemplate1.xlsx` |
| `303_tutor` | Listen-and-repeat training with the tutor, with interim production tests | `tutorpinyin.xlsx`, `loopTemplate1.xlsx` |
| `304_1_post_test1` | Post-test, part 1: attentive listening | `ERPstimuli.xlsx` |
| `304_2_post_test2` | Post-test, part 2: tone identification and tone production tests | `audio.xlsx`, `loopTemplate1.xlsx` |

## Stimulus and condition files

| File | Content |
|---|---|
| `tutorpinyin.xlsx` | The 63 training syllables (pinyin with tone marks) that the tutor produced in the listen-and-repeat training |
| `loopTemplate1.xlsx` | Items for the tone production tests (`TEST_LIST`) |
| `audio.xlsx` | The 72 tone identification trials: audio file and correct response (1 = rising, 2 = falling-rising, 3 = falling) |
| `audio_clip/` | The 72 identification stimuli: 6 syllables × 3 tones × 4 speakers. File names follow `syllable + tone (_speaker).wav`, e.g. `bo3_2.wav` |
| `ERPstimuli.xlsx`, `a1.wav`–`a11.wav` | Attentive-listening sequence (stimulus order, ISI, EEG trigger codes) and its sound files |
| `11.png`–`19.png`, `Tone_information_audio.wav` | Tone instruction slides (in Finnish) and their narration |
| `eye_open.wav`, `eye_closed.wav`, `pre_pretest.wav`, `question1–4.wav`, `test1.wav`, `test3.wav`, `after_pretest.wav`, `ready.wav`, `Sametest.wav`, `interactive.wav`, `non_interactive.wav` | Spoken task instructions (in Finnish) |

The tutor training itself was delivered live over video, so there are no pre-recorded training videos.

