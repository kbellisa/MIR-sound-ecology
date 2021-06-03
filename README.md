# MIR-sound-ecology

This research was conducted for a dissertation. The aim was to extract features using conceptual ideas from Music Information Retrieval systems for a sound ecology problem.

This project requires use of R, sonic batch annotator, and access to Sonic Visualiser.  It is necessary to install vamp plug-in files (https://www.vamp-plugins.org/download.html) to extract the features using Unix command line after installation of sonic batch annotator (https://vamp-plugins.org/sonic-annotator/).

Once you have selected the features and have the vamp plug-in, then you can perform this task in R. This study did not exclude outliers but results would be greatly improved by removing outliers [sounds that were extreme like thunder]. As real world data contains a mix of unpredictable loud sounds, we used all measurements.

Music Information Retrieval Project:  ORDER OF SCRIPTS

1. METRICSsoundBank.R (retrieve site population) export Rdata [note: or use any acoustic data source]

2. RandomSample_UndisturbedCollection.R (retrieve undisturbed sites from site population from within a restricted database) [note: this may only be performed with credentials -- substitute with any acoustic data source]

3. FeatureDataSets subfolder EXTRACT FEATURES IN UNIX (use Terminal on Mac) Run Sonic Annotator Scripts to generate feature extractions on random sample all feature results zipped

4. TrainingPhase folder Metrics_FeaturesSummary_Final.R (generate matrix of all features, select subsample from all features with equal proportion for each feature, then kruskal test to compare subsample to sample)

5. LDA subfolder discriminant analysis / ground truth model

6. WindowObservationLengths [note: this allows to compare lengths of recordings to select optimal duration]

7. AutomatedCSV.m
CorrelationPlot.R — SEM and Merge performed in MATLAB — Observation mean values for each window length period performed in MATLAB
read values.csv [note: this series of steps uses Matlab to generate comparison of recording lengths - once length is determined, it is a fixed value in the code]

8. TESTING PHASE TestingPhase folder test_data_v2 = predictor set TestFeatures_180 = all features
Discriminant Analysis on data removing outliers test_data_MDS_LD.R

9. Quadratic Discriminant Analysis on data including outliers in QDA subfolder SoundMetric_QDA.R

