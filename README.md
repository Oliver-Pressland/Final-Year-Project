# Automated feature detection and classification of Electrocardiographic Time Series data - Final Year Project for Brunel University London

Myocardial infarction is often correlated to the rise of heart conditions that may potentially be minor in nature, but can be spotted using data called an Electrocardiogram. These are irregular heart rhythms known as arrhythmia. Fortunately, heart rhythms can be monitored and diagnosed using a graph called an Electrocardiogram (ECG). As a part of my final year project at Brunel University London, I created a machine learning pipeline that would extract key features from ECGs in order to develop a solution that classifies different heart conditions. This pipeline has been developed in the form of a Jupyter Notebook.

Using an open source dataset consisting of 17 heart conditions, a waterfall approach was used to develop software in Python to discover and extract significant features, in particular the height and width of ECG components. These features were then analysed, and then trained on several supervised learning models. This resulted in the best model achieving 97% accuracy, 83% sensitivity, 84% positive predictive value and 98% specificity.

The data used for this project is a Mendeley Dataset published by Pawel Plawiak (Plawiak 2017) which is under a CC BY 4.0 license (Creative Commons — Attribution 4.0 International — CC BY 4.0, n.d.). The dataset contains 17 conditions and 1000 fragments. 

The directory "notebooks/libraries" contains the algorithms written for the project.

Plawiak, P. (2017), ‘Ecg signals (1000 fragments)’.

Creative Commons — Attribution 4.0 International — CC BY 4.0, (n.d.). URL: https://creativecommons.org/licenses/by/4.0/
