# collectivedetection

Contains codes to produce the results in the paper:
Davidson, J.D., Sosna, M.M.G., Twomey, C.R., Sridhar, V.H., Leblanc, S.P., Couzin, I.D., 2021. Collective detection based on visual information in animal groups. Journal of The Royal Society Interface 18, 20210142. https://doi.org/10.1098/rsif.2021.0142

Associated data is available at:  
https://datadryad.org/stash/share/UH5t3WTf-PJPjzxmgDk1V2B7HwcsJiz1Nk6UVrSuO6A

## Note:  Only the files '4 - Data figures.ipynb' and '5 - Model.ipynb' are needed to reproduce results in the paper, if the simplified data and processed results are used (zip file saved_data_and_results.zip in the Dryad archive). 

# Files to reproduce figures and results in the paper
## 4 - Data figures.ipynb
Code to make all figures in the paper except Fig 7 (which is for the model)

## 5 - Model.ipynb
Code to fit model and make Figure 7, which is associated with the model and data comparison

## Analytical Approx.nb
Mathematica code to do the series expansion for the analytical approximation to model detection coverage

## detectionfunctions.py, raycastingfunctions.py
Contain functions used in other files


# Files to process data, calculate detection results, and save in a simpler format:
## 1 - Process data.ipynb
Reads in .h5 data files and saves simpler format '\*basicdata-v3.pklz', which can then be read in with, for example:
```python
[positions,orientations,tailcoords,lefteye,righteye,groupcentroid,groupheading,groupheadingxy,groupstates,pixelwidth,pixelheight,rotcoords] = pickle.load(gzip.open('10-fish/0066/basicdata-v3.pklz','rb'))
```

## 2 - Detection-run simulation.ipynb
Runs collective detection algorithm and saves for each trial, saving in the format (for example):
```python
pickle.dump([pointdist,allseen,allseen_high,allseen_low,allrelangle,processing_skipvalue],gzip.open('10-fish/0066/detection-probabilistic-perfish.pklz',"wb"),protocol=4)
```

## 3 - Process Detection results.ipynb
Processes detection and tracking results and saves in condensed form.  This saves the 4 files found in saved_data_and_results.zip in the Dryad archive, which can be read in as follows:
```python
[grid_allseen, grid_allseen_problow, grid_allseen_probhigh, grid_allseen_blind,processing_skipvalue] = pickle.load(gzip.open('detectionresults.pklz','rb'))   
[grid_frontbackdist,grid_sidesidedist,grid_groupstates,grid_groupheading,grid_orientations,grid_tailcoords,grid_positions,grid_groupcentroid,grid_lefteye,grid_righteye] = pickle.load(gzip.open('data-grid.pklz','rb'))
[degreebins,distbins,d1D,exampledata] = pickle.load(open('distributions+exampledata.pkl','rb'))
[indiv_mean,indiv_std,group_mean,group_std] = pickle.load(open('ind-groupmeans.pkl','rb'))
```
