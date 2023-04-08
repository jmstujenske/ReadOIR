# Read and convert Olympus OIR files
These scripts are used for reading Olympus .oir files into matlab or converting to tiff.
This script was edited by J. M. Stujenske in June 2022 to increase the speed of the functions. This consists of minimal edits to the original code, but speeds up reading oir files substantially.
oir2tiff.m written in August 2022 to directly convert oirs to tiffs. Requires the tiff read/write repository at https://github.com/jmstujenske/Matlab_FastTiffReadWrite .

# oir2stdData
MATLAB script for converting Olympus .oir image files into MATLAB readable images and metadata. <br>
These scripts were tested at MATLAB 2018a for Windows and Linux
## Usage
Default usage <br>
```
[~,stdData]=oir2stdData(pathToFile); 
```
pathToFile: string for file location. <br>
 stdData is struct containing images and metadata <br>
 stdData.Image: cell array: each cell contains xyzt movie <br>
 stdData.Metadata: struct contains metadata of movie(s) <br>
 <i> caution: Opening files can need huge memory. </i><br>

If you want to save data as a series of small files (<1GB), use,  <br>
```
output_list=oir2stdData(pathToFile,0,0);
```

## Trouble shoot
If an error occurred, try,<br>
```
[~,stdData]=oir2stdData(pathToFile,1,1); <br>
```
or
```
output_list=oir2stdData(pathToFile,0,1); 
```
Setting the last flag as 1 make opening speed bit slower but sometimes resolve errors.
