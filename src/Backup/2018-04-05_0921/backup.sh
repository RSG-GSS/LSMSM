##### User: Hirotaka Miura (B1HXM10).               
#####	Position: Research Analytics Associate.         
#####	Organization: Federal Reserve Bank of New York.
##### 07/30/2017: Modified.
#####	07/25/2017: Previously modified.
#####	07/25/2017: Created.
#####	Description: 
#####		- Script to backup files.
#####	Modifications:
#####		07/30/2017:
#####			- Add command to copy files without extension.

#!/bin/bash

##### Retrieve and set datetime variable.
#datetime=`date +%Y-%m-%d_%H%M%S`
datetime=`date +%Y-%m-%d_%H%M`
##### Set destination folder path.
dpath="Backup/$datetime"
##### Generate destination folder.
mkdir $dpath
##### Backup files with extension.
cp -p *.* $dpath
##### Backup files without extension.
##### See https://unix.stackexchange.com/questions/213674/copy-all-files-that-have-no-extension.
shopt -s extglob
cp -p !(*.*) $dpath
