# ChromIntNet
Converts ChIA-PET paired format to Cytoscape compatible SIF and edge attribute files

## Input format
Paired form from ChIA-PET
```
Name chrom1 start1 end1 chrom2 start2 end2 score
```

## Usage
```bash
mkdir <outputFolder>
python2.7 ChromIntNet.py --cluster-prefix <pref> <input> <outputFolder>/
```