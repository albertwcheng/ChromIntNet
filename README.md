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

## Output format
[Cytoscape](https://cytoscape.org/) compatible SIF file and edge attributes (tab delimited)
In the <outputFolder>/
SIF file clusterNet.sif
```
<pref>_<N1> ChIAPETLoop <pref>_<N2>
```
Edge attribute file clusterNet.edgeAttrs.txt where interaction strength = sum of scores for all ChIA-PET pairs with tags in each cluster.
```
<pref>_<N1> (ChIAPETLoop) <pref>_<N2> <interaction_strength>
```
Interval file for the clusters clusters.bed where numTags is the number of tags clustered together
```
<chrom> <start> <end> <pref>_<N> <numTags>
```
Tag interval file nodes.bed where name = <pref>_<N>/<Name> where Name is from the original ChIA-PET input

```
<chrom> <start> <end> <pref>_<N>/<Name>
```