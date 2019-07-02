# DeepLife_16S
python code for analyzing diversity within the Census of Deep Life (CoDL) 16S rRNA amplicon dataset

## Authors
Kenneth J. Locey, Diné College

## Files
**radMetrics.py**: Reads .biom files and generates an output file of the following metrics for each sample:

* Evenness:  
	* Variance of the sample abundances  
	* Smith and Wilson's evenness index (Evar)    
	* Simpson's evenness  
	* EQ evenness  
	* O evenness    
	* Camargo's evenness  
	* Nee's evenness
	* Pielou's evenness  
	* Heip's evenness  

* Dominance
	* Berger-Parker
	* Simpson's dominance
	* Maximum abundance
	* McNaughton's dominance

* Rarity
	* skewnness of the abundance distribution
	* log-modulo skewness of the abundance distribution 
	* percent of abundances that are singletons
	* percent of abundances with less than 0.1% of total sample abundance

* Richness
	* Preston's alpha and richness, from Curtis and Sloan (2002).
	* Chao1
	* ACE
	* Jackknife1
	* Jackknife2
	* Margalef's index
	* Menhinick's index

**
