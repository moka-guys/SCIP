# SCIP

## About
The SCIP python app was developed for the Sickle Cell in Pregnancy (SCIP) project, which aims to provide NIPD for pregnancies that may be affected by Sickle Cell Disease or Beta Thalassemia. 

The SCIP python app consists of the main script scip.py which coordinates the entire analysis and several accessory scripts containing supportive functions. The SCIP python app has three inputs: (1) The HBB mpileup file which contains count data on 208 highly polymorphic SNPs which are used to estimate the paternal and therefore foetal fraction of DNA in the sample. 

With these inputs, the SCIP python app identifies paternal 'informative' alleles, calculates the median foetal fraction, determines the maternal genetic background for each SCD allele and uses these information together to calculate the upper and lower boundaries for the modified sequential probability ratio test analysis, predict the foetal genotype, and output a report containing the sample name, predictions for the HbS, HbC, HbD and HbE alleles, an overall genotype, informative SNP information and interactive graphs (see figure 8). 

## Usage
`python ./code/scip.py [output html file] [mpileup file path] [sickle cell alleles mpileup file path]`

## Docker image build
`make build`

To build the docker image from scratch, using no cached layers:

`make cleanbuild` 

## Docker image upload
To upload the built docker image to DockerHub:

`make push`

## Test
This runs against the two files with all or partially missing positions

`make test`