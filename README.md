# SCIP

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