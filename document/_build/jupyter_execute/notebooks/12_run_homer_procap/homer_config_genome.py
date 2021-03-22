Hi Alex, I was following the HOMER analysis scripts from this tutorial http://homer.ucsd.edu/homer/ngs/csRNAseq/index.html

I am running into some problems and I was wondering if you could give me a hand.

At several points, the script requires HOMER to download the genome information from UCSC into the HOMER data folder. Therefore, I got the permission denied.

I did try to feed in the scripts with the genome data I downloaded myself, but the output results seems weird to me.

Is it possible for you to download the genome data using HOMER so that I can trouble shoot the problems I encountered by comparing the results I have with the results I will run using the default configuration in the HOMER document?

I have documented the scripts to download the genome info in this attached notebook.

### Note: currently there is only mouse genome under the homer data
!ls /data/reddylab/software/homer/data/genomes/

!ls /data/reddylab/software/homer/data/genomes/mm10

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
export PATH=/data/reddylab/software/homer/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/config_list.txt \
    <<'EOF'
#!/bin/bash

### this script is to check the available genome that can be downloaded using HOMER
perl /data/reddylab/software/homer/configureHomer.pl -list

EOF

!cat /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/config_list.txt

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/config_hg.txt \
    <<'EOF'
#!/bin/bash

### this script download the genome information under the folder
### /data/reddylab/software/homer/data/genomes/
### below I am trying to download the hg19 and hg38 genome information
perl /data/reddylab/software/homer/configureHomer.pl -install hg19
perl /data/reddylab/software/homer/configureHomer.pl -install hg38

EOF

!cat /data/reddylab/software/homer/config.txt







%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
export PATH=/data/reddylab/software/homer/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/config_hg.txt \
    <<'EOF'
#!/bin/bash

### this script is to check the available genome that can be downloaded using HOMER
perl /data/reddylab/software/homer/configureHomer.pl -list

### this script download the genome information under the folder
### /data/reddylab/software/homer/data/genomes/
### below I am trying to download the hg19 and hg38 genome information
perl /data/reddylab/software/homer/.//configureHomer.pl -install hg19
perl /data/reddylab/software/homer/.//configureHomer.pl -install hg38

EOF

### Note: currently there is only mouse genome under the homer data
!ls /data/reddylab/software/homer/data/genomes/
!ls /data/reddylab/software/homer/data/genomes/mm10
!ls -1 /data/reddylab/software/homer/data/genomes/hg38
!ls -1 /data/reddylab/software/homer/data/genomes/hg19

cat /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/config_hg.txt

cat /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/config_h38.txt

