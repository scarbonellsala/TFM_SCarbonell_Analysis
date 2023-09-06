#!/bin/bash -ex

# -C for color
# -M nim number of counts to include
# -P color palete file


function sushi() {
    local coords=$1
    local min=$2
    local gene=$3
    ~/bin/ggsashimi/ggsashimi.py -b ggsashimi.bams.LRonly.tsv -c ${coords}  -P ggsashimi.colorPalete.tsv -C 1 -M $min -F pdf -o ${gene}.${coords}.${min}.LRonly &
    ~/bin/ggsashimi/ggsashimi.py -b ggsashimi.bams.SRonly.tsv -c ${coords}  -P ggsashimi.colorPalete.tsv -C 1 -M $min -F pdf -o ${gene}.${coords}.${min}.SRonly  &
    ~/bin/ggsashimi/ggsashimi.py -b ggsashimi.bams.LRandSRcombined.tsv -c ${coords}  -P ggsashimi.colorPalete6.tsv -C 1 -M $min -F pdf -o ${gene}.${coords}.${min}.LRandSRcombined &
   wait
}


sushi chr10:112431442-112435000 6 ZDHHC6

sushi chr10:133401590-133402148 6 MTG1

sushi chr10:133400724-133403014 5 MTG1

sushi chr10:133400724-133403014 10 MTG1

sushi chr10:133400000-133403014 10 MTG1

sushi chr10:133402000-133403014 10 MTG1

sushi chr10:133401000-133403014 10 MTG1

sushi chr10:133401090-133403014 10 MTG1

sushi chr10:133392157-133424520 10 MTG1

sushi chr10:133392157-133424520 20 MTG1