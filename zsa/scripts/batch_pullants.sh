#$ -S /bin/bash
#$ -V
#$ -cwd

ARGS=`pull_args.py $*`
ANT_LIST=49,41,47,19,29,28,34,51,10,3,25,48,24,55,27,57,9,58,1,4,17,13,56,59,22,61,35,18,5,32,30,23
echo pull_antpols.py -p xx -a "($ANT_LIST)_($ANT_LIST)" $ARGS
pull_antpols.py -p xx -a "($ANT_LIST)_($ANT_LIST)" $ARGS
