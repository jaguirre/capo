PREFIX='..'
EVEN_FILES=${PREFIX}'/even/sep0,2/frf_new/*I.uvGAL'
ODD_FILES=${PREFIX}'/odd/sep0,2/frf_new/*I.uvGAL'
EVEN_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${EVEN_FILES[@]}`
ODD_FILES=`lst_select.py -C psa6622_v003 --ra=4_10 ${ODD_FILES[@]}`
SEP='0,2'

for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,2,40)))"` ; do
    mkdir inject_sep${SEP}_${inject}
    echo SIGNAL_LEVEL=${inject}
    
    #~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject} ${EVEN_FILES} ${ODD_FILES}
    #echo "~/capo/pspec_pipeline/pspec_cov_v003_sigloss.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject}" ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt
    
    #noise only
    ~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject} --noise_only ${EVEN_FILES} ${ODD_FILES}
    echo "~/capo/pspec_pipeline/sigloss_sim.py --window=none -a cross -p I -c 110_130 -C psa6622_v003 -b 20 -i ${inject} --noise_only" ${EVEN_FILES} ${ODD_FILES} > inject_sep${SEP}_${inject}/notes.txt

    mv *bootsigloss*.npz inject_sep${SEP}_${inject}/.
done