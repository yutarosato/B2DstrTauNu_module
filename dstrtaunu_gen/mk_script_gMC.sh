#! /bin/tcsh -f

set MDL    = `echo $1 | cut -d "/" -f 1`
set EVENT  = `echo $1 | cut -d "/" -f 2`
set STREAM = `echo $1 | cut -d "/" -f 3`
set TAG    = `echo $1 | cut -d "/" -f 4`

set expNo  = `echo $2 | cut -d "/" -f 1`
set runNoi = `echo $2 | cut -d "/" -f 2`
set runNof = `echo $2 | cut -d "/" -f 3`
      
set FLAG = $3
set RUN  = $4
set CNT  = $5

set NAME = `printf "e%03dc%03ds%02d" $expNo $CNT $STREAM`

set OUT   = 'script/'$EVENT'/'$FLAG'_'$EVENT'_'$NAME$TAG'.sh'
set LOG   =    'log/'$EVENT'/'$FLAG'_'$EVENT'_'$NAME$TAG'.log'
set HBK   =    'hbk/'$EVENT'/'$FLAG'_'$EVENT'_'$NAME$TAG'.hbk'

mkdir -p 'script/'$EVENT
#****************************************
      
echo '#! /bin/tcsh -f\n' > $OUT

echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu_gen' >> $OUT
echo 'setenv BELLE_MESSAGE_LEVEL INFO' >> $OUT
echo 'setenv USE_GRAND_REPROCESS_DATA 1' >> $OUT
echo 'source /sw/belle/local/etc/cshrc_general\n' >> $OUT

echo 'echo    $LD_LIBRARY_PATH'     >> $OUT
echo 'echo    $BASF_MODULE_DIR'     >> $OUT
echo 'echo    $PANTHER_TABLE_DIR\n' >> $OUT

echo "mkdir -p    log/"$EVENT >> $OUT
echo "mkdir -p    hbk/"$EVENT >> $OUT

echo 'basf << EOF > ' $LOG'\n' >> $OUT

echo 'nprocess set 0\n' >> $OUT

echo 'path create main'     >> $OUT
echo 'path create analysis\n' >> $OUT

echo 'path add_module main fix_mdst'         >> $OUT
echo 'path add_condition main >:0:analysis'  >> $OUT
echo 'path add_condition main <=:0:KILL'     >> $OUT
echo 'path add_module analysis' $MDL         >> $OUT
echo 'path add_condition analysis >:0:EXIT'  >> $OUT
echo 'path add_condition analysis <=:0:KILL\n' >> $OUT

echo 'initialize' >> $OUT
echo 'histogram define '$HBK'\n' >> $OUT

if( $expNo < 30 ) then
  @ STREAM += 10
endif

~syutaro/script/mk_data_list_gmc.sh $expNo $runNoi $runNof $EVENT $STREAM >> $OUT

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm ' $HBK >> $OUT

echo $expNo $runNoi $runNof $CNT

chmod 755 $OUT

echo 'bsub -q l '$OUT >> $RUN

