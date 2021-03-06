#! /bin/tcsh -f

set MDL    = `echo $1 | cut -d "/" -f 1`
set EVENT  = `echo $1 | cut -d "/" -f 2`
set STREAM = `echo $1 | cut -d "/" -f 3`

set FILE = $2
      
set FLAG = $3
set CNT  = $4

set NAME   = `basename $FILE .mdst | cut -d "-" -f 5`


set OUT   = 'script_single/'$EVENT'/'$FLAG'_'$EVENT'_'$NAME'.sh'
set LOG   =    'log_single/'$EVENT'/'$FLAG'_'$EVENT'_'$NAME'.log'
set HBK   =    'hbk_single/'$EVENT'/'$FLAG'_'$EVENT'_'$NAME'.hbk'
set MDST  =   'mdst_single/'$EVENT'/'$FLAG'_'$NAME'.mdst'
#set MDST  =  'index_single/'$EVENT'/'$FLAG'_'$NAME'.index'

mkdir -p 'script_single/'$EVENT
#****************************************
      
echo '#! /bin/tcsh -f\n' > $OUT

echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/skim_dstrtaunu_single' >> $OUT
echo 'setenv BELLE_MESSAGE_LEVEL INFO' >> $OUT
echo 'setenv USE_GRAND_REPROCESS_DATA 1' >> $OUT
echo 'setenv FPDA_TIMEOUT 0' >> $OUT
echo 'source /sw/belle/local/etc/cshrc_general\n' >> $OUT # KEKCC
#echo 'source /belle/local/etc/cshrc_general\n' >> $OUT # NAGOYA

echo 'echo    $LD_LIBRARY_PATH'     >> $OUT
echo 'echo    $BASF_MODULE_DIR'     >> $OUT
echo 'echo    $PANTHER_TABLE_DIR\n' >> $OUT

echo "mkdir -p    log_single/"$EVENT >> $OUT
echo "mkdir -p    hbk_single/"$EVENT >> $OUT
echo "mkdir -p   mdst_single/"$EVENT >> $OUT
#echo "mkdir -p  index_single/"$EVENT >> $OUT

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

echo 'module put_parameter ' $MDL 'flag_SkimFile\\1'    >> $OUT
echo 'module put_parameter ' $MDL 'SkimFileName\\'$MDST >> $OUT
echo 'module put_parameter ' $MDL 'flag_single\\1\n'    >> $OUT

echo '' >> $OUT

echo 'initialize' >> $OUT
echo 'histogram define '$HBK'\n' >> $OUT

echo 'table save mdst_all'        >> $OUT
echo 'table save evtcls_all'      >> $OUT
echo 'table save evtvtx_all'      >> $OUT
echo 'table save gsim_rand'       >> $OUT
echo 'table save hepevt_all'      >> $OUT
echo 'table save mctype'          >> $OUT
echo 'table save level4_all'      >> $OUT
echo 'table save bgtbl_info'      >> $OUT
echo 'table save dattof_trgl0'    >> $OUT
echo 'table save reccdc_timing\n' >> $OUT

echo 'process_event '$FILE >> $OUT

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm '$HBK >> $OUT

echo $FILE

chmod 755 $OUT
