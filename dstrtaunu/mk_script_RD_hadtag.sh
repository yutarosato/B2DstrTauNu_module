#! /bin/tcsh -f

set MDL  = $1
set FILE = $2
set FLAG = $3
set CNT  = $4

set NAME   = `basename $FILE .mdst | cut -d "_" -f 2`

set OUT   = 'script_hadtag/rd/'$FLAG'_'$NAME'.sh'
set LOG   =    'log_hadtag/rd/'$FLAG'_'$NAME'.log'
set HBK   =    'hbk_hadtag/rd/'$FLAG'_'$NAME'.hbk'

mkdir -p 'script_hadtag/rd'
#****************************************
      
echo '#! /bin/tcsh -f\n' > $OUT

echo "setenv EKPTURBO_DIR /group/belle/ekp-turbo/data # for KEKCC"                              >> $OUT # FULLRECON
echo "#setenv EKPTURBO_DIR /panfs/sgt1/bfactory/belle_data/caseB/ekp-turbo/data # for NAGOYA\n" >> $OUT # FULLRECON

echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu' >> $OUT
echo 'setenv BELLE_MESSAGE_LEVEL INFO' >> $OUT
echo 'setenv USE_GRAND_REPROCESS_DATA 1' >> $OUT
echo 'setenv FPDA_TIMEOUT 0' >> $OUT
echo 'source /sw/belle/local/etc/cshrc_general\n' >> $OUT # KEKCC
#echo 'source /belle/local/etc/cshrc_general\n' >> $OUT # NAGOYA


echo 'echo    $LD_LIBRARY_PATH'     >> $OUT
echo 'echo    $BASF_MODULE_DIR'     >> $OUT
echo 'echo    $PANTHER_TABLE_DIR\n' >> $OUT

echo "mkdir -p    log_hadtag/rd" >> $OUT
echo "mkdir -p    hbk_hadtag/rd" >> $OUT

echo 'basf << EOF > ' $LOG'\n' >> $OUT

echo 'nprocess set 0\n' >> $OUT

echo 'path create main'       >> $OUT
echo 'path create brecon'     >> $OUT
echo 'path create analysis\n' >> $OUT

echo 'path add_module main fix_mdst'           >> $OUT
echo 'path add_condition main >:0:brecon'      >> $OUT
echo 'path add_condition main <=:0:KILL'       >> $OUT
echo 'path add_module brecon ekpturbo'         >> $OUT
echo 'path add_condition brecon >:0:analysis'  >> $OUT
echo 'path add_condition brecon <=:0:KILL'     >> $OUT
echo 'path add_module analysis' $MDL           >> $OUT
echo 'path add_condition analysis >:0:EXIT'    >> $OUT
echo 'path add_condition analysis <=:0:KILL\n' >> $OUT

echo 'module put_parameter ' $MDL 'flag_hadtag\\1\n'    >> $OUT # FULLRECON

echo 'initialize' >> $OUT
echo 'histogram define '$HBK'\n' >> $OUT

echo 'process_event '$FILE >> $OUT

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm '$HBK >> $OUT

echo $FILE

chmod 755 $OUT
