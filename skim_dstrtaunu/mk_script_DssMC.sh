#! /bin/tcsh -f

set MDL       = $1
set FLAG      = $2
set EXP       = $3
set CNT       = $4
set CNT_START = $5
set CNT_END   = $6
set list      = $7

set OUT  = 'script/'$FLAG'/'$FLAG'_e0'$EXP'_c'`printf "%03d" $CNT`'.sh'
set LOG  =    'log/'$FLAG'/'$FLAG'_e0'$EXP'_c'`printf "%03d" $CNT`'.log'
set HBK  =    'hbk/'$FLAG'/'$FLAG'_e0'$EXP'_c'`printf "%03d" $CNT`'.hbk'
set MDST =   'mdst/'$FLAG'/'$FLAG'_e0'$EXP'_c'`printf "%03d" $CNT`'.mdst'
#set MDST =  'index/'$FLAG'/'$FLAG'_e0'$EXP'_c'`printf "%03d" $CNT`'.index'
     
mkdir -p 'script/'$FLAG
#****************************************
      
echo '#! /bin/tcsh -f\n' > $OUT

echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/skim_dstrtaunu' >> $OUT
#echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu' >> $OUT
echo 'setenv BELLE_MESSAGE_LEVEL INFO'   >> $OUT
echo 'setenv USE_GRAND_REPROCESS_DATA 1' >> $OUT
echo 'setenv FPDA_TIMEOUT 0' >> $OUT
echo 'source /sw/belle/local/etc/cshrc_general\n' >> $OUT # KEKCC
#echo 'source /belle/local/etc/cshrc_general\n' >> $OUT # NAGOYA


echo 'echo    $LD_LIBRARY_PATH'     >> $OUT
echo 'echo    $BASF_MODULE_DIR'     >> $OUT
echo 'echo    $PANTHER_TABLE_DIR\n' >> $OUT

echo "mkdir -p    log/"$FLAG >> $OUT
echo "mkdir -p    hbk/"$FLAG >> $OUT
echo "mkdir -p   mdst/"$FLAG >> $OUT

echo 'basf << EOF > ' $LOG'\n' >> $OUT

echo 'nprocess set 0\n' >> $OUT
  
echo 'path create main'     >> $OUT
echo 'path create analysis\n' >> $OUT

echo 'path add_module main fix_mdst'           >> $OUT
echo 'path add_condition main >:0:analysis'    >> $OUT
echo 'path add_condition main <=:0:KILL'       >> $OUT
echo 'path add_module analysis' $MDL           >> $OUT
echo 'path add_condition analysis >:0:EXIT'    >> $OUT
echo 'path add_condition analysis <=:0:KILL\n' >> $OUT

echo 'module put_parameter ' $MDL 'flag_SkimFile\\1'        >> $OUT
echo 'module put_parameter ' $MDL 'SkimFileName\\'$MDST'\n' >> $OUT

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

sed -n "${CNT_START},${CNT_END}p" ${list} | awk '{print "process_event "$1}' >> $OUT

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm '$HBK >> $OUT 

echo "$EXP $CNT $CNT_START $CNT_END $OUT"

chmod 755 $OUT

