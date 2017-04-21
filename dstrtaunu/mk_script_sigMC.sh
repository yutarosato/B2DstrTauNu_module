#! /bin/tcsh -f

set DIR    = $1
set MDL    = `echo $2 | cut -d "/" -f 1`
set TAG    = `echo $2 | cut -d "/" -f 2`
set EXP   = $3
set FLAG   = $4
set MODE   = $5
      
set  OUT = 'script/'$MODE'/'$FLAG'_'$MODE'_'$EXP$TAG'.sh'
set  LOG =    'log/'$MODE'/'$FLAG'_'$MODE'_'$EXP$TAG'.log'
set  HBK =    'hbk/'$MODE'/'$FLAG'_'$MODE'_'$EXP$TAG'.hbk'
#set MDST =   'mdst/'$MODE'/'$FLAG'_'$MODE'_'$EXP$TAG'.mdst'
     
mkdir -p 'script/'$MODE
#****************************************
      
echo '#! /bin/tcsh -f\n' > $OUT

echo 'setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu' >> $OUT
echo 'setenv BELLE_MESSAGE_LEVEL INFO' >> $OUT
echo 'setenv USE_GRAND_REPROCESS_DATA 1' >> $OUT
echo 'setenv FPDA_TIMEOUT 0' >> $OUT
echo 'source /sw/belle/local/etc/cshrc_general\n' >> $OUT # KEKCC
#echo 'source /belle/local/etc/cshrc_general\n' >> $OUT # NAGOYA


echo 'echo    $LD_LIBRARY_PATH'     >> $OUT
echo 'echo    $BASF_MODULE_DIR'     >> $OUT
echo 'echo    $PANTHER_TABLE_DIR\n' >> $OUT

echo "mkdir -p    log/"$MODE >> $OUT
echo "mkdir -p    hbk/"$MODE >> $OUT
#echo "mkdir -p   mdst/"$MODE >> $OUT

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

#echo 'module put_parameter ' $MDL 'flag_SkimFile\\1'        >> $OUT
#echo 'module put_parameter ' $MDL 'SkimFileName\\'$MDST'\n' >> $OUT

echo '' >> $OUT

echo 'initialize' >> $OUT
echo 'histogram define '$HBK'\n' >> $OUT

set CNT = 0
foreach f ($DIR"evtgen_exp_"$EXP"_"*.mdst)
#foreach f ($DIR"evtgen_exp_"$EXP"_"*b0*.mdst) # onlys B0B0bar events are used for sig_kpi/nor_kpi with large samples(set??)
echo 'process_event '$f >> $OUT
@ CNT++
end

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm '$HBK >> $OUT 

echo $CNT $OUT

chmod 755 $OUT