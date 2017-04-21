#! /bin/tcsh -f

set DIR    = $1
set MDL    = `echo $2 | cut -d "/" -f 1`
set PAR    = `echo $2 | cut -d "/" -f 2`
set STREAM = `echo $2 | cut -d "/" -f 3`
set EXP    = $3
set FLAG   = $4
      
#set  OUT = "script/${PAR}/${FLAG}_${PAR}_${EXP}_s${STREAM}.sh"
#set  LOG =    "log/${PAR}/${FLAG}_${PAR}_${EXP}_s${STREAM}.log"
#set  HBK =    "hbk/${PAR}/${FLAG}_${PAR}_${EXP}_s${STREAM}.hbk"

set  OUT = "script/${PAR}/${FLAG}_${PAR}_${EXP}_s0${STREAM}.sh"
set  LOG =    "log/${PAR}/${FLAG}_${PAR}_${EXP}_s0${STREAM}.log"
set  HBK =    "hbk/${PAR}/${FLAG}_${PAR}_${EXP}_s0${STREAM}.hbk"
     
mkdir -p "script/${PAR}"
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

echo "mkdir -p    log/${PAR}" >> $OUT
echo "mkdir -p    hbk/${PAR}" >> $OUT

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

echo 'module put_parameter ' $MDL 'flag_DststMC\\2' >> $OUT

echo '' >> $OUT

echo 'initialize' >> $OUT
echo 'histogram define '$HBK'\n' >> $OUT

set CNT = 0
foreach f (${DIR}"/evtgen_exp_"${EXP}"_"*".mdst")
echo 'process_event '$f >> $OUT
@ CNT++
end

echo '\nterminate' >> $OUT
echo 'EOF\n' >> $OUT

echo 'h2root '$HBK '&& rm '$HBK >> $OUT 

echo $CNT $OUT

chmod 755 $OUT
