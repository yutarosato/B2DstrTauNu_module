#! /bin/tcsh -f

setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu
setenv BELLE_MESSAGE_LEVEL INFO
setenv USE_GRAND_REPROCESS_DATA 1
source /sw/belle/local/etc/cshrc_general

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR

basf << EOF >  test.log

nprocess set 0

path create main
path create analysis

path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis DSTRTAUNU
path add_condition analysis >:0:EXIT
path add_condition analysis <=:0:KILL


initialize
histogram define test.hbk

process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/sig_kpi/mdst_setA/evtgen_exp_31__dtokpi_bp_sig__A_-0.mdst 100
process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/sig_kpi/mdst_setA/evtgen_exp_31__dtokpi_b0_sig__A_-0.mdst 100
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-uds&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-charm&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-mixed&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-charged&dt=on_resonance&bl=caseB&dv=zfserv&st=0

terminate
EOF

h2root test.hbk && rm test.hbk

