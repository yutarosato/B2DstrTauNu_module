#! /bin/tcsh -f

setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu
setenv BELLE_MESSAGE_LEVEL INFO
setenv USE_GRAND_REPROCESS_DATA 1
setenv FPDA_TIMEOUT 0
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

module put_parameter  DSTRTAUNU flag_single\1


initialize
histogram define test.hbk

process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000007.mdst 50
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000012.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000013.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000017.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000024.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000025.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000026.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_single/mdst_single/mdst15//mixed/gMC_e000055r000027.mdst

terminate
EOF

h2root test.hbk && rm test.hbk
