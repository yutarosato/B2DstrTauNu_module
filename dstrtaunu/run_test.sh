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

initialize
histogram define test.hbk

#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/mixed/gMC_mixed_e027c004s00.mdst 0
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/mixed/gMC_mixed_e055c026s00.mdst 0
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/mixed/gMC_mixed_e063c013s00.mdst 0
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/mixed/gMC_mixed_e055c006s01.mdst 0
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/mixed/gMC_mixed_e037c054s02.mdst 0
process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu/mdst/mdst20/mixed/gMC_mixed_e037c053s02.mdst 100


terminate
EOF

h2root test.hbk && rm test.hbk
