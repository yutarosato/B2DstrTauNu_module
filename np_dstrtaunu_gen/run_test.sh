#! /bin/tcsh -f

setenv MY_TOP_DIR $HOME/dstrtaunu/modules/np_dstrtaunu_gen
setenv BELLE_MESSAGE_LEVEL INFO
setenv USE_GRAND_REPROCESS_DATA 1
setenv FPDA_TIMEOUT 0
source /sw/belle/local/etc/cshrc_general

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR

mkdir -p    log/0p0
mkdir -p    hbk/0p0
basf << EOF >  test.log

nprocess set 0

path create main
path create analysis

path add_module main DSTRTAUNU
path add_condition main >:0:EXIT
path add_condition main <=:0:KILL


initialize
histogram define test.hbk

#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np/evtgen/mdst/0p0/evtgen_exp_31__BSemiTauonic_BSTD_2HDMTYPE2-0p0_Charged_Dstaunu__A_-0.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np/evtgen/mdst/0p0/evtgen_exp_31__BSemiTauonic_BSTD_2HDMTYPE2-0p0_Charged_Dtaunu__A_-0.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np/evtgen/mdst/0p0/evtgen_exp_31__BSemiTauonic_BSTD_2HDMTYPE2-0p0_Mixed_Dstaunu__A_-0.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np/evtgen/mdst/0p0/evtgen_exp_31__BSemiTauonic_BSTD_2HDMTYPE2-0p0_Mixed_Dtaunu__A_-0.mdst
process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np/evtgen/mdst_double/0p0/evtgen_exp_31__BSemiTauonic_BSTD_2HDMTYPE2-0p0_Mixed_Dstaunu_DoubleSemileptonicDecay__A_-0.mdst 1000
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np/evtgen/mdst_double/0p0/evtgen_exp_31__BSemiTauonic_BSTD_2HDMTYPE2-0p0_Charged_Dtaunu_DoubleSemileptonicDecay__A_-0.mdst 1000

terminate
EOF

h2root test.hbk && rm test.hbk
