#! /bin/tcsh -f

setenv MY_TOP_DIR $HOME/dstrtaunu/modules/dstrtaunu
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

path add_module main fix_mdst
path add_condition main >:0:analysis
path add_condition main <=:0:KILL
path add_module analysis DSTRTAUNU
path add_condition analysis >:0:EXIT
path add_condition analysis <=:0:KILL


initialize
histogram define test.hbk

process_event /home/belle/syutaro/Store/dstrtaunu/mcprod_np_nagoya/store/BSemiTauonic_BSTD_2HDMTYPE2-0p0_Mixed_Dstaunu/mdst/0//evtgen_exp_31_BSemiTauonic_BSTD_2HDMTYPE2-0p0_Mixed_Dstaunu-0.mdst 100

terminate
EOF

h2root test.hbk && rm test.hbk
