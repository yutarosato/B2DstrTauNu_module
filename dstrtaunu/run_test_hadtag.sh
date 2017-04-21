#! /bin/tcsh -f

setenv EKPTURBO_DIR /group/belle/ekp-turbo/mc_stream_0 # for KEKCC
#setenv EKPTURBO_DIR /panfs/sgt1/bfactory/belle_data/caseB/ekp-turbo/mc_stream_0 # for NAGOYA

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
path create brecon
path create analysis

path add_module main fix_mdst
path add_condition main >:0:brecon
path add_condition main <=:0:KILL
path add_module brecon ekpturbo
path add_condition brecon >:0:analysis
path add_condition brecon <=:0:KILL
path add_module analysis DSTRTAUNU
path add_condition analysis >:0:EXIT
path add_condition analysis <=:0:KILL

module put_parameter  DSTRTAUNU flag_hadtag\1

initialize
histogram define test.hbk

process_event /home/belle/syutaro/Store/dstrtaunu/modules/skim_dstrtaunu_hadtag/mdst_hadtag/mdst15//mixed/gMC_mixed_e031c001s00.mdst 100

terminate
EOF

h2root test.hbk && rm test.hbk
