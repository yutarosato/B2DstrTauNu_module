#! /bin/tcsh -f

setenv MY_TOP_DIR $HOME/dstrtaunu/modules/skim_dstrtaunu_single
setenv BELLE_MESSAGE_LEVEL INFO
setenv USE_GRAND_REPROCESS_DATA 1
setenv FPDA_TIMEOUT 0
source /sw/belle/local/etc/cshrc_general

echo    $LD_LIBRARY_PATH
echo    $BASF_MODULE_DIR
echo    $PANTHER_TABLE_DIR

mkdir -p    log_single/mixed
mkdir -p    hbk_single/mixed
mkdir -p   mdst_single/mixed
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

module put_parameter  DSTRTAUNU flag_SkimFile\1
module put_parameter  DSTRTAUNU SkimFileName\test.mdst
module put_parameter  DSTRTAUNU flag_single\1


initialize
histogram define test.hbk

table save mdst_all
table save evtcls_all
table save evtvtx_all
table save gsim_rand
table save hepevt_all
table save mctype
table save level4_all
table save bgtbl_info
table save dattof_trgl0
table save reccdc_timing

process_event /group/belle/bdata_b/mcprod/dat/e000055/evtgen/mixed/00/all/0127/on_resonance/00/evtgen-mixed-00-all-e000055r000007-b20090127_0910.mdst 100

terminate
EOF

h2root test.hbk && rm test.hbk
