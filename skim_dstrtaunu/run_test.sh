#! /bin/tcsh -f

setenv MY_TOP_DIR $HOME/dstrtaunu/modules/skim_dstrtaunu
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

module put_parameter DSTRTAUNU flag_SkimFile\1
module put_parameter DSTRTAUNU SkimFileName\test.mdst


initialize
histogram define test.hbk

process_event /hsm/belle/bdata2/users/huschle/DssMC/gsim/mdst/evtgen_exp_25_Dss-0.mdst 100000
#process_event /hsm/belle/bdata2/users/huschle/DssMC/gsim/mdst/evtgen_exp_25_Dss-1.mdst

#process_url http://bweb3/montecarlo.php?ex=27&rs=1507&re=1507&ty=evtgen-charged&dt=on_resonance&bl=caseB&dv=zfserv&st=10
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/nor_kpi_dstst/mdst_setAA/evtgen_exp_31__dststdtokpi_b0_nor__AA_-0.mdst 100
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/nor_kpi_dstpipi/mdst_setA/evtgen_exp_31__dstpipidtokpi_b0_nor__AA_-0.mdst 100

#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e013c008s00.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e013c008s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e013c008s02.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e013c008s03.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e013c008s04.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e013c008s05.mdst

#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c001s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c002s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c003s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c004s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c005s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c006s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c007s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c008s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c009s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c010s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c011s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c012s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c013s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c014s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c015s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c016s01.mdst
#process_event /home/belle/syutaro/Store/dstrtaunu/modules/dstrtaunu/mdst/mdst6/mixed/gMC_mixed_e031c017s01.mdst

#process_url http://bweb3/montecarlo.php?ex=31&rs=138&re=383&ty=evtgen-mixed&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_event /home/belle/syutaro/dstrtaunu/modules/skim_dstrtaunu/index/hz/mixed/gMC_mixed_e031c001s00.index

#process_url http://bweb3/montecarlo.php?ex=55&rs=40&re=40&ty=evtgen-mixed&dt=on_resonance&bl=caseB&dv=zfserv&st=5

#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/sig_kpi/mdst_setA/evtgen_exp_31__dtokpi_bp_sig__A_-0.mdst 5000
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/sig_kpi/mdst_setA/evtgen_exp_31__dtokpi_bp_sig__A_-0.mdst 5000
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/sig_kpi/mdst_setA/evtgen_exp_47__dtokpi_b0_sig__A_-6.mdst 1000
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/nor_kpi/mdst_setA/evtgen_exp_07__dtokpi_b0_nor__A_-0.mdst 2000
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/nor_kpi/mdst_setA/evtgen_exp_31__dtokpi_b0_nor__A_-1.mdst 1000

#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/nor_kpi_dpi0/mdst_setA/evtgen_exp_31__dtokpi_dpi0_b0_nor__A_-0.mdst 50000
#process_event /home/belle/syutaro/Store/dstrtaunu/mcprod/gsim/nor_kpi_dpi0/mdst_setA/evtgen_exp_31__dtokpi_dpi0_bp_nor__A_-0.mdst 10

#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-uds&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-charm&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-mixed&dt=on_resonance&bl=caseB&dv=zfserv&st=0
#process_url http://bweb3/montecarlo.php?ex=31&rs=231&re=231&ty=evtgen-charged&dt=on_resonance&bl=caseB&dv=zfserv&st=0

terminate
EOF

h2root test.hbk && rm test.hbk


