#! /bin/tcsh -f

if( $#argv != 2 )then
    echo " Usage : $0 [EVENT][STREAM]"
    echo '[EVENT ] = {uds, charm, mixed, charged}'
    echo '[STREAM] = {0-5(9)}'
    exit 1
endif

#****************************************************************************************************************
set DIR    = '~/ewp/skim/gmc/right/index/'
set MDL    = 'DSTRTAUNU'
set TAG    = ''    # _NAME
set FLAG   = 'gMC'
#set TYPE   = 'on_resonance' # [on_resonance, continuum, 5S_scan, 5S_onresonance, 3S_scan, 2S_scan, 1S_scan]
set EVENT  = $1 # [uds, charm, mixed, charmed]
set STREAM = $2 # [ 0-5 ]
set list_data = `awk '{print $1}'  ~syutaro/script/data_list_on_resonance.txt`  
#****************************************************************************************************************
set RUN   = 'run_'$MDL'_'$FLAG'_'$EVENT'_'$STREAM$TAG'.sh'
#****************************************************************************************************************

rm -f $RUN

# $1 : $DIR                                 # location of mdst file
# $2 : $SET($MDL/$TAG)               # module name and etc...
# $3 : $f                                   # mdst file name for sigMC
# $4 : $DATA($FLAG)                         # MC flag
# $5 : $RUN                                 # run script file name

#  ************************ MC ***********************
echo '[ ' $MDL ' ' $FLAG ' ' $TAG ' ]' 

set CNT=0
set POST=0
		   
foreach f($list_data)
     set expNo  = `echo $f | cut -d "/" -f 1`
     if( $expNo != $POST ) then
       @ CNT = 0
     endif
     @ POST = $expNo
     @ CNT += 1
     echo -n $CNT ' '
     ./mk_script_gMC.sh $MDL/$EVENT/$STREAM/$TAG $f $FLAG $RUN $CNT
  end
end

chmod 744 $RUN
