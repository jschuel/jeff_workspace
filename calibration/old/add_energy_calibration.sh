if [ $# -lt 1 ]
then
    echo "./update_ntuple_all.sh <run number>"
    echo "  OR"
    echo "./update_ntuple_all.sh <first run number> <last run number>"
    echo ""
    exit 1
fi

first_run=$1
last_run=$first_run

if [ $# == 2 ]
then
    last_run=$2
fi

iRun=$first_run

while [ $iRun -le $last_run ]
do

    root -l -q -b 'fast_energy_calibration.C('\"iiwi\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"honu\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"kohola\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"nene\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"tako\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"humu\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"elepaio\"', '$iRun')'
    root -l -q -b 'fast_energy_calibration.C('\"palila\"', '$iRun')'                                                                                    

    let iRun=$iRun+1
done
