# Re-run one treatment chamber to correct mistakes, add output variables, etc.
ens_id=01440

exeroot=${E3SM_ROOT}/output/UQ_20240104_US-SPR_ICB1850CNRDCTCBC_ad_spinup/bld

cd ${E3SM_ROOT}/output/UQ/UQ_20240104_US-SPR_ICB20TRCNPRDCTCBC/g${ens_id}

treatments=('TAMB' 'T0.00' 'T2.25' 'T4.50' 'T6.75' 'T9.00' \
            'T0.00CO2' 'T2.25CO2' 'T4.50CO2' 'T6.75CO2' 'T9.00CO2')
plots=(07 06 20 13 8 17 19 11 4 16 10)

# Get the number of elements
length=${#treatments[@]}

plot_prev=10

# Iterate over the indices
for ((i=0; i<length; i++)); do
    trmt="${treatments[i]}"
    plot="${plots[i]}"

    sed -i s/plot${plot_prev}/plot${plot}/g lnd_in

    ${exeroot}/e3sm.exe >& e3sm_log_${trmt}.txt

    mv *.h1.*.nc ./${trmt}
    mv *.h2.*.nc ./${trmt}

    plot_prev=${plot}
done