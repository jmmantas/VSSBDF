#!/bin/bash

# Clean up and compile the program nice
make clean
make

# Check if compilation was successful
if [ $? -ne 0 ]; then
    echo "Compilation failed. Exiting."
    exit 1
fi

# Define the dt values and other parameters
#Tol_values=(1e-3 1e-4 1e-5  1e-6 1e-7  )
splitting_types=(0 1)
#problem_ids=(0 1 2 3 4)
problem_ids=(3)
order_values=( 1 2 3 4)
#order_values=( 1)
N=210
dt=0.001
dtRK=0.0001
tf=10
NumExp=5
Tol=1e-2  # This is 10^-6 in bash syntax

declare -A problem_names=(
    [0]="1D Advection-Diffusion"
    [1]="Combustion"
    [2]="1D Brusselator"
    [3]="2D Heat"
    [4]="CUSP"
    [5]="SimpleAD1"
)

declare -A splitting_names=(
    [0]="physics"
    [1]="Jacobian"
)



if [ ! -d "results" ]; then
    mkdir results
fi

# Loop over problem_id values
for problem_id in "${problem_ids[@]}"
do
    subfolder="results/${problem_names[$problem_id]}"
    if [ ! -d "$subfolder" ]; then
        mkdir "$subfolder"
    fi

    for order in "${order_values[@]}"
    do
        #for splittingtype in "${splitting_types[@]}"
        #do
           # output_filename="$subfolder/combined_output_for_${problem_names[$problem_id]}_${splitting_names[$splittingtype]}_Order_$order.txt"
            output_filename="$subfolder/${problem_names[$problem_id]}_$order.txt"
            # Optionally clear the file before writing to it (especially useful if the script is re-run)
            > "$output_filename"
            
            #for Tol in "${Tol_values[@]}"
           # do
                # Outer loop to run everything three times
                for run in {1..3}
                do
               
                   # echo "==== Starting Run $run ====" >> "$output_filename"
                  #  echo "==== order= $order, N=$N, splittingType=$splittingtype, dt=$dt, Final Time $tf,  Tol=$Tol ====" >> "$output_filename"
                    # Appending to the file using >>
                    ./IVP_Solver $problem_id $order $N $tf $dtRK $dt $NumExp $Tol  >> "$output_filename" 2>&1
                done
            #done
       # done
    done
done
