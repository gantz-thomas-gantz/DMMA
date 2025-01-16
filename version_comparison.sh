#!/bin/bash

# Set the base URL
BASE_URL="https://ppar.tme-crypto.fr/jan.marxen/"

# Range of problem sizes
PROBLEM_SIZES=(25 26 27 28 29 30)  # You can modify this range as needed

# Versions to test
VERSIONS=("original" "version1" "version3" "version4" "version5")  # Define all versions here

# Loop through each problem size
nodes=1

# Prepare the weak_scaling.txt file to store the results
echo "Problem Size, Number of Nodes, ${VERSIONS[@]}" > versions_weak_scaling.txt

for problem_size in "${PROBLEM_SIZES[@]}"; do
    # Download the text file for the current problem size
    URL="${BASE_URL}${problem_size}"
    wget -q -O challenge.txt "$URL"

    # Extract parameters n, C0, and C1 from the text file
    n=$(grep -oP "(?<=--n )\d+" challenge.txt)

    # Extract C0 and C1 correctly (handle the two elements in each tuple)
    C0=$(grep -oP "(?<=C0 = \().*(?=\))" challenge.txt | tr -d '[:space:]')
    C1=$(grep -oP "(?<=C1 = \().*(?=\))" challenge.txt | tr -d '[:space:]')

    # Check if the values are extracted properly
    if [ -z "$n" ] || [ -z "$C0" ] || [ -z "$C1" ]; then
        echo "Failed to extract parameters for problem size $problem_size."
        continue
    fi

    # Split C0 and C1 into their respective elements (using comma as delimiter)
    C0_SECOND=$(echo $C0 | cut -d',' -f1)
    C0_FIRST=$(echo $C0 | cut -d',' -f2)
    C1_SECOND=$(echo $C1 | cut -d',' -f1)
    C1_FIRST=$(echo $C1 | cut -d',' -f2)

    # Print the extracted parameters for verification
    echo "Problem size: $problem_size"
    echo "n: $n, C0: ($C0_FIRST, $C0_SECOND), C1: ($C1_FIRST, $C1_SECOND)"

    # Create hostfile
    INPUT_FILE="$OAR_NODE_FILE"
    HOST_FILE="mpi_hostfile"

    # Deduplicate and count slots for each node
    awk '{slots[$1]++} END {for (host in slots) print host " slots=" slots[host]}' "$INPUT_FILE" > "$HOST_FILE"

    # Initialize a string to hold execution times for the current problem size
    EXEC_TIMES=""

    # Inner loop to run each version and measure execution time
    for VERSION in "${VERSIONS[@]}"; do
        # Run the MPI program and time it
        START_TIME=$(date +%s)
        mpiexec --hostfile "$HOST_FILE" -n "$nodes" build/$VERSION --n "$n" --C0 "$C0_FIRST$C0_SECOND" --C1 "$C1_FIRST$C1_SECOND"
        END_TIME=$(date +%s)

        # Calculate the execution time
        EXEC_TIME=$((END_TIME - START_TIME))
        echo "Execution time for problem size $problem_size, version $VERSION: $EXEC_TIME seconds."

        # Append the execution time to the EXEC_TIMES string
        EXEC_TIMES="$EXEC_TIMES, $EXEC_TIME"
    done

    # Log the result into weak_scaling.txt
    echo "$problem_size, $nodes$EXEC_TIMES" >> versions_weak_scaling.txt

    # Clean up the downloaded file
    rm challenge.txt

    # Double the number of nodes for the next iteration
    nodes=$((nodes * 2))
done

