#!/bin/bash

# Generate deduplicated host file
INPUT_FILE="$OAR_NODE_FILE"
HOST_FILE="mpi_hostfile"

# Deduplicate and count slots for each node
awk '{slots[$1]++} END {for (host in slots) print host " slots=" slots[host]}' "$INPUT_FILE" > "$HOST_FILE"

# Overwrite the output file to start fresh
OUTPUT_FILE="stress_test_output.txt"
> "$OUTPUT_FILE"  # This ensures the file is empty at the start of each run

# Run the MPI application in the background
#mpiexec --hostfile "$HOST_FILE" --verbose -n 8 ./v4 --n 31 --C0 0e3b626cd41b8a62 --C1 e3637e86c789bc99
mpiexec -n 32 ./v4 --n 31 --C0 0e3b626cd41b8a62 --C1 e3637e86c789bc99

#mpiexec --hostfile "$HOST_FILE" --verbose -n 4 valgrind --tool=massif --massif-out-file=massif.out.%p ./v4 --n 28 --C0 868246f0f3ec2059 --C1 61e6c7d55f0881ea
#mpiexec --hostfile "$HOST_FILE" --verbose -n 1 valgrind --tool=massif --time-unit=ms --detailed-freq=10 --massif-out-file=massif.out.%p ./v4 --n 32 --C0 4d0d80193cc34050 --C1 8aa1e483e1dc62e6

#mpiexec --hostfile "$HOST_FILE" --verbose -n 2 valgrind --tool=massif --time-unit=ms --detailed-freq=5 --massif-out-file=massif.out.%p ./v4 --n 28 --C0 868246f0f3ec2059 --C1 61e6c7d55f0881ea




