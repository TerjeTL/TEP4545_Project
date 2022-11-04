#!/bin/bash
cd "${0%/*}"

let file_name="force_coeffs_default.csv"

#### ARGUMENT HANDLING ####

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -n|--name)
      file_name="$2"
      shift # past argument
      shift # past value
      ;;
    --default)
      DEFAULT=YES
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

###########################


# Uncomment header and save as csv
sed -e 's/# T/T/g' ./postProcessing/forceCoeffs_object/0/coefficient.dat > ../../../Data/${file_name}

# Change dir
cd ../../../Data/

# Reformat csv with format: comment=#, sep="\t"
sed -i -e 's/ //g' ${file_name}
