#!/bin/bash
cd "${0%/*}"

#### DEFAULT CONFIG ####

num_cpu=8
clean_run=true
num_lt=5.0
num_cell_density="2.25e3"
num_dt="3e-4"
num_wall_distance=2.0

########################


#### DEFS ######

no_prompt=false

################


#### ARGUMENT HANDLING ####

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -y|-Y)
      no_prompt=true
      shift # past argument (no value)
      ;;
    -p|--processors)
      num_cpu=$2
      shift # past argument
      shift # past value
      ;;
    -dt|--time_step)
      num_dt=$2
      shift # past argument
      shift # past value
      ;;
    -lt)
      num_lt=$2
      shift
      shift
      ;;
    -dh|--cell_density)
      num_cell_density=$2
      shift
      shift
      ;;
    -wh|--wall_distance)
      num_wall_distance=$2
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

# Copy mesh file to run directory (nie - using blockmesh here)
#cp ./../../../../Meshing/flat_plate_in_channel.msh .

if [ "$clean_run" = true ]; then
  # Remove previous results FOOKIN LIES
  openfoam2206 foamListTimes -rm
  rm -r ./postProcessing
  
  # Change parameters
  sed -i "s/^LT_ratio.*$/LT_ratio ${num_lt};/" ./system/blockMeshDict
  sed -i "s/^CellRes.*$/CellRes ${num_cell_density};/" ./system/blockMeshDict
  sed -i "s/^DistancePlateWall.*$/DistancePlateWall ${num_wall_distance};/" ./system/blockMeshDict

  # Generate mesh
  openfoam2206 blockMesh

  # Check mesh
  openfoam2206 checkMesh
  if [ "$no_prompt" = false ]; then

    read -r -p "Continue? [y/N] " response
    if [[ "$response" =~ ^([yY][eE][sS]|[yY])$ ]]
    then
      continue
    else
      exit 0
    fi
  fi
  
  # Split mesh for parallelization
  let num_cpu_x=num_cpu/2
  sed -i "s/^numberOfSubdomains.*$/numberOfSubdomains\t${num_cpu};/" ./system/decomposeParDict
  sed -i "s/^\tn.*$/\tn (${num_cpu_x} 2 1);/" ./system/decomposeParDict

  openfoam2206 decomposePar -force

  sed -i "s/^deltaT.*$/deltaT\t${num_dt};/" ./system/controlDict
fi


# Run icoFoam in parallel
openfoam2206 mpirun --allow-run-as-root -np "$num_cpu" icoFoam -parallel
