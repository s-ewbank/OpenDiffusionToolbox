#!/bin/bash
#SBATCH -p normal
#SBATCH -c 1
#SBATCH --cpus-per-task=4
#SBATCH --time 48:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=snewbank@stanford.edu
#SBATCH --output=/scratch/users/snewbank/slurm_output//slurm-%j.out
#

config_path=""
script=""
dir=""

while [[ "$#" -gt 0 ]]; do
  case $1 in
    --config) config_path="$2"; shift ;;
    --script) script="$2"; shift ;;
    --dir) d="$2"; shift ;;
    *) echo "Unknown parameter: $1" >&2; exit 1 ;;
  esac
  shift
done

if [[ -z "$config_path" || -z "$script" ]]; then
  echo "Usage: $0 --config <path/to/config.txt> --script <path/to/script_to_submit_job_for.sh> "
  exit 1
fi

filedir="$(cd "$(dirname "$BASH_SOURCE")" && pwd)"
source ${filedir}/../READ_CONFIG.sh $config_path

singularity exec $container_path bash $script --config $config_path --dir $d