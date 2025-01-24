#!/bin/bash

# Check if a config file was provided
if [[ $# -lt 1 ]]; then
  echo "Usage: bash read.sh <config_path>"
  exit 1
fi

# Read the config file path from the argument
config_path="$1"

# Ensure the file exists
if [[ ! -f "$config_path" ]]; then
  echo "Error: Config file '$config_path' not found!"
  exit 1
fi

# Parse the config file
while IFS='=' read -r key value; do
  # Remove whitespace around key and value
  key=$(echo "$key" | xargs)
  value=$(echo "$value" | xargs)

  # Skip empty lines and comments
  if [[ -n "$key" && ! "$key" =~ ^# ]]; then
    # Export the variable dynamically
    export "$key"="$value"
  fi
done < "$config_path"

# Handle arrays separately
if [[ -n "$volumes" ]]; then
  IFS=',' read -r -a volumes <<< "$volumes"
fi

if [[ -n "$special_ids" ]]; then
  IFS=',' read -r -a special_ids <<< "$special_ids"
fi

bash ${scripts_dir}/ODTM-misc_startup.sh

# Output the variables (optional for debugging)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Loaded variables:"
echo "  - rootdir: $rootdir"
echo "  - groups_file: $groups_file"
echo "  - scripts_dir: $scripts_dir"
echo "  - container_path: $container_path"
echo "  - mailtype: $mailtype"
echo "  - mailuser: $mailuser"
echo "  - slurm_out_path: $slurm_out_path"
echo "  - organism: $organism"
echo "  - volumes: ${volumes[@]}"
echo "  - do_mppca: $do_mppca"
echo "  - do_slice: $do_slice"
echo "  - special_ids: ${special_ids[@]}"
