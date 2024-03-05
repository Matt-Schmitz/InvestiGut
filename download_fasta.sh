#!/bin/bash

mkdir -p data


urls=(
    "https://zenodo.org/records/8398776/files/Investigut_Proteins.fa.xz_0?download=1"
    "https://zenodo.org/records/8404062/files/Investigut_Proteins.fa.xz_1?download=1"
)


download_url() {
    local url="$1"
    local output_file="data/$(basename "${url%%\?*}")"
    local retries=3 
    local delay=5

    for (( i=1; i<=$retries; i++ )); do
        echo "Downloading $url (attempt $i)..."
        if wget -q --show-progress -c "$url" -O "$output_file"; then
            echo "Download successful: $output_file"
            return 0
        else
            echo "Download failed: $output_file"
            echo "Retrying in $delay seconds..."
            sleep "$delay"
        fi
    done

    echo "Failed to download $url after $retries attempts."
    return 1
}

for url in "${urls[@]}"; do
    download_url "$url"
done

cat "data/Investigut_Proteins.fa.xz_0" "data/Investigut_Proteins.fa.xz_1" | xzcat | diamond makedb --in - -d data/investigut