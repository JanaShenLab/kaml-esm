#!/bin/bash
#
#                         ,--,
#           ,-.----.   ,---.'|             ____  ,-.----.
#           \::::/  \  |:::| :           ,'  , `.\::::/  \    ,----..
#           |::::    \ ::::: |        ,-+-,.' _ ||::::    \  /:::/   \
#           |:::|  .\ :|:::' :     ,-+-. ;   ,:|||:::|  .\ :|::::     :
#           .::::  |: |;:::; '    ,--.'|'   |::;|.::::  |: |.:::|  ;. /
#           |:::|   \ :':::| |__ |:::|  ,', |::':|:::|   \ :.:::; /--`
#           |:::: .   /|:::| :.'||:::| /::| |::|||:::: .   /;:::| ;  __
#           ;:::| |`-' '::::    ;':::| :::| :::|,;:::| |`-' |:::: |.' .'
#           |:::| ;    |:::|  ./ ;:::. |::; |--' |:::| ;    .:::| '_.' :
#           ::::' |    ;:::: ;   |:::: |::| ,    ::::' |    ':::; ::\  |
#           ::::: :    |:::,/    |:::: '::|/     ::::: :    ':::| '/  .'
#           |:::| :    '---'     ;:::| |`-'      |:::| :    |::::    /
#           `---'.|              |:::;/          `---'.|     \:::\ .'
#             `---`              '---'             `---`      `---`
#
#                      Protein Language Model Playground
#
#                              Guy W. Dayhoff II
#            Shchool of Pharmacy, University of Maryland Baltimore
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Check for required arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <conda|venv>"
    exit 1
fi

TOOL="$1"
PYTHON_V=3.10.12
SETUP_DIR="$(pwd)"
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../bin" && pwd)"

# Function to modify venv activation script
modify_venv_activation() {
    local env_name="$1"
    local activate_script="$env_name/bin/activate"

    echo "Modifying venv activation script: $activate_script"

    # Add environment variables & script path
    {
        echo -e "\n# Added by setup_envs.sh"
        echo "export PATH=\"$SCRIPTS_DIR:\$PATH\""
        #echo "export PLMPG_FORGE_TOKEN=\"$SETUP_DIR/forge_token.txt\""
    } >> "$activate_script"

    # Add alias masking based on env_name
    if [ "$env_name" == "plmpg-esm2" ]; then
        {
            echo ""
            echo "# Alias masking for plmpg-esm2"
            echo "alias plmpg-esm3-forge-fold='echo ERROR: esm3 scripts"
            echo "    only available in plmpg-esm3 environment'"
            echo "alias plmpg-esm3-local-fold='echo ERROR: esm3 scripts"
            echo "    only available in plmpg-esm3 environment'"
            echo "alias plmpg-esmC-extractor='echo ERROR: esmC extractor"
            echo "    only available in plmpg-esm3 environment'"
        } >> "$activate_script"
    elif [ "$env_name" == "plmpg-esm3" ]; then
        {
            echo ""
            echo "# Alias masking for plmpg-esm3"
            echo "alias plmpg-esm2-extractor='echo ERROR: esm2 extractor"
            echo "    only available in plmpg-esm2 environment'"
        } >> "$activate_script"
    fi

    # Cleanup on deactivate: unalias commands and call original deactivate.
    {
        echo ""
        echo "# Cleanup: Unalias commands when deactivating"
        echo "unalias plmpg-esm3-forge-fold 2>/dev/null"
        echo "unalias plmpg-esm3-local-fold 2>/dev/null"
        echo "unalias plmpg-esmC-extractor 2>/dev/null"
        echo "unalias plmpg-esm2-extractor 2>/dev/null"
        echo ""
        echo "# Ensure original deactivate() is called after cleanup"
        echo "if [[ \"\$(type -t deactivate)\" == \"function\" ]]; then"
        echo "    # Rename original deactivate to _original_deactivate"
        echo "    eval \"\$(declare -f deactivate | sed 's|^deactivate[[:space:]]*()|_original_deactivate()|')\""
        echo "    function deactivate() {"
        echo "        unalias plmpg-esm3-forge-fold 2>/dev/null"
        echo "        unalias plmpg-esm3-local-fold 2>/dev/null"
        echo "        unalias plmpg-esmC-extractor 2>/dev/null"
        echo "        unalias plmpg-esm2-extractor 2>/dev/null"
        echo "        _original_deactivate"
        echo "    }"
        echo "fi"
    } >> "$activate_script"
}

# Function to create and set up a Conda environment
setup_conda_env() {
    local env_name="$1"
    local requirements_file="$2"

    echo "Setting up Conda environment: $env_name"
    conda create --yes --name "$env_name" python=${PYTHON_V} || {
        echo "Error creating Conda env"
        exit 1
    }

    conda_env_path=$(conda env list | grep "^$env_name\s" | awk '{print $2}')
    [[ -z "$conda_env_path" ]] && {
        echo "Could not locate Conda env path"
        exit 1
    }

    echo "Installing dependencies in $env_name"
    conda run --name "$env_name" pip install --upgrade pip
    conda run --name "$env_name" pip install -r "$requirements_file"

    echo "Conda environment '$env_name' setup complete."
}

# Function to create and set up a Python venv
setup_venv_env() {
    local env_name="$1"
    local requirements_file="$2"

    echo "Setting up Python venv: $env_name"
    python3 -m venv "$env_name" || {
        echo "Error creating venv"
        exit 1
    }

    modify_venv_activation "$env_name"

    echo "Installing dependencies in venv: $env_name"
    source "$env_name/bin/activate"
    pip install --upgrade pip
    pip install -r "$requirements_file"
    deactivate

    echo "Venv environment '$env_name' setup complete."
}

# Main script logic: loop over *_env.txt files
for file in *_env.txt; do
    env_name="${file%_env.txt}"
    case "$TOOL" in
        conda)
            setup_conda_env "$env_name" "$file"
            ;;
        venv)
            setup_venv_env "$env_name" "$file"
            ;;
        *)
            echo "Unsupported tool: $TOOL"
            exit 1
            ;;
    esac
done

echo "All environments set up successfully."

