#!/bin/bash

# Cross-platform package updater for macOS and Linux
# Updates apt-packages.txt with actual installed versions from Docker build output

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APT_PACKAGES_FILE="$SCRIPT_DIR/apt-packages.txt"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color

usage() {
    echo -e "${MAGENTA}APT Packages Version Updater${NC}"
    echo -e "${MAGENTA}=============================${NC}"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -f, --file FILE    Parse Docker build log from file"
    echo "  -i, --interactive  Run in interactive mode"
    echo "  -h, --help        Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                          # Parse from stdin"
    echo "  $0 -f build.log            # Parse from file"
    echo "  $0 -i                       # Interactive mode"
    echo ""
    echo "  # Capture build output and update packages:"
    echo "  docker build ... 2>&1 | tee build.log | $0"
}

parse_docker_output() {
    local input_source="$1"
    local updates=()
    local in_update_section=false
    
    while IFS= read -r line; do
        if [[ "$line" =~ "PACKAGE VERSION UPDATES DETECTED:" ]]; then
            in_update_section=true
            continue
        fi
        
        if $in_update_section; then
            if [[ "$line" =~ ^=+ ]]; then
                if [[ "$line" =~ "RECOMMENDED ACTION:" ]]; then
                    continue
                else
                    in_update_section=false
                    continue
                fi
            fi
            
            if [[ "$line" =~ ^([a-zA-Z0-9\.-]+)=(.+)$ ]]; then
                local package_name="${BASH_REMATCH[1]}"
                local version="${BASH_REMATCH[2]}"
                updates+=("$package_name=$version")
                echo -e "${GREEN}Found update: $package_name = $version${NC}" >&2
            fi
        fi
    done < "$input_source"
    
    printf '%s\n' "${updates[@]}"
}

update_apt_packages_file() {
    local -a updates=("$@")
    
    if [[ ! -f "$APT_PACKAGES_FILE" ]]; then
        echo -e "${RED}Error: apt-packages.txt not found at: $APT_PACKAGES_FILE${NC}" >&2
        exit 1
    fi
    
    if [[ ${#updates[@]} -eq 0 ]]; then
        echo -e "${BLUE}No package version updates found in the input.${NC}"
        return 0
    fi
    
    # Create associative array for updates
    declare -A update_map
    for update in "${updates[@]}"; do
        local package_name="${update%%=*}"
        local version="${update#*=}"
        update_map["$package_name"]="$version"
    done
    
    # Create backup
    local backup_file="$APT_PACKAGES_FILE.backup.$(date +%Y%m%d_%H%M%S)"
    cp "$APT_PACKAGES_FILE" "$backup_file"
    echo -e "${CYAN}Backup created: $backup_file${NC}"
    
    # Update file
    local temp_file=$(mktemp)
    local updated=false
    
    while IFS= read -r line; do
        if [[ "$line" =~ ^([a-zA-Z0-9\.-]+)= ]]; then
            local package_name="${BASH_REMATCH[1]}"
            if [[ -n "${update_map[$package_name]}" ]]; then
                local new_version="${update_map[$package_name]}"
                echo "$package_name=$new_version" >> "$temp_file"
                echo -e "${YELLOW}Updating $package_name to version $new_version${NC}"
                updated=true
            else
                echo "$line" >> "$temp_file"
            fi
        else
            echo "$line" >> "$temp_file"
        fi
    done < "$APT_PACKAGES_FILE"
    
    if $updated; then
        mv "$temp_file" "$APT_PACKAGES_FILE"
        echo -e "${GREEN}Updated $APT_PACKAGES_FILE successfully!${NC}"
    else
        rm "$temp_file"
        echo -e "${BLUE}No updates were needed.${NC}"
    fi
}

# Main execution
main() {
    local build_log_file=""
    local interactive=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -f|--file)
                build_log_file="$2"
                shift 2
                ;;
            -i|--interactive)
                interactive=true
                shift
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}" >&2
                usage >&2
                exit 1
                ;;
        esac
    done
    
    echo -e "${MAGENTA}APT Packages Version Updater${NC}"
    echo -e "${MAGENTA}=============================${NC}"
    
    local input_source
    local cleanup_temp=false
    
    if [[ -n "$build_log_file" ]]; then
        if [[ ! -f "$build_log_file" ]]; then
            echo -e "${RED}Error: Build log file not found: $build_log_file${NC}" >&2
            exit 1
        fi
        input_source="$build_log_file"
        echo -e "${CYAN}Reading from build log file: $build_log_file${NC}"
    elif $interactive; then
        echo -e "${CYAN}Enter Docker build output (Ctrl+D when finished):${NC}"
        input_source=$(mktemp)
        cat > "$input_source"
        cleanup_temp=true
    else
        echo -e "${CYAN}Reading from stdin...${NC}"
        input_source=$(mktemp)
        cat > "$input_source"
        cleanup_temp=true
    fi
    
    # Parse updates
    local -a updates
    readarray -t updates < <(parse_docker_output "$input_source")
    
    # Cleanup temp file if created
    if $cleanup_temp; then
        rm -f "$input_source"
    fi
    
    if [[ ${#updates[@]} -eq 0 ]]; then
        echo -e "${YELLOW}No package version updates found in the build output.${NC}"
        exit 0
    fi
    
    echo -e "\n${CYAN}Found ${#updates[@]} package updates:${NC}"
    for update in "${updates[@]}"; do
        echo "  $update"
    done
    
    if $interactive; then
        echo -e "\n${CYAN}Review the updates above.${NC}"
        read -p "Do you want to apply these updates to apt-packages.txt? (y/N): " confirm
        if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
            echo -e "${YELLOW}Updates cancelled by user.${NC}"
            exit 0
        fi
    fi
    
    update_apt_packages_file "${updates[@]}"
}

# Run main function with all arguments
main "$@"
