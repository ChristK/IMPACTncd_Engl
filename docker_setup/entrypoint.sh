#!/bin/bash
# -----------------------------------------------------------------------------
# entrypoint.sh - Docker Container User Identity Management
# -----------------------------------------------------------------------------
#
# PURPOSE:
# This entrypoint script handles dynamic user creation for non-root Docker execution.
# It allows the container to run with the same user ID and group ID as the host user,
# which is crucial for avoiding file permission issues when mounting directories
# between the host and container.
#
# HOW IT WORKS:
# 1. Gets user information from environment variables:
#    - USER_ID and GROUP_ID: the numeric IDs from the host
#    - USER_NAME and GROUP_NAME: the actual username and group name from the host
#    - Falls back to defaults (dockeruser/dockergroup) if not provided
#
# 2. Handles root execution:
#    - If no USER_ID is specified or if running as root (UID 0), executes command directly
#
# 3. Creates group if needed:
#    - Checks if a group with the specified GROUP_ID exists
#    - If not, creates a new group with the name GROUP_NAME and ID GROUP_ID
#
# 4. Creates user if needed:
#    - Checks if a user with the specified USER_ID exists
#    - If not, creates a new user with:
#      * Username: USER_NAME
#      * UID: USER_ID
#      * Primary group: GROUP_ID
#      * Home directory: /home/$USER_NAME
#      * Shell: /bin/bash
#    - Ensures the home directory has correct ownership
#
# 5. Executes the command as the specified user:
#    - Uses 'gosu' to switch to the non-root user and execute the original command
#
# WHY THIS MATTERS:
# Without this script, the container would run as root or a generic user, causing
# files created inside the container to have incorrect ownership when viewed from
# the host system. This script ensures that files created in the container have
# the same ownership as if they were created directly on the host by the actual user.
#
# This is particularly important for this project since directories are mounted
# between the host and container, and output files should be owned by the actual
# user account rather than root or a generic dockeruser.
#
# USAGE:
# This script is automatically called as the Docker ENTRYPOINT when the container
# starts. The setup_dev_docker_env.sh script passes the necessary environment variables:
#   -e USER_ID="$(id -u)"
#   -e GROUP_ID="$(id -g)"
#   -e USER_NAME="$(whoami)"
#   -e GROUP_NAME="$(id -gn)"
#
# -----------------------------------------------------------------------------

# Get the user and group IDs from environment variables passed by Docker
USER_ID=${USER_ID:-$(id -u)}
GROUP_ID=${GROUP_ID:-$(id -g)}
USER_NAME=${USER_NAME:-"dockeruser"}
GROUP_NAME=${GROUP_NAME:-"dockergroup"}

echo "Entrypoint debug: USER_ID=$USER_ID, GROUP_ID=$GROUP_ID, USER_NAME=$USER_NAME, GROUP_NAME=$GROUP_NAME"

# If no specific user is requested or running as root, just execute the command
if [ -z "$USER_ID" ] || [ "$USER_ID" = "0" ]; then
    echo "Running as root or no USER_ID specified"
    exec "$@"
fi

# Create group if it doesn't exist
if ! getent group "$GROUP_ID" > /dev/null 2>&1; then
    echo "Creating group with ID $GROUP_ID"
    # If group name already exists but with different ID, use a safe fallback
    if getent group "$GROUP_NAME" > /dev/null 2>&1; then
        echo "Group name $GROUP_NAME already exists, using fallback"
        GROUP_NAME="dockergroup_${GROUP_ID}"
    fi
    groupadd -g "$GROUP_ID" "$GROUP_NAME"
    echo "Created group: $GROUP_NAME with ID $GROUP_ID"
else
    echo "Group with ID $GROUP_ID already exists"
fi

# Create user if it doesn't exist
if ! getent passwd "$USER_ID" > /dev/null 2>&1; then
    echo "Creating user with ID $USER_ID"
    # Sanitize username (replace spaces and special characters with underscores)
    SAFE_USER_NAME=$(echo "$USER_NAME" | sed 's/[^a-zA-Z0-9]/_/g' | sed 's/__*/_/g' | sed 's/^_\|_$//g')
    
    # Ensure we have a valid username
    if [ -z "$SAFE_USER_NAME" ]; then
        SAFE_USER_NAME="dockeruser"
    fi
    
    echo "Sanitized username: $SAFE_USER_NAME"
    
    # Get the actual group name that was created (in case it was modified above)
    ACTUAL_GROUP_NAME=$(getent group "$GROUP_ID" | cut -d: -f1)
    echo "Using group name: $ACTUAL_GROUP_NAME"
    
    # Create user with home directory, using the actual group name
    useradd -u "$USER_ID" -g "$ACTUAL_GROUP_NAME" -m -s /bin/bash "$SAFE_USER_NAME"
    
    # Ensure home directory has correct ownership
    chown "$USER_ID:$GROUP_ID" "/home/$SAFE_USER_NAME"
    echo "Created user: $SAFE_USER_NAME with home /home/$SAFE_USER_NAME"
else
    echo "User with ID $USER_ID already exists"
fi

# Execute the command as the specified user
echo "Executing command as user $USER_ID:$GROUP_ID"
exec gosu "$USER_ID:$GROUP_ID" "$@"
