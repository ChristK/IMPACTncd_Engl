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
# 5. Grants supplementary groups (if ADDITIONAL_GIDS is provided):
#    - ADDITIONAL_GIDS is a space-separated list of numeric host GIDs
#      (the setup scripts pass `id -G`)
#    - Each GID is registered in the container (existing GIDs are reused;
#      missing ones get a placeholder group) and added to the user
#    - This is what makes group-based permissions on shared bind mounts work
#      (e.g. the team's shared synthpop cache): `gosu uid:gid` alone strips
#      supplementary groups
#
# 6. Executes the command as the specified user:
#    - Uses `gosu "$USER_ID"` (no :gid) so the user's FULL group list from
#      /etc/group applies — primary plus the supplementary groups from step 5
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
# starts. The setup_*_docker_env.sh scripts pass the necessary environment variables:
#   -e USER_ID="$(id -u)"
#   -e GROUP_ID="$(id -g)"
#   -e USER_NAME="$(whoami)"
#   -e GROUP_NAME="$(id -gn)"
#   -e ADDITIONAL_GIDS="$(id -G)"  # optional: supplementary host GIDs (step 5);
#                                  # omit for the old single-group behavior
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

# Resolve the container-side username for the target UID (needed below; the
# user is guaranteed to exist at this point — created above or pre-existing).
TARGET_USER_NAME=$(getent passwd "$USER_ID" | cut -d: -f1)

# If the (pre-existing) user's primary group differs from the requested
# GROUP_ID, align it — `gosu "$USER_ID"` below uses the passwd entry's primary
# group, whereas the previous `gosu uid:gid` form forced GROUP_ID.
if [ -n "$TARGET_USER_NAME" ]; then
    CURRENT_GID=$(getent passwd "$USER_ID" | cut -d: -f4)
    if [ "$CURRENT_GID" != "$GROUP_ID" ]; then
        echo "Aligning primary group of $TARGET_USER_NAME to GID $GROUP_ID"
        usermod -g "$GROUP_ID" "$TARGET_USER_NAME"
    fi
fi

# -----------------------------------------------------------------------------
# Supplementary groups (ADDITIONAL_GIDS: space-separated numeric host GIDs,
# passed by the setup scripts as `id -G`).
#
# WHY: `gosu uid:gid` grants ONLY that single group — supplementary groups are
# stripped. That silently breaks group-based permissions on bind mounts (e.g.
# a shared synthpop cache directory owned by a shared host group): the host
# user is in the group, but the container process is not. So we register each
# host GID on the container user and exec `gosu "$USER_ID"` (no :gid) instead,
# which applies the user's full group list from /etc/group (initgroups).
#
# Only numeric GIDs are accepted (defense against a malformed/hostile value —
# this block runs as root). Unset/empty ADDITIONAL_GIDS keeps the previous
# single-group behavior, so older setup scripts remain fully compatible.
# -----------------------------------------------------------------------------
if [ -n "${ADDITIONAL_GIDS:-}" ] && [ -n "$TARGET_USER_NAME" ]; then
    for gid in $ADDITIONAL_GIDS; do
        case "$gid" in
            ''|*[!0-9]*)
                echo "Skipping invalid supplementary GID '$gid' (numeric GIDs only)"
                continue
                ;;
        esac
        # Never grant the root group — `id -G` can legitimately contain 0 on
        # some admin hosts, but group 0 in the container is not something a
        # shared-cache mount should ever need. Numeric -eq also catches "00".
        if [ "$gid" -eq 0 ]; then
            echo "Skipping GID 0 (root group is never granted)"
            continue
        fi
        # Primary group is already set; skip it.
        [ "$gid" = "$GROUP_ID" ] && continue
        # Create a placeholder group if no container group has this GID yet
        # (the kernel checks numeric GIDs; the name is irrelevant).
        if ! getent group "$gid" > /dev/null 2>&1; then
            groupadd -g "$gid" "hostgrp_$gid"
        fi
        usermod -aG "$gid" "$TARGET_USER_NAME"
    done
    echo "Supplementary groups for $TARGET_USER_NAME: $(id -G "$TARGET_USER_NAME")"
fi

# Execute the command as the specified user. `gosu "$USER_ID"` (without :gid)
# resolves the user's passwd entry and applies primary + supplementary groups.
# GUARD: if user creation failed above (e.g. the sanitized host username
# collides with a system account baked into the image, or groupadd rejected an
# exotic host group name), there is NO passwd entry for USER_ID — and
# `gosu <uid>` for a passwd-less UID falls back to gid=0 (ROOT group). In that
# case degrade to the old explicit uid:gid form: correct primary group, no
# supplementary groups (matching the pre-ADDITIONAL_GIDS behavior).
if [ -n "$TARGET_USER_NAME" ]; then
    echo "Executing command as user $USER_ID (primary GID $GROUP_ID)"
    exec gosu "$USER_ID" "$@"
else
    echo "WARNING: no passwd entry for UID $USER_ID (user creation failed above);" >&2
    echo "         supplementary groups unavailable — running as $USER_ID:$GROUP_ID" >&2
    exec gosu "$USER_ID:$GROUP_ID" "$@"
fi
